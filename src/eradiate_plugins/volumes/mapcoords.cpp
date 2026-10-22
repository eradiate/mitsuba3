#include <mitsuba/core/properties.h>
#include <mitsuba/core/transform.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/eradiate/extremum.h>


#include <drjit/util.h>

#include <set>
#include <type_traits>
#include <vector>


NAMESPACE_BEGIN(mitsuba)

/*
 * MapCoords
 * 
 * This file defines a static library to write and chain a set of coordinates
 * mappings. A mapping is defined as an object able to perform a transformation
 * of a point in a 3D coord system to another through an "apply" method. This library
 * provides builtin mappings and tools to combine them or maintain compatitibility
 * with other Mitsuba interfaces.
 *
 */

// Helper namespace to perform tests on the mapping functions
// We implement a concept through static SFINAE, until the codebase is ported
// to C++20.
namespace mapchecker {

    // Defines the return type of a call to the .apply method of the Mapping
    // templated type. The .apply method must accept a single parameter of type
    // Mapping::Point3f.
    // The method must be const, and the parameter is passed by const reference
    template<typename Mapping>
    using mapping_apply_expr_t = decltype(std::declval<const Mapping&>().apply(std::declval<const typename Mapping::Point3f&>()));


    // The method .apply is not present or has incorrect parameters
    template<typename Mapping, typename = void>
    struct inspect_apply_p3f { 
    private:
        struct checker_fallback_t { }; // Dummy fallback type
    public:
        using type = checker_fallback_t;
    };

    // The method .apply is present with the correct parameters.
    // If the compiler can form this struct from the template Mapping, then the
    // .apply method as defined by mapping_apply_expr_t exists and its return
    // type is stored in the ::type using.
    template<typename Mapping>
    struct inspect_apply_p3f<Mapping, std::void_t<mapping_apply_expr_t<Mapping>>> { 
        using type = mapping_apply_expr_t<Mapping>;
    };

    // Final check expression
    template<typename Mapping>
    constexpr bool has_apply_p3f_v =
        std::is_same_v<typename inspect_apply_p3f<Mapping>::type, typename Mapping::Point3f>;

    // Assert helper
    template<typename Mapping>
    static constexpr void assert_is_mapping() {
        static_assert(has_apply_p3f_v<Mapping>,
            "A Mapping type must possess an endomorphic Point3f -> Point3f "
            ".apply method. The signature must be: "
            "Mapping::Point3f apply(const Mapping::Point3f&) const"
        );
        static_assert(std::is_constructible_v<Mapping, const Properties&>,
            "A Mapping type must be constructible from scenee properties of "
            "type: const Properties&"
        );
    }

    // Optional .to_string() support: a Mapping is not required to define
    // one, so MappingChain::to_string() can still describe its mappings
    // when some (or all) of them don't.
    template<typename Mapping, typename = void>
    struct has_to_string : std::false_type { };

    template<typename Mapping>
    struct has_to_string<Mapping, std::void_t<decltype(std::declval<const Mapping&>().to_string())>>
        : std::true_type { };

    template<typename Mapping>
    constexpr bool has_to_string_v = has_to_string<Mapping>::value;

    template<typename Mapping>
    std::string mapping_to_string(const Mapping &m) {
        if constexpr (has_to_string_v<Mapping>)
            return m.to_string();
        else
            return "<mapping with no to_string()>";
    }

} // end namespace mapchecker


template <typename Float, typename Spectrum, typename Mapping>
struct mapping_traits {

    /*
     * Mapping traits
     *
     * Define the mapping specific behaviors
     *
     * This traits struct currently exposes default min and max functions for
     * which some mappings may require a specific implementation. min returns
     * inf, max returns -inf, thus the default behavior of a mapping on the
     * mapped volume is a noop.
     *
     * Extremum behavior can be adapted to a specific behavior as well.
     */

    MI_IMPORT_TYPES();

    static inline ScalarFloat max(const Mapping&) { return -dr::Infinity<ScalarFloat>; }
    static inline ScalarFloat min(const Mapping&) { return  dr::Infinity<ScalarFloat>; }
    static inline std::pair<Float, Float> extremum(
        const Mapping&,
        BoundingBox3f,
        Float min, Float max
    ) {
        return { min, max };
    }
};


template <typename Float, typename Spectrum, typename Backend>
struct backend_traits {

    /*
     * Data storage backend traits
     *
     * Defines the interface to a storage type for the actual volume data. By
     * default the backend type is assumed to follow the 
     * mitsuba::Volume<Float, Spectrum> interface.
     *
     * This traits effectively provides support to arbitrary storage types at
     * the minimal cost of writing a case-specific specialization of the
     * traits, without extra virtual or dynamic interface.
     *
     */

    MI_IMPORT_TYPES();

    static inline UnpolarizedSpectrum eval(const Backend &backend, const Interaction3f &it, Mask active) {
        return backend.eval(it, active);
    }

    static inline Float eval_1(const Backend &backend, const Interaction3f &it, Mask active) {
        return backend.eval_1(it, active);
    }

    static std::pair<Float, Float> extremum(const Backend &backend, BoundingBox3f bbox, Mask local) {
        return backend.extremum(bbox, local);
    }

    static ScalarVector3i resolution(const Backend &backend) { return backend.resolution(); }

    static ScalarFloat max(const Backend &backend) { return backend.max(); }

    static ScalarFloat min(const Backend &backend) { return backend.min(); }
};


template <typename Float, typename Spectrum, typename VolumeBackend, typename... MapTypes>
class MappingChain final : public Volume<Float, Spectrum> {
public:

    /*
     * Mapping chain template class
     *
     * Structure a set of Mappings from world coordinates to the volume backend
     * coordinates. Accepts an arbitrary number of mappings through a variadic
     * MapTypes template parameter.
     *
     * Mappings are constructed one after each other in the order of occurence
     * in the MapTypes variadic template parameter. They are passed the 
     * properties of the current constructor to maintain a compact dict
     * interface.
     *
     * Mappings are applied in chain, in the order of occurence in the MapTypes 
     * variadic template parameter.
     *
     * The MappingChain performs a transformation from world coordinates
     * inverting the Volume<Float, Spectrum> interface's "to_world" attribute
     * before applying the chain of mappings.
     *
     * Extremum properties are defined using the target backend coordinates.
     */

    MI_IMPORT_BASE(Volume, m_to_local, m_bbox)
    MI_IMPORT_TYPES(VolumeGrid, ExtremumStructure)

    using BackendDispatch = backend_traits<Float, Spectrum, VolumeBackend>;

    template <typename M>
    using MappingDispatch = mapping_traits<Float, Spectrum, M>;

    // Fulfill the Mapping contract
    MappingChain(const Properties &props)
    : Base(props), m_mappings(MapTypes(props)...) {

        (mapchecker::assert_is_mapping<MapTypes>(), ...);

        m_backend = props.get_volume<VolumeBackend>("volume", 1.f);
    }

    Point<Float, 3> apply(const Point<Float, 3>& p) const {
        Point<Float, 3> result = m_to_local * p;
        std::apply([&](const auto &...mappings) {
            ((result = mappings.apply(result)), ...);
        }, m_mappings);
        return result;
    }

    // Volume<Float, Spectrum> implementation
    inline Interaction3f map_interaction(const Interaction3f &it) const {
        Interaction3f it2 = it;
        it2.p = apply(it.p);
        return it2;
    }

    UnpolarizedSpectrum eval(const Interaction3f &it, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::TextureEvaluate, active);
        return BackendDispatch::eval(*m_backend, map_interaction(it), active);
    }

    Float eval_1(const Interaction3f &it, Mask active = true) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::TextureEvaluate, active);
        return BackendDispatch::eval_1(*m_backend, map_interaction(it), active);
    }

    typename Base::PinGuard pin() const override { return m_backend->pin(); }

    void traverse(TraversalCallback *cb) override {
        cb->put("volume", m_backend.get(), ParamFlags::NonDifferentiable);
        Base::traverse(cb);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "MappingChain[" << std::endl
            << "  to_local = " << string::indent(m_to_local, 13) << "," << std::endl
            << "  bbox = " << string::indent(m_bbox) << "," << std::endl
            << "  mappings = [" << std::endl;
        std::apply([&](const auto &...mappings) {
            ((oss << "    " << string::indent(mapchecker::mapping_to_string(mappings), 4)
                  << "," << std::endl), ...);
        }, m_mappings);
        oss << "  ]," << std::endl
            << "  volume = " << string::indent(m_backend) << std::endl
            << "]";
        return oss.str();
    }

    // Extremum structures support: extremum(), resolution(), min() and max()
    std::pair<Float, Float> extremum(BoundingBox3f bbox, Mask local) const override {
        if (dr::any(!local))
            NotImplementedError("MappingChain only supports local bounds");
        BoundingBox3f queried = bbox;
        bbox.clip(BoundingBox3f(Point3f(0.f), Point3f(1.f)));
        auto [min, max] = BackendDispatch::extremum(*m_backend, bbox, local);
        std::apply([&](const auto &...mappings) {
            (([&] {
                using MappingType = std::decay_t<decltype(mappings)>;
                std::tie(min, max) = MappingDispatch<MappingType>::extremum(mappings, queried, min, max);
            }()), ...);
        }, m_mappings);
        return {min, max};
    }

    ScalarVector3i resolution() const override { return BackendDispatch::resolution(*m_backend); }

    ScalarFloat max() const override {
        ScalarFloat result = BackendDispatch::max(*m_backend);
        std::apply([&](const auto &...mappings) {
            (([&] {
                using MappingType = std::decay_t<decltype(mappings)>;
                result = dr::maximum(result, MappingDispatch<MappingType>::max(mappings));
            }()), ...);
        }, m_mappings);
        return result;
    }

    ScalarFloat min() const override {
        ScalarFloat result = BackendDispatch::min(*m_backend);
        std::apply([&](const auto &...mappings) {
            (([&] {
                using MappingType = std::decay_t<decltype(mappings)>;
                result = dr::minimum(result, MappingDispatch<MappingType>::min(mappings));
            }()), ...);
        }, m_mappings);
        return result;
    }

private:
    std::tuple<MapTypes...> m_mappings;
    ref<VolumeBackend> m_backend;
};


// Spherical coords
template <typename Float, typename Spectrum>
struct CubeShellMapping {
public:
    MI_IMPORT_TYPES();

    CubeShellMapping(const Properties &props)
        : m_rmin(props.get<ScalarFloat>("rmin", 0.f)),
          m_rmax(props.get<ScalarFloat>("rmax", 1.f)),
          m_fillmin(props.get<ScalarFloat>("fillmin", 0.f)),
          m_fillmax(props.get<ScalarFloat>("fillmax", 0.f)) {
        if (m_rmin > m_rmax)
            Throw("rmin must be lower than rmax!");
    }

    // Actual core logic of the plugin
    Point3f apply(const Point3f &p) const {
        Float r = dr::norm(p);
        return Point3f(
            (r - m_rmin) / (m_rmax - m_rmin),
            dr::acos(p.z() / r) * dr::InvPi<ScalarFloat>,
            dr::atan2(p.y(), p.x()) * dr::InvTwoPi<ScalarFloat> + 0.5f
        );
    }

    ScalarFloat fillmin() const { return m_fillmin; }
    ScalarFloat fillmax() const { return m_fillmax; }

    std::string to_string() const {
        std::ostringstream oss;
        oss << "CubeShellMapping[rmin=" << m_rmin << ", rmax=" << m_rmax
            << ", fillmin=" << m_fillmin << ", fillmax=" << m_fillmax << "]";
        return oss.str();
    }

private:
    ScalarFloat m_rmin, m_rmax, m_fillmin, m_fillmax;
};

template <typename Float, typename Spectrum>
struct mapping_traits<Float, Spectrum, CubeShellMapping<Float, Spectrum>> {
public:
    MI_IMPORT_TYPES()

    // Require support for fillmin/max
    static ScalarFloat max(const CubeShellMapping<Float, Spectrum> &m) { return m.fillmax(); }
    static ScalarFloat min(const CubeShellMapping<Float, Spectrum> &m) { return m.fillmin(); }

    // Extremum implementation requires mapping specific fillmin/max behavior
    static std::pair<Float, Float>
    extremum(const CubeShellMapping<Float, Spectrum> &m, BoundingBox3f bbox, Float min, Float max) {
        Mask below_rmin = bbox.min.x() < 0.f;
        dr::masked(max, below_rmin) = dr::maximum(max, Float(m.fillmin()));
        dr::masked(min, below_rmin) = dr::minimum(min, Float(m.fillmin()));

        Mask above_rmax = bbox.max.x() > 1.f;
        dr::masked(max, above_rmax) = dr::maximum(max, Float(m.fillmax()));
        dr::masked(min, above_rmax) = dr::minimum(min, Float(m.fillmax()));

        return {min, max};
    }
};


/**
 * Irregular guided Z mapping
 * 
 * Remaps an arbitrary, strictly-increasing sequence of cell boundaries in 
 * [0, 1] read from a nested "boundaries_z" volume, one point per entry.
 *
 * Locating which cell a query falls into is done with a guide table
 * instead of a direct binary search: a little memory is spent at construction
 * on a coarse index that points to buckets of actual boundaries data.
 *
 * Each boundary is assigned to one of a fixed set of buckets by a cheap,
 * order-preserving function of its value. The buckets end up ordered the same
 * way the boundaries themselves are. A query searches only its own bucket,
 * identified using the same cheap function. The vast majority of the data
 * point is thus ruled out.
 *
 * The indexing function used is the integer reinterpretation of a value's
 * IEEE 754 floating-point representation: for positive values, this
 * reinterpretation increases monotonically with the value itself. Discarding
 * its lowest bits then yields a bucket whose width doubles each time the
 * magnitude of the values it covers doubles. This benefits any data whose fine
 * detail sits among its smaller values and coarser detail among its larger
 * ones.
 *
 */

template <typename Float, typename Spectrum>
struct IrregularGuidedZMapping {
public:
    MI_IMPORT_TYPES()
    using FloatStorage  = DynamicBuffer<Float>;
    using UInt32Storage = DynamicBuffer<UInt32>;

    IrregularGuidedZMapping(const Properties &props) {
        const TensorXf &bt = props.get_any<TensorXf>("boundaries_z");
        auto n_points = (uint32_t) bt.size();
        if (n_points < 2)
            Throw("IrregularGuidedZMapping: 'boundaries_z' must have at least 2 "
                  "points (got %d)", (int) n_points);

        std::vector<ScalarFloat> b(n_points);
        for (uint32_t i = 0; i < n_points; ++i)
            b[i] = dr::slice(bt.array(), i);

        for (uint32_t i = 1; i < n_points; ++i)
            if (!(b[i] > b[i - 1]))
                Throw("IrregularGuidedZMapping: 'boundaries_z' must be strictly "
                      "increasing (entry %d: %f <= entry %d: %f)",
                      i, b[i], i - 1, b[i - 1]);
        if (b.front() != 0.f || b.back() != 1.f)
            Throw("IrregularGuidedZMapping: 'boundaries_z' must start at 0 and "
                  "end at 1 (got %f .. %f)", b.front(), b.back());

        m_n_cells = n_points - 1;
        m_boundaries = dr::load<FloatStorage>(b.data(), b.size());

        // Float32 bit pattern of the first entry: deliberately truncated,
        // used only as a bucketing hint.
        m_bits0 = bits_of(b[0]);

        // Define the default shift: the smallest shift that keeps the table
        // size on the same order as the actual number of boundaries.
        ScalarUInt32 bits1 = bits_of((ScalarFloat) 1.f);
        uint32_t default_shift = 0;
        while ((bits1 >> default_shift) > n_points)
            ++default_shift;
        m_shift = props.get<uint32_t>("shift", default_shift);

        std::vector<uint32_t> bin(n_points);
        for (uint32_t j = 0; j < n_points; ++j)
            bin[j] = bucket(b[j], m_bits0, m_shift);

        uint32_t kmax = bin.back();
        std::vector<uint32_t> t(kmax + 2, 0);
        for (uint32_t k = 0, j = 0; k <= kmax + 1; ++k) {
            while (j < n_points && bin[j] < k) ++j;
            t[k] = j;
        }
        m_max_run = 0;
        for (uint32_t k = 0; k <= kmax; ++k)
            m_max_run = dr::maximum(m_max_run, t[k + 1] - t[k]);
        m_table = dr::load<UInt32Storage>(t.data(), t.size());
    }

    /// Float32 bit pattern of a value, reinterpreted as an integer.
    template <typename T>
    static dr::uint32_array_t<T> bits_of(T x) {
        return dr::reinterpret_array<dr::uint32_array_t<T>>(dr::float32_array_t<T>(x));
    }

    /// Bucket index of a value: its bits, relative to bits0, shifted.
    template <typename T>
    static dr::uint32_array_t<T> bucket(T x, ScalarUInt32 bits0, uint32_t shift) {
        return (bits_of(x) - dr::uint32_array_t<T>(bits0)) >> shift;
    }

    Point3f apply(const Point3f &p) const {
        Float x = p.z();
        UInt32 k = bucket(x, m_bits0, m_shift);

        UInt32 t_lo = dr::gather<UInt32>(m_table, k, true);
        UInt32 t_hi = dr::gather<UInt32>(m_table, k + 1u, true);
        UInt32 cell = dr::select(t_lo > 0u, t_lo - 1u, UInt32(0u));
        UInt32 hi   = dr::select(t_hi > 0u, t_hi - 1u, UInt32(0u));

        for (uint32_t step = 0; step < m_max_run; ++step) {
            Mask advance = (cell < hi) &&
                (dr::gather<Float>(m_boundaries, cell + 1u, true) <= x);
            cell = dr::select(advance, cell + 1u, cell);
        }

        Float lo = dr::gather<Float>(m_boundaries, cell, true);
        Float hi_b = dr::gather<Float>(m_boundaries, cell + 1u, true);
        Float frac = dr::clip((x - lo) / dr::maximum(hi_b - lo, 1e-8f), 0.f, 1.f);

        Float inv_n = 1.f / (Float) m_n_cells;
        Float r = (Float(cell) + frac) * inv_n;
        Float cell_hi = (Float(cell) + 1.f) * inv_n;
        Point3f result = p;
        // Avoid rounding issues
        Float margin = dr::minimum(8.f * dr::Epsilon<Float> * cell_hi, 0.5f * inv_n);
        result.z() = dr::minimum(r, cell_hi - margin);
        return result;
    }

    std::string to_string() const {
        std::ostringstream oss;
        oss << "IrregularGuidedZMapping[n_cells=" << m_n_cells
            << ", shift=" << m_shift << ", max_run=" << m_max_run << "]";
        return oss.str();
    }

private:
    FloatStorage  m_boundaries;
    UInt32Storage m_table;
    uint32_t m_n_cells    = 0;
    ScalarUInt32 m_bits0  = 0;
    uint32_t m_shift      = 0;
    uint32_t m_max_run    = 0;
};


/**
 * Irregular Search Mapping
 * 
 * Irregular-grid mapping via a straightforward binary search
 * (dr::binary_search.
 * Remaps an arbitrary, strictly-increasing sequence of cell boundaries in
 * [0, 1] read from a "boundaries_x"/"boundaries_y"/"boundaries_z"
 * volume, chosen by the Axis template parameter.
 * Targets regularly spaced [0, 1] texture coordinate, on that one axis.
 * The other two axes pass through unchanged. Axis is fixed at compile
 * time (0 = x, 1 = y, 2 = z).
 */

template <typename Float, typename Spectrum, size_t Axis>
struct IrregularSearchMapping {
public:
    static_assert(Axis < 3,
        "IrregularSearchMapping: Axis must be 0 (x), 1 (y) or 2 (z)");

    MI_IMPORT_TYPES()
    using FloatStorage = DynamicBuffer<Float>;

    IrregularSearchMapping(const Properties &props) {
        static constexpr const char *key =
            (Axis == 0) ? "boundaries_x" :
            (Axis == 1) ? "boundaries_y" : "boundaries_z";

        const TensorXf &bt = props.get_any<TensorXf>(key);
        auto n_points = (uint32_t) bt.size();
        if (n_points < 2)
            Throw("IrregularSearchMapping: '%s' must have at least 2 "
                  "points (got %d)", key, (int) n_points);

        std::vector<ScalarFloat> b(n_points);
        for (uint32_t i = 0; i < n_points; ++i)
            b[i] = dr::slice(bt.array(), i);

        for (uint32_t i = 1; i < n_points; ++i)
            if (!(b[i] > b[i - 1]))
                Throw("IrregularSearchMapping: '%s' must be strictly "
                      "increasing (entry %d: %f <= entry %d: %f)",
                      key, i, b[i], i - 1, b[i - 1]);
        if (b.front() != 0.f || b.back() != 1.f)
            Throw("IrregularSearchMapping: '%s' must start at 0 and end "
                  "at 1 (got %f .. %f)", key, b.front(), b.back());

        m_n_cells = n_points - 1;
        m_boundaries = dr::load<FloatStorage>(b.data(), b.size());
    }

    Point3f apply(const Point3f &p) const {
        Float x = p[Axis];

        UInt32 cell = searchsorted<UInt32>(m_boundaries, 0u, m_n_cells + 1u, x, true);

        Float lo   = dr::gather<Float>(m_boundaries, cell, true);
        Float hi_b = dr::gather<Float>(m_boundaries, cell + 1u, true);
        Float frac = dr::clip((x - lo) / dr::maximum(hi_b - lo, 1e-8f), 0.f, 1.f);

        Float inv_n    = 1.f / (Float) m_n_cells;
        Float r        = (Float(cell) + frac) * inv_n;
        Float cell_hi  = (Float(cell) + 1.f) * inv_n;
        // Avoid rounding issues
        Float margin = dr::minimum(8.f * dr::Epsilon<Float> * cell_hi, 0.5f * inv_n);

        Point3f result = p;
        result[Axis] = dr::minimum(r, cell_hi - margin);
        return result;
    }

    std::string to_string() const {
        static constexpr const char *axis_name =
            (Axis == 0) ? "x" : (Axis == 1) ? "y" : "z";
        std::ostringstream oss;
        oss << "IrregularSearchMapping[axis=" << axis_name
            << ", n_cells=" << m_n_cells << "]";
        return oss.str();
    }

private:
    /// Specialized call to dr::binary_search for readability
    template <typename Index, typename ArrT, typename BoundT, typename Value, typename MaskT>
    Index searchsorted(const ArrT &arr, BoundT start, BoundT len,
                        Value query, MaskT active) const {
        Index lo = dr::binary_search<Index>(
            start, start + len,
            [&arr, active, query](Index i) DRJIT_INLINE_LAMBDA {
                return active && (dr::gather<Value>(arr, i, active) < query);
            });
        return dr::clip(lo, start + 1u, start + len - 1u) - 1u;
    }

    FloatStorage m_boundaries;
    uint32_t m_n_cells = 0;
};


// A generic irregular Z mapping plugin
// This class can use MI_DECLARE_CLASS and will thus be properly registered for logging functions typically.
// Since MappingChain parent specialization is already a full plugin with a virtual parent interface, this 
// doesn't change the number of vcalls at runtime.
template <typename Float, typename Spectrum>
class IrregularZCoords final : public MappingChain<
        Float, Spectrum, Volume<Float, Spectrum>, IrregularSearchMapping<Float, Spectrum, 2>> {
public:
    using Base = MappingChain<
        Float, Spectrum, Volume<Float, Spectrum>, IrregularSearchMapping<Float, Spectrum, 2>>;
    using Base::Base;  // inherit MappingChain's constructor as-is
    MI_DECLARE_CLASS(IrregularZCoords)
};

// Functionaly equivalent to the sphericalcoords plugin
// This using cannot implement the MI_DECLARE_CLASS but is otherwise sufficient
template <typename Float, typename Spectrum>
using SphericalCoordsChain = MappingChain<
    Float, Spectrum,
    Volume<Float, Spectrum>,
    CubeShellMapping<Float, Spectrum>
>;

// A generic irregular mapping plugin in the 3 dimensions
// This time we declare a mitsuba-like-ish macro to get the best conciseness
// and the MI_DECLARE_CLASS
#define ERDMI_DECLARE_MAPPING_CHAIN(Name, ...)                                \
    template <typename Float, typename Spectrum>                          \
    class Name final : public MappingChain<Float, Spectrum, __VA_ARGS__> { \
    public:                                                                \
        using Base = MappingChain<Float, Spectrum, __VA_ARGS__>;           \
        using Base::Base;                                                  \
        MI_DECLARE_CLASS(Name)                                             \
    };
ERDMI_DECLARE_MAPPING_CHAIN(
    IrregularCoords,
    Volume<Float, Spectrum>, 
    IrregularSearchMapping<Float, Spectrum, 0>,
    IrregularSearchMapping<Float, Spectrum, 1>,
    IrregularSearchMapping<Float, Spectrum, 2>
)

// An irregular Z mapping optimized for an exponential-ish positive profile
// Reusing the macro
ERDMI_DECLARE_MAPPING_CHAIN(
    IrregularGuidedZCoords,
    Volume<Float, Spectrum>, 
    IrregularGuidedZMapping<Float, Spectrum>
);

// A spherical shell and irregular Z profile
template <typename Float, typename Spectrum>
using SphericalIrregularZCoords = MappingChain<
    Float, Spectrum, 
    Volume<Float, Spectrum>,
    CubeShellMapping<Float, Spectrum>,
    IrregularSearchMapping<Float, Spectrum, 2>
>;

// Only one plugin can be exposed per compiled unit, this is a limitation.
// We should thus split this module into a header library imported by several
// thin plugins implementing only the mapping logic and compatible between each
// other by construction
MI_EXPORT_PLUGIN(IrregularGuidedZCoords)

// Otherwise we can write our own machinery for this, if we absolutely want the
// entire library to be held in one file.
#if 0
#define MI_EXPORT_PLUGIN_BEGIN()                                      \
    extern "C" MI_EXPORT void init_plugin(std::string_view name,      \
                                          PluginRegisterFn fn) {

#define MI_EXPORT_PLUGIN_VARIANT(PluginName, Name)                    \
    if (name == PluginName) { MI_REGISTER_PLUGIN(fn, name, Name); return; }

#define MI_EXPORT_PLUGIN_END()                                        \
        Throw("Unknown plugin variant requested: '%s'",               \
              std::string(name).c_str());                             \
    }

MI_EXPORT_PLUGIN_BEGIN()
    MI_EXPORT_PLUGIN_VARIANT("irregularzvolume",       IrregularZCoords)
    MI_EXPORT_PLUGIN_VARIANT("irregularcoordsvolume",  IrregularCoords)
    MI_EXPORT_PLUGIN_VARIANT("irregularguidedzvolume", IrregularGuidedZCoords)
MI_EXPORT_PLUGIN_END()

/* CMake side
function(add_multi_plugin TARGET SOURCE)
  add_plugin(${TARGET} ${SOURCE})
  foreach(NAME ${ARGN})
    add_custom_command(TARGET ${TARGET} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy
              $<TARGET_FILE:${TARGET}>
              ${MI_BINARY_DIR}/plugins/${NAME}$<TARGET_FILE_SUFFIX:${TARGET}>)
    install(FILES $<TARGET_FILE:${TARGET}>
            DESTINATION ${CMAKE_INSTALL_BINDIR}/plugins
            RENAME ${NAME}$<TARGET_FILE_SUFFIX:${TARGET}>)
  endforeach()
endfunction()

add_multi_plugin(irregularzvolume mapcoords.cpp
  irregularcoordsvolume
  irregularguidedzvolume
)
*/
#endif

NAMESPACE_END(mitsuba)
