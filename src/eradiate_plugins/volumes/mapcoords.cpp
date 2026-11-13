#include <mitsuba/core/properties.h>
#include <mitsuba/core/string.h>
#include <mitsuba/core/transform.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/eradiate/extremum.h>

#include <drjit/util.h>

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


template <typename Float, typename Spectrum, typename Backend>
struct backend_traits {

    /*
     * Data storage backend traits
     *
     * Defines the interface to a storage type for the actual volume data. By
     * default the backend type is assumed to follow the
     * mitsuba::Volume<Float, Spectrum> interface.
     *
     */

    MI_IMPORT_TYPES();

    using PinGuard = typename Volume<Float, Spectrum>::PinGuard;
    using Storage  = ref<Backend>;

    static Storage construct(const Properties &props) {
        return props.get_volume<Backend>("volume", 1.f);
    }

    static PinGuard pin(const Storage &storage, const Volume<Float, Spectrum> *) {
        return storage->pin();
    }

    static void traverse(Storage &storage, TraversalCallback *cb) {
        cb->put("volume", storage.get(), ParamFlags::NonDifferentiable);
    }

    static std::string describe(const Storage &storage) {
        return string::indent(storage);
    }

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

    static std::pair<ScalarFloat, ScalarFloat> minmax(const Backend &backend) {
        return { backend.min(), backend.max() };
    }

    static Point3f tex_wrap(const Storage &, const Point3f &p) { return p; }

};


/**
 * Specialization for a Tensor backend
 *
 * This backend uses plain tensors to hold the data and performs no
 * interpolation. It is meant for volumes whose true data is piecewise-constant
 * per cell (one value per layer, no continuous field to filter).
 *
 */

template <typename Float, typename Spectrum>
struct backend_traits<Float, Spectrum, dr::Tensor<DynamicBuffer<Float>>> {
    MI_IMPORT_TYPES();
    using Backend = dr::Tensor<DynamicBuffer<Float>>;
    using PinGuard = typename Volume<Float, Spectrum>::PinGuard;

    struct Storage {
        Backend value;
        dr::WrapMode wrap_mode;
        Storage() : value(), wrap_mode(dr::WrapMode::Clamp) { }
        Storage(Backend &&v, dr::WrapMode wm) : value(std::move(v)), wrap_mode(wm) { }
        Backend &operator*() const { return const_cast<Backend &>(value); }
    };

    static Storage construct(const Properties &props) {
        std::string wrap_mode_str = props.get<std::string>("wrap_mode", "clamp");
        dr::WrapMode wrap_mode = dr::WrapMode::Clamp;
        if (wrap_mode_str == "repeat")
            wrap_mode = dr::WrapMode::Repeat;
        else if (wrap_mode_str == "mirror")
            wrap_mode = dr::WrapMode::Mirror;
        else if (wrap_mode_str != "clamp")
            Throw("mapcoordsvolume: unknown wrap_mode \"%s\"", wrap_mode_str);

        Storage storage(Backend(props.get_any<TensorXf>("data")), wrap_mode);
        if (props.has_property("volume"))
            (void) props.get<ref<Object>>("volume");
        return storage;
    }

    static PinGuard pin(const Storage &, const Volume<Float, Spectrum> *owner) {
        return PinGuard(owner);
    }

    static void traverse(Storage &storage, TraversalCallback *cb) {
        cb->put("data", storage.value, ParamFlags::Differentiable);
    }

    static std::string describe(const Storage &storage) {
        const char *wrap_mode_str =
            storage.wrap_mode == dr::WrapMode::Repeat ? "repeat" :
            storage.wrap_mode == dr::WrapMode::Mirror ? "mirror" : "clamp";
        std::ostringstream oss;
        oss << "FlatTensor[dimensions=" << resolution(*storage)
            << ", wrap_mode=\"" << wrap_mode_str << "\"]";
        return oss.str();
    }

    static inline UnpolarizedSpectrum eval(const Backend &backend, const Interaction3f &it, Mask active) {
        return eval_1(backend, it, active);
    }

    static inline Float eval_1(const Backend &backend, const Interaction3f &it, Mask active) {
        ScalarVector3i res = resolution(backend);
        Vector3i coord = dr::clip(
            dr::floor2int<Vector3i>(it.p * Vector3f(res)), Vector3i(0), Vector3i(res) - 1);
        UInt32 idx = dr::fmadd(dr::fmadd(coord.z(), res.y(), coord.y()), res.x(), coord.x());
        return dr::gather<Float>(backend.array(), idx, active);
    }

    static std::pair<Float, Float> extremum(const Backend &backend, BoundingBox3f bbox, Mask local) {
        if (dr::any(!local))
            NotImplementedError("Flat array backend only supports local bounds");

        bbox.clip(BoundingBox3f(Point3f(0.f), Point3f(1.f)));
        Mask valid = bbox.valid();

        Vector3i res = resolution(backend);
        Vector3i voxel_min = dr::maximum(dr::floor(bbox.min * Vector3f(res)), Vector3i(0));
        Vector3i voxel_max = dr::minimum(dr::floor(bbox.max * Vector3f(res)), res - 1);

        auto z0 = voxel_min.z();
        Float min_val = dr::Infinity<Float>, max_val = -dr::Infinity<Float>;

        std::tie(z0, min_val, max_val) = dr::while_loop(
            std::make_tuple(z0, min_val, max_val),
            [valid, voxel_max](const auto &z, const auto&, const auto&) {
                return valid && (z <= voxel_max.z());
            },
            [&backend, valid, voxel_min, voxel_max, res](auto &z, auto &min_val, auto &max_val) {
                auto y0 = voxel_min.y();
                std::tie(y0, min_val, max_val) = dr::while_loop(
                    std::make_tuple(y0, min_val, max_val),
                    [valid, voxel_max](const auto &y, const auto&, const auto&) {
                        return valid && (y <= voxel_max.y());
                    },
                    [&backend, valid, voxel_min, voxel_max, res, z](auto &y, auto &min_val, auto &max_val) {
                        auto x0 = voxel_min.x();
                        std::tie(x0, min_val, max_val) = dr::while_loop(
                            std::make_tuple(x0, min_val, max_val),
                            [valid, voxel_max](const auto &x, const auto&, const auto&) {
                                return valid && (x <= voxel_max.x());
                            },
                            [&backend, valid, res, z, y](auto &x, auto &min_val, auto &max_val) {
                                UInt32 idx = dr::fmadd(dr::fmadd(UInt32(z), UInt32(res.y()), UInt32(y)),
                                                        UInt32(res.x()), UInt32(x));
                                Float val = dr::gather<Float>(backend.array(), idx, valid);
                                min_val = dr::minimum(min_val, val);
                                max_val = dr::maximum(max_val, val);
                                x += 1;
                            }
                        );
                        y += 1;
                    }
                );
                z += 1;
            }
        );

        min_val = dr::select(valid, min_val, Float(0.f));
        max_val = dr::select(valid, max_val, Float(0.f));
        return { min_val, max_val };
    }

    static ScalarVector3i resolution(const Backend &backend) {
        return { (int) backend.shape(2), (int) backend.shape(1), (int) backend.shape(0) };
    }

    static std::pair<ScalarFloat, ScalarFloat> minmax(const Backend &backend) {
        const auto &arr = dr::detach(backend.array());
        return { (ScalarFloat) dr::min_nested(arr), (ScalarFloat) dr::max_nested(arr) };
    }

    static Float wrap_axis(dr::WrapMode mode, Float x) {
        if (mode == dr::WrapMode::Clamp)
            return dr::clip(x, 0.f, 1.f);
        if (mode == dr::WrapMode::Repeat)
            return x - dr::floor(x);
        Float m = x - 2.f * dr::floor(x * 0.5f);
        return dr::select(m <= 1.f, m, 2.f - m);
    }

    static Point3f tex_wrap(const Storage &storage, const Point3f &p) {
        dr::WrapMode mode = storage.wrap_mode;
        return Point3f(wrap_axis(mode, p.x()), wrap_axis(mode, p.y()), wrap_axis(mode, p.z()));
    }

};

template <typename Float, typename Spectrum, typename VolumeBackend, typename... MapTypes>
class MappingChain : public Volume<Float, Spectrum> {
public:

    /*
     * Mapping chain template class
     *
     * Structure a set of Mappings from world coordinates to the volume backend
     * coordinates. Accepts an arbitrary number of mappings through a variadic
     * MapTypes template parameter.
     *
     * Mappings are applied in chain, in the order of occurence in the MapTypes
     * variadic template parameter.
     *
     * Mappings are constructed one after each other in the order of occurence
     * in the MapTypes variadic template parameter. They are passed the
     * properties of the current constructor to maintain a compact dict
     * interface.
     *
     * The MappingChain performs a transformation from world coordinates
     * applying the Volume<Float, Spectrum> interface's "to_local" attribute
     * before applying the chain of mappings.
     *
     * Extremum properties are defined using the target backend coordinates.
     */

    MI_IMPORT_BASE(Volume, m_to_local, m_bbox)
    MI_IMPORT_TYPES()

    using BackendDispatch = backend_traits<Float, Spectrum, VolumeBackend>;

    using BackendStorage = typename BackendDispatch::Storage;

    MappingChain(const Properties &props)
    : Base(props), m_mappings(MapTypes(props)...), m_backend(BackendDispatch::construct(props)) { }

    typename Base::PinGuard pin() const override {
        return BackendDispatch::pin(m_backend, this);
    }

    VolumeParametrization<ScalarFloat> parametrization() const override {
        const auto &first_mapping = std::get<0>(m_mappings);
        VolumeParametrization<ScalarFloat> param;
        param.to_world = m_to_local.inverse();
        param.flag = first_mapping.coord_flag();
        param.uv_range = first_mapping.coord_range(); // The first mapping decides the chain's
                                                      // external coordinates system.
        return param;
    }

    void traverse(TraversalCallback *cb) override {
        BackendDispatch::traverse(m_backend, cb);
        Base::traverse(cb);
    }

    std::tuple<Point<Float, 3>, Float, Float> apply(const Point<Float, 3>& p, Mask active) const {
        Point<Float, 3> result = m_to_local * p;
        Float extra_min = dr::Infinity<ScalarFloat>, extra_max = -dr::Infinity<ScalarFloat>;
        std::apply([&](const auto &...mappings) {
            ((
                [&] {
                    auto [next, e_min, e_max] = mappings.apply(result, active);
                    result = next;
                    extra_min = dr::minimum(extra_min, e_min);
                    extra_max = dr::maximum(extra_max, e_max);
                }()
            ), ...);
        }, m_mappings);
        result = BackendDispatch::tex_wrap(m_backend, result);
        return { result, extra_min, extra_max };
    }

    inline std::tuple<Interaction3f, Float, Float> map_interaction(const Interaction3f &it, Mask active) const {
        auto [p, e_min, e_max] = apply(it.p, active);
        return { Interaction3f(it.t, it.time, it.wavelengths, p, it.n), e_min, e_max };
    }

    UnpolarizedSpectrum eval(const Interaction3f &it, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::TextureEvaluate, active);
        auto [mapped_it, e_min, e_max] = map_interaction(it, active);
        UnpolarizedSpectrum value = BackendDispatch::eval(*m_backend, mapped_it, active);
        return dr::select(e_min <= e_max, e_min, value);
    }

    Float eval_1(const Interaction3f &it, Mask active = true) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::TextureEvaluate, active);
        auto [mapped_it, e_min, e_max] = map_interaction(it, active);
        Float value = BackendDispatch::eval_1(*m_backend, mapped_it, active);
        return dr::select(e_min <= e_max, e_min, value);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "MappingChain[" << std::endl
            << "  to_local = " << string::indent(m_to_local, 13) << "," << std::endl
            << "  bbox = " << string::indent(m_bbox) << "," << std::endl
            << "  mappings = [" << std::endl;
        std::apply([&](const auto &...mappings) {
            ((oss << "    " << string::indent(mappings.to_string(), 4)
                  << "," << std::endl), ...);
        }, m_mappings);
        oss << "  ]," << std::endl << "  volume = ";
        oss << BackendDispatch::describe(m_backend);
        oss << std::endl << "]";
        return oss.str();
    }

    // Extremum structures support: extremum(), resolution(), min() and max()

    // The query region is narrowed through every mapping in the chain, in
    // the same order apply() chains point evaluation, before the backend is
    // consulted. Each mapping may also contribute an "extra" min/max pair for
    // inputs that bypass the backend entirely (e.g. a constant fill value
    // outside some sub-domain). Those contributions accumulate across the
    // chain and are folded in after the backend query.
    std::pair<Float, Float> extremum(BoundingBox3f bbox, Mask local) const override {
        if (dr::any(!local))
            NotImplementedError("MappingChain only supports local bounds");

        bbox.clip(BoundingBox3f(Point3f(0.f), Point3f(1.f)));

        Float extra_min = dr::Infinity<ScalarFloat>, extra_max = -dr::Infinity<ScalarFloat>;
        std::apply([&](const auto &...mappings) {
            ((
                [&] {
                    auto [region, e_min, e_max] = mappings.extremum(bbox);
                    bbox = region;
                    extra_min = dr::minimum(extra_min, e_min);
                    extra_max = dr::maximum(extra_max, e_max);
                }()
            ), ...);
        }, m_mappings);

        auto [min, max] = BackendDispatch::extremum(*m_backend, bbox, local);
        return { dr::minimum(min, extra_min), dr::maximum(max, extra_max) };
    }

    ScalarVector3i resolution() const override { return BackendDispatch::resolution(*m_backend); }

    ScalarFloat max() const override {
        ScalarFloat result = BackendDispatch::minmax(*m_backend).second;
        std::apply([&](const auto &...mappings) {
            ((result = dr::maximum(result, mappings.max())), ...);
        }, m_mappings);
        return result;
    }

    ScalarFloat min() const override {
        ScalarFloat result = BackendDispatch::minmax(*m_backend).first;
        std::apply([&](const auto &...mappings) {
            ((result = dr::minimum(result, mappings.min())), ...);
        }, m_mappings);
        return result;
    }

private:
    std::tuple<MapTypes...> m_mappings;
    BackendStorage m_backend;
};


/**
 * \brief Given a cell already located along a strictly-increasing boundary
 * array, finish mapping a query x into this axis's output [0,1] fraction.
 *
 */
template <typename Float>
Float irregular_fraction_from_bounds(
    Float inv_n, dr::uint32_array_t<Float> cell, Float x,
    Float lo, Float hi_b
) {
    Float frac = dr::clip((x - lo) / dr::maximum(hi_b - lo, 1e-8f), 0.f, 1.f);

    Float r       = (Float(cell) + frac) * inv_n;
    Float cell_hi = (Float(cell) + 1.f) * inv_n;
    Float margin  = dr::minimum(8.f * dr::Epsilon<Float> * cell_hi, 0.5f * inv_n);
    return dr::minimum(r, cell_hi - margin);
}

template <typename Float>
Float irregular_cell_to_fraction(
    const DynamicBuffer<Float> &boundaries, uint32_t n_cells,
    dr::uint32_array_t<Float> cell, Float x, dr::mask_t<Float> active
) {
    Float lo   = dr::gather<Float>(boundaries, cell, active);
    Float hi_b = dr::gather<Float>(boundaries, cell + 1u, active);
    return irregular_fraction_from_bounds<Float>(1.f / (Float) n_cells, cell, x, lo, hi_b);
}


/**
 * Irregular guided mapping
 *
 * Remaps an arbitrary, strictly-increasing sequence of cell boundaries in
 * [0, 1] read from a "boundaries_x"/"boundaries_y"/"boundaries_z" volume,
 * chosen by the Axis template parameter, one point per entry. The other
 * two axes pass through unchanged. Axis is fixed at compile time (0 = x,
 * 1 = y, 2 = z).
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

template <typename Float, typename Spectrum, size_t Axis>
struct IrregularGuidedMapping {
public:
    static_assert(Axis < 3,
        "IrregularGuidedMapping: Axis must be 0 (x), 1 (y) or 2 (z)");

    MI_IMPORT_TYPES()
    using FloatStorage   = DynamicBuffer<Float>;
    using UInt32Storage  = DynamicBuffer<UInt32>;
    using UInt64Storage  = DynamicBuffer<UInt64>;
    using ScalarUIntBits = dr::uint_array_t<ScalarFloat>;

    IrregularGuidedMapping(const Properties &props) {
        static constexpr const char *key =
            (Axis == 0) ? "boundaries_x" :
            (Axis == 1) ? "boundaries_y" : "boundaries_z";

        const TensorXf &bt = props.get_any<TensorXf>(key);
        auto n_points = (uint32_t) bt.size();
        if (n_points < 2)
            Throw("IrregularGuidedMapping: '%s' must have at least 2 "
                  "points (got %d)", key, (int) n_points);

        std::vector<ScalarFloat> b(n_points);
        for (uint32_t i = 0; i < n_points; ++i)
            b[i] = dr::slice(bt.array(), i);

        for (uint32_t i = 1; i < n_points; ++i)
            if (!(b[i] > b[i - 1]))
                Throw("IrregularGuidedMapping: '%s' must be strictly "
                      "increasing (entry %d: %f <= entry %d: %f)",
                      key, i, b[i], i - 1, b[i - 1]);
        if (b.front() != 0.f || b.back() != 1.f)
            Throw("IrregularGuidedMapping: '%s' must start at 0 and "
                  "end at 1 (got %f .. %f)", key, b.front(), b.back());

        m_n_cells = n_points - 1;
        m_inv_n = 1.f / (Float) m_n_cells;
        m_boundaries = dr::load<FloatStorage>(b.data(), b.size());

        m_bits0 = bits_of(b[0]);

        auto build_table = [&](uint32_t shift, std::vector<uint32_t> &t) -> uint32_t {
            std::vector<uint32_t> bin(n_points);
            for (uint32_t j = 0; j < n_points; ++j)
                bin[j] = (uint32_t) bucket(b[j], m_bits0, shift);

            uint32_t kmax = bin.back();
            t.assign(kmax + 2, 0);
            for (uint32_t k = 0, j = 0; k <= kmax + 1; ++k) {
                while (j < n_points && bin[j] < k) ++j;
                t[k] = j;
            }
            uint32_t max_run = 0;
            for (uint32_t k = 0; k <= kmax; ++k)
                max_run = dr::maximum(max_run, t[k + 1] - t[k]);
            return max_run;
        };

        constexpr uint32_t bit_width = sizeof(ScalarUIntBits) * 8;
        ScalarUIntBits bits_last = bits_of(b.back());
        uint64_t max_table_size = (uint64_t) n_points * 64u;

        uint32_t best_shift = 0, best_max_run = n_points;
        std::vector<uint32_t> scratch;
        for (int s = (int) bit_width - 1; s >= 0; --s) {
            uint32_t shift = (uint32_t) s;
            uint64_t projected_size = (uint64_t) ((bits_last - m_bits0) >> shift) + 2u;
            if (projected_size > max_table_size)
                break;

            uint32_t max_run = build_table(shift, scratch);
            if (max_run < best_max_run) {
                best_max_run = max_run;
                best_shift = shift;
            }
        }

        m_shift = best_shift;

        std::vector<uint32_t> t;
        m_max_run = build_table(m_shift, t);

        std::vector<uint64_t> packed(t.size() - 1);
        for (size_t k = 0; k < packed.size(); ++k)
            packed[k] = (uint64_t) t[k + 1] << 32 | (uint64_t) t[k];
        m_table64 = dr::load<UInt64Storage>(packed.data(), packed.size());

        if constexpr (!dr::is_jit_v<Float>)
            m_boundaries_ptr = m_boundaries.data();

        std::cout << "IrregularGuidedMapping: shift=" << m_shift
                  << ", max boundaries per bucket=" << m_max_run << std::endl;
    }

    template <typename T>
    static dr::uint_array_t<T> bits_of(T x) {
        return dr::reinterpret_array<dr::uint_array_t<T>>(x);
    }

    template <typename T>
    static dr::uint_array_t<T> bucket(T x, ScalarUIntBits bits0, uint32_t shift) {
        return (bits_of(x) - dr::uint_array_t<T>(bits0)) >> shift;
    }

    /// Guide table + binary search lookup, with no coherence assumption.
    UInt32 lookup_cell(Float x, Mask active) const {
        x = dr::clip(x, Float(0.f), dr::OneMinusEpsilon<Float>);

        UInt32 k = UInt32(bucket(x, m_bits0, m_shift));

        UInt64 packed = dr::gather<UInt64>(m_table64, k, active);
        UInt32 t_lo = UInt32(packed);
        UInt32 t_hi = UInt32(packed >> 32);
        UInt32 cell = dr::select(t_lo > 0u, t_lo - 1u, UInt32(0u));
        UInt32 hi = dr::minimum(t_hi, m_n_cells);

        UInt32 found = dr::binary_search<UInt32>(
            cell, hi + 1u,
            [&](UInt32 i) DRJIT_INLINE_LAMBDA {
                if constexpr (!dr::is_jit_v<Float>)
                    return m_boundaries_ptr[i] <= x;
                else
                    return dr::gather<Float>(m_boundaries, i, active) <= x;
            });

        hi = dr::maximum(hi, cell + 1u);
        return dr::clip(found, cell + 1u, hi) - 1u;
    }

    std::tuple<Point3f, Float, Float> apply(const Point3f &p, Mask active) const {
        Float x = p[Axis];
        UInt32 cell = lookup_cell(x, active);
        Point3f result = p;
        if constexpr (!dr::is_jit_v<Float>) {
            Float lo = m_boundaries_ptr[cell];
            Float hi = m_boundaries_ptr[cell + 1u];
            result[Axis] = irregular_fraction_from_bounds<Float>(m_inv_n, cell, x, lo, hi);
        } else {
            result[Axis] = irregular_cell_to_fraction(m_boundaries, m_n_cells, cell, x, active);
        }
        return { result, dr::Infinity<ScalarFloat>, -dr::Infinity<ScalarFloat> };
    }

    Float locate_scalar(Float x) const {
        return irregular_cell_to_fraction(m_boundaries, m_n_cells, lookup_cell(x, true), x, true);
    }

    std::tuple<BoundingBox3f, Float, Float> extremum(BoundingBox3f region) const {
        Float lo = locate_scalar(region.min[Axis]);
        Float hi = locate_scalar(region.max[Axis]);
        region.min[Axis] = dr::minimum(lo, hi);
        region.max[Axis] = dr::maximum(lo, hi);
        return { region, dr::Infinity<ScalarFloat>, -dr::Infinity<ScalarFloat> };
    }

    ScalarFloat max() const { return -dr::Infinity<ScalarFloat>; }
    ScalarFloat min() const { return  dr::Infinity<ScalarFloat>; }

    VolumeCoordFlag coord_flag() const { return VolumeCoordFlag::Grid; }
    ScalarBoundingBox3f coord_range() const {
        return ScalarBoundingBox3f(ScalarPoint3f(0.f), ScalarPoint3f(1.f));
    }

    std::string to_string() const {
        static constexpr const char *axis_name =
            (Axis == 0) ? "x" : (Axis == 1) ? "y" : "z";
        std::ostringstream oss;
        oss << "IrregularGuidedMapping[axis=" << axis_name
            << ", n_cells=" << m_n_cells
            << ", shift=" << m_shift << ", max_run=" << m_max_run << "]";
        return oss.str();
    }

private:
    FloatStorage  m_boundaries;
    UInt64Storage m_table64;
    const ScalarFloat *m_boundaries_ptr = nullptr;
    uint32_t m_n_cells    = 0;
    Float m_inv_n         = 0.f;
    ScalarUIntBits m_bits0 = 0;
    uint32_t m_shift      = 0;
    uint32_t m_max_run    = 0;
};


/**
 * Irregular Search Mapping
 * 
 * Irregular-grid mapping via a straightforward binary search
 * (dr::binary_search).
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

    std::tuple<Point3f, Float, Float> apply(const Point3f &p, Mask active) const {
        Point3f result = p;
        if constexpr (Axis == 2) {
            Float x = p[Axis];
            UInt32 cell = searchsorted<UInt32>(m_boundaries, 0u, m_n_cells + 1u, x, active);
            result[Axis] = irregular_cell_to_fraction(m_boundaries, m_n_cells, cell, x, active);
        } else {
            result[Axis] = locate_scalar(p[Axis], active);
        }
        return { result, dr::Infinity<ScalarFloat>, -dr::Infinity<ScalarFloat> };
    }

    Float locate_scalar(Float x, Mask active = true) const {
        UInt32 cell = searchsorted<UInt32>(m_boundaries, 0u, m_n_cells + 1u, x, active);
        return irregular_cell_to_fraction(m_boundaries, m_n_cells, cell, x, active);
    }

    std::tuple<BoundingBox3f, Float, Float> extremum(BoundingBox3f region) const {
        Float lo = locate_scalar(region.min[Axis]);
        Float hi = locate_scalar(region.max[Axis]);
        region.min[Axis] = dr::minimum(lo, hi);
        region.max[Axis] = dr::maximum(lo, hi);
        return { region, dr::Infinity<ScalarFloat>, -dr::Infinity<ScalarFloat> };
    }

    ScalarFloat max() const { return -dr::Infinity<ScalarFloat>; }
    ScalarFloat min() const { return  dr::Infinity<ScalarFloat>; }

    VolumeCoordFlag coord_flag() const { return VolumeCoordFlag::Grid; }
    ScalarBoundingBox3f coord_range() const {
        return ScalarBoundingBox3f(ScalarPoint3f(0.f), ScalarPoint3f(1.f));
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

    FloatStorage m_boundaries;
    uint32_t m_n_cells = 0;
};


/**
 * Spherical mapping
 *
 * Remaps a Cartesian point directly to normalized spherical
 * coordinates: x = (r - rmin) / (rmax - rmin), y = theta / pi,
 * z = phi / (2 pi) + 1/2, with r the distance from the origin, theta
 * the zenith angle from +Z and phi the azimuth from the (+X, +Z)
 * plane.
 *
 * Points with r outside [rmin, rmax] fall outside this mapping's domain.
 * Rather than a per-sample fallback, their contribution is folded into
 * the fillmin/fillmax "extra" min/max pair extremum() already supports
 * for any mapping that bypasses the backend over part of its domain.
 */

template <typename Float, typename Spectrum>
struct SphericalMapping {
public:
    MI_IMPORT_TYPES()

    SphericalMapping(const Properties &props) {
        m_rmin = props.get<ScalarFloat>("rmin", 0.f);
        m_rmax = props.get<ScalarFloat>("rmax", 1.f);
        if (m_rmin > m_rmax)
            Throw("SphericalMapping: 'rmin' must be lower than 'rmax'");

        m_fillmin = props.get<ScalarFloat>("fillmin", 0.f);
        m_fillmax = props.get<ScalarFloat>("fillmax", 0.f);
    }

    std::tuple<Point3f, Float, Float> apply(const Point3f &p, Mask /*active*/) const {
        Float r = dr::norm(p);
        Point3f result(
            (r - m_rmin) / (m_rmax - m_rmin),
            dr::acos(p.z() / r) * dr::InvPi<ScalarFloat>,
            dr::atan2(p.y(), p.x()) * dr::InvTwoPi<ScalarFloat> + .5f
        );

        Float extra_min = dr::Infinity<ScalarFloat>, extra_max = -dr::Infinity<ScalarFloat>;
        Mask below_rmin = r < m_rmin, above_rmax = r > m_rmax;
        extra_min = dr::select(below_rmin, Float(m_fillmin), extra_min);
        extra_max = dr::select(below_rmin, Float(m_fillmin), extra_max);
        extra_min = dr::select(above_rmax, Float(m_fillmax), extra_min);
        extra_max = dr::select(above_rmax, Float(m_fillmax), extra_max);

        return { result, extra_min, extra_max };
    }

    // Paired with extremum_spherical, the region is already expressed in 
    // this mapping's (r, theta, phi) output space. No forward mapping needed here.
    std::tuple<BoundingBox3f, Float, Float> extremum(BoundingBox3f region) const {
        Float extra_min = dr::Infinity<ScalarFloat>, extra_max = -dr::Infinity<ScalarFloat>;
        Mask below_rmin = region.min.x() < 0.f, above_rmax = region.max.x() > 1.f;
        extra_min = dr::select(below_rmin, Float(m_fillmin), extra_min);
        extra_max = dr::select(below_rmin, Float(m_fillmin), extra_max);
        extra_min = dr::select(above_rmax, Float(m_fillmax), extra_min);
        extra_max = dr::select(above_rmax, Float(m_fillmax), extra_max);

        region.min = dr::clip(region.min, 0.f, 1.f);
        region.max = dr::clip(region.max, 0.f, 1.f);

        return { region, extra_min, extra_max };
    }

    ScalarFloat max() const { return dr::maximum(m_fillmin, m_fillmax); }
    ScalarFloat min() const { return dr::minimum(m_fillmin, m_fillmax); }

    VolumeCoordFlag coord_flag() const { return VolumeCoordFlag::Spherical; }
    ScalarBoundingBox3f coord_range() const {
        return ScalarBoundingBox3f(ScalarPoint3f(m_rmin, 0.f, 0.f), ScalarPoint3f(m_rmax, 1.f, 1.f));
    }

    std::string to_string() const {
        std::ostringstream oss;
        oss << "SphericalMapping[rmin=" << m_rmin << ", rmax=" << m_rmax
            << ", fillmin=" << m_fillmin << ", fillmax=" << m_fillmax << "]";
        return oss.str();
    }

private:
    ScalarFloat m_rmin, m_rmax, m_fillmin, m_fillmax;
};


/**!
.. _volume-mapcoordsvolume:

Irregular coordinate mapping (:monosp:`mapcoordsvolume`)
--------------------------------------------------------

.. pluginparameters:: Common parameters

 * - mapping
   - |string|
   - Selects the remapping strategy. Three values exist: ``guided``,
     ``search`` and ``spherical``. ``guided`` builds a bit-pattern
     guide table at construction, tuned for exponential-like Z
     profiles. ``search`` performs a binary search over the Z
     boundaries instead; it produces identical results to
     ``guided``. ``spherical`` remaps Cartesian points to normalized
     spherical coordinates first, then applies an irregular guided
     remap to the radial axis.
   - —

 * - backend
   - |string|
   - Selects where the remapped density data lives. Two values
     exist: ``volume`` and ``flat``. ``volume`` wraps a nested
     volume plugin and is the default. ``flat`` owns a flat tensor
     directly, with no interpolation. ``mapping = "spherical"``
     only supports ``flat``.
   - —

 * - volume
   - |volume|
   - Nested volume plugin supplying the density data, addressed in
     the remapped coordinates. Required for the ``volume`` backend.
     Unused for the ``flat`` backend.
   - —

 * - data
   - |tensor|
   - Flat 3D tensor of density data, addressed in the remapped
     coordinates with one value per cell and no interpolation.
     Required for the ``flat`` backend. Unused for the ``volume``
     backend.
   - ✓

 * - wrap_mode
   - |string|
   - Boundary handling applied after remapping. Valid values are
     ``clamp``, ``repeat`` and ``mirror``, with ``clamp`` as
     default. Only valid for the ``flat`` backend. Setting it for
     the ``volume`` backend raises an error, since the nested volume
     plugin handles its own wrapping in that case.
   - —

This plugin implements several coordinate mapping strategies from world
frame to local data frame, selected at construction time using the
``mapping`` parameter. The data can be stored as a 3D tensor or a more
generic Mitsuba volume, selected using the ``backend`` parameter.
Each mapping uses a different set of parameters, described below.

.. pluginparameters:: Search/Guided mapping-specific parameters

 * - boundaries_z
   - |tensor|
   - Strictly increasing sequence of Z cell boundaries in [0, 1],
     with at least 2 points. Must start at 0 and end at 1. Defines
     the irregular Z grid. Required for ``guided`` and ``search``
     mappings. Unused for ``spherical``.
   - —

.. pluginparameters:: Spherical mapping-specific parameters

 * - boundaries_x
   - |tensor|
   - Strictly increasing sequence of radial cell boundaries in
     [0, 1], with at least 2 points. Must start at 0 and end at 1.
     Defines the irregular radial grid. Required for ``spherical``,
     unused otherwise.
   - —

 * - fillmin
   - |float|
   - Constant value substituted, at both ``eval()`` and ``extremum()``
     time, for points such that :math:`r < r_\mathrm{min}`. Only
     used by ``spherical``. Default: 0
   - —

 * - fillmax
   - |float|
   - Constant value substituted, same convention as ``fillmin``, for
     points such that :math:`r_\mathrm{max} < r`. Only used by
     ``spherical``. Default: 0
   - —

 * - rmin
   - |float|
   - Radius for the inner limit of the spherical shell, relative to
     the local [0, 1]³ cube's center. The cube's own face sits at
     r = 1. Only used by ``spherical``. Default: 0
   - —

 * - rmax
   - |float|
   - Radius for the outer limit of the spherical shell, same
     convention as ``rmin``. Only used by ``spherical``.
     Default: 1
   - —

*/

template <typename Float, typename Spectrum>
class MapCoordsVolume final : public Volume<Float, Spectrum> {
public:
    MI_IMPORT_BASE(Volume)

    MapCoordsVolume(const Properties &props) : Base(props), m_props(props) {
        m_mapping = props.get<std::string>("mapping");
        m_backend = props.get<std::string>("backend", "volume");
        props.mark_queried("mapping");
        props.mark_queried("backend");

        if (m_mapping == "guided" || m_mapping == "search") {
            if (m_backend == "volume") {
                if (!props.has_property("volume"))
                    Throw("MapCoordsVolume: backend=\"volume\" requires a \"volume\" property");
                props.mark_queried("volume");
            } else if (m_backend == "flat") {
                if (!props.has_property("data"))
                    Throw("MapCoordsVolume: backend=\"flat\" requires a \"data\" property");
                props.mark_queried("data");
            } else {
                Throw("MapCoordsVolume: unknown \"backend\" = \"%s\" (expected "
                      "\"volume\" or \"flat\")", m_backend);
            }

            if (!props.has_property("boundaries_z"))
                Throw("MapCoordsVolume: missing required \"boundaries_z\" property");
            props.mark_queried("boundaries_z");
        } else if (m_mapping == "spherical") {
            if (m_backend != "flat")
                Throw("MapCoordsVolume: mapping=\"spherical\" only supports backend=\"flat\"");

            if (!props.has_property("data"))
                Throw("MapCoordsVolume: backend=\"flat\" requires a \"data\" property");
            props.mark_queried("data");

            if (!props.has_property("boundaries_x"))
                Throw("MapCoordsVolume: mapping=\"spherical\" missing required \"boundaries_x\" property");
            props.mark_queried("boundaries_x");

            for (const char *key : { "rmin", "rmax", "fillmin", "fillmax" })
                if (props.has_property(key))
                    props.mark_queried(key);
        } else {
            Throw("MapCoordsVolume: unknown \"mapping\" = \"%s\" (expected "
                  "\"guided\", \"search\" or \"spherical\")", m_mapping);
        }

        if (props.has_property("wrap_mode")) {
            if (m_backend != "flat")
                Throw("MapCoordsVolume: \"wrap_mode\" only applies to backend=\"flat\"");
            props.mark_queried("wrap_mode");
        }
    }

    std::vector<ref<Object>> expand() const override {
        using VolumeBackend  = Volume<Float, Spectrum>;
        using FlatBackend    = dr::Tensor<DynamicBuffer<Float>>;
        using Guided         = IrregularGuidedMapping<Float, Spectrum, 2>;
        using Search         = IrregularSearchMapping<Float, Spectrum, 2>;
        using Spherical      = SphericalMapping<Float, Spectrum>;
        using GuidedRadius   = IrregularGuidedMapping<Float, Spectrum, 0>;

        ref<Object> result;

        if (m_mapping == "guided") {
            if (m_backend == "volume")
                result = (Object *) new MappingChain<Float, Spectrum, VolumeBackend, Guided>(m_props);
            else
                result = (Object *) new MappingChain<Float, Spectrum, FlatBackend, Guided>(m_props);
        } else if (m_mapping == "search") {
            if (m_backend == "volume")
                result = (Object *) new MappingChain<Float, Spectrum, VolumeBackend, Search>(m_props);
            else
                result = (Object *) new MappingChain<Float, Spectrum, FlatBackend, Search>(m_props);
        } else {
            result = (Object *) new MappingChain<Float, Spectrum, FlatBackend, Spherical, GuidedRadius>(m_props);
        }

        return { result };
    }

    MI_DECLARE_CLASS(MapCoordsVolume)

protected:
    Properties m_props;
    std::string m_mapping, m_backend;
};

MI_EXPORT_PLUGIN(MapCoordsVolume)

NAMESPACE_END(mitsuba)
