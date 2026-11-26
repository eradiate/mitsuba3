#include <mitsuba/core/properties.h>
#include <mitsuba/core/string.h>
#include <mitsuba/core/transform.h>
#include <mitsuba/render/eradiate/extremum.h>
#include <mitsuba/render/volume.h>

#include <drjit/util.h>

#include <vector>

NAMESPACE_BEGIN(mitsuba)

/*
 * MapCoords
 *
 * This file defines a static library to write and chain a set of coordinates
 * mappings. A mapping is defined as an object able to perform a transformation
 * of a point in a 3D coord system to another through an "apply" method. This
 * library provides builtin mappings and tools to combine them or maintain
 * compatitibility with other Mitsuba interfaces.
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

    static PinGuard pin(const Storage &storage,
                        const Volume<Float, Spectrum> *) {
        return storage->pin();
    }

    static void traverse(Storage &storage, TraversalCallback *cb) {
        cb->put("volume", storage.get(), ParamFlags::Differentiable);
    }

    static std::string describe(const Storage &storage) {
        return string::indent(storage);
    }

    static inline UnpolarizedSpectrum
    eval(const Backend &backend, const Interaction3f &it, Mask active) {
        return backend.eval(it, active);
    }

    static inline Float eval_1(const Backend &backend, const Interaction3f &it,
                               Mask active) {
        return backend.eval_1(it, active);
    }

    static std::pair<Float, Float> extremum(const Backend &backend,
                                            BoundingBox3f bbox, Mask local) {
        return backend.extremum(bbox, local);
    }

    static ScalarVector3i dimensions(const Backend &backend) {
        return backend.resolution();
    }

    static std::pair<ScalarFloat, ScalarFloat> minmax(const Backend &backend) {
        return { backend.min(), backend.max() };
    }

    static Point3f tex_wrap(dr::WrapMode, const Point3f &p) { return p; }

    static VolumeParametrization<ScalarFloat>
    parametrization(const Storage &storage, dr::WrapMode) {
        return storage->parametrization();
    }
};

/**
 * Specialization for a Tensor storage backend
 *
 * This backend uses plain tensors to hold the data and performs no
 * interpolation. It is meant for volumes whose true data is piecewise-constant
 * per cell (one value per layer, no continuous field to filter).
 *
 * This backend is single-channel only.
 *
 */

template <typename Float, typename Spectrum>
struct backend_traits<Float, Spectrum, dr::Tensor<DynamicBuffer<Float>>> {
    MI_IMPORT_TYPES();
    using Backend  = dr::Tensor<DynamicBuffer<Float>>;
    using PinGuard = typename Volume<Float, Spectrum>::PinGuard;

    struct Storage {
        Backend value;
        Storage() = default;
        Storage(Backend &&v) : value(std::move(v)) {}
        Backend &operator*() const { return const_cast<Backend &>(value); }
    };

    static Storage construct(const Properties &props) {
        if constexpr (dr::size_v<UnpolarizedSpectrum> != 1)
            Throw("MapCoordsVolume: backend=\"flat\" holds single-channel data "
                  "and requires a variant with a single spectral channel "
                  "(got %zu channels)",
                  (size_t) dr::size_v<UnpolarizedSpectrum>);
        return Storage(Backend(props.get_any<TensorXf>("data")));
    }

    static PinGuard pin(const Storage &, const Volume<Float, Spectrum> *owner) {
        return PinGuard(owner);
    }

    static void traverse(Storage &storage, TraversalCallback *cb) {
        cb->put("data", storage.value, ParamFlags::Differentiable);
    }

    static std::string describe(const Storage &storage) {
        std::ostringstream oss;
        oss << "FlatTensor[dimensions=" << dimensions(*storage) << "]";
        return oss.str();
    }

    static inline UnpolarizedSpectrum
    eval(const Backend &backend, const Interaction3f &it, Mask active) {
        return eval_1(backend, it, active);
    }

    static inline Float eval_1(const Backend &backend, const Interaction3f &it,
                               Mask active) {
        ScalarVector3i res = dimensions(backend);
        Vector3i coord = dr::clip(dr::floor2int<Vector3i>(it.p * Vector3f(res)),
                                  Vector3i(0), Vector3i(res) - 1);
        UInt32 idx     = dr::fmadd(dr::fmadd(coord.z(), res.y(), coord.y()),
                                   res.x(), coord.x());
        return dr::gather<Float>(backend.array(), idx, active);
    }

    static std::pair<Float, Float> extremum(const Backend &backend,
                                            BoundingBox3f bbox, Mask local) {
        if (dr::any(!local))
            NotImplementedError(
                "Flat array backend only supports local bounds");

        bbox.clip(BoundingBox3f(Point3f(0.f), Point3f(1.f)));
        Mask valid = bbox.valid();

        Vector3i res = dimensions(backend);
        Vector3i voxel_min =
            dr::clip(dr::floor2int<Vector3i>(bbox.min * Vector3f(res)),
                     Vector3i(0), res - 1);
        Vector3i voxel_max =
            dr::minimum(dr::floor(bbox.max * Vector3f(res)), res - 1);

        auto z0       = voxel_min.z();
        Float min_val = dr::Infinity<Float>, max_val = -dr::Infinity<Float>;

        std::tie(z0, min_val, max_val) = dr::while_loop(
            std::make_tuple(z0, min_val, max_val),
            [valid, voxel_max](const auto &z, const auto &, const auto &) {
                return valid && (z <= voxel_max.z());
            },
            [&backend, valid, voxel_min, voxel_max, res](auto &z, auto &min_val,
                                                         auto &max_val) {
                auto y0                        = voxel_min.y();
                std::tie(y0, min_val, max_val) = dr::while_loop(
                    std::make_tuple(y0, min_val, max_val),
                    [valid, voxel_max](const auto &y, const auto &,
                                       const auto &) {
                        return valid && (y <= voxel_max.y());
                    },
                    [&backend, valid, voxel_min, voxel_max, res,
                     z](auto &y, auto &min_val, auto &max_val) {
                        auto x0                        = voxel_min.x();
                        std::tie(x0, min_val, max_val) = dr::while_loop(
                            std::make_tuple(x0, min_val, max_val),
                            [valid, voxel_max](const auto &x, const auto &,
                                               const auto &) {
                                return valid && (x <= voxel_max.x());
                            },
                            [&backend, valid, res, z, y](auto &x, auto &min_val,
                                                         auto &max_val) {
                                UInt32 idx = dr::fmadd(
                                    dr::fmadd(UInt32(z), UInt32(res.y()),
                                              UInt32(y)),
                                    UInt32(res.x()), UInt32(x));
                                Float val = dr::gather<Float>(backend.array(),
                                                              idx, valid);
                                min_val   = dr::minimum(min_val, val);
                                max_val   = dr::maximum(max_val, val);
                                x += 1;
                            });
                        y += 1;
                    });
                z += 1;
            });

        min_val = dr::select(valid, min_val, Float(0.f));
        max_val = dr::select(valid, max_val, Float(0.f));
        return { min_val, max_val };
    }

    static ScalarVector3i dimensions(const Backend &backend) {
        return { (int) backend.shape(2), (int) backend.shape(1),
                 (int) backend.shape(0) };
    }

    static std::pair<ScalarFloat, ScalarFloat> minmax(const Backend &backend) {
        const auto &arr = dr::detach(backend.array());
        return { (ScalarFloat) dr::min_nested(arr),
                 (ScalarFloat) dr::max_nested(arr) };
    }

    static Float wrap_axis(dr::WrapMode mode, Float x) {
        if (mode == dr::WrapMode::Clamp)
            return dr::clip(x, 0.f, 1.f);
        if (mode == dr::WrapMode::Repeat)
            return x - dr::floor(x);
        Float m = x - 2.f * dr::floor(x * 0.5f);
        return dr::select(m <= 1.f, m, 2.f - m);
    }

    static Point3f tex_wrap(dr::WrapMode mode, const Point3f &p) {
        return Point3f(wrap_axis(mode, p.x()), wrap_axis(mode, p.y()),
                       wrap_axis(mode, p.z()));
    }

    static VolumeParametrization<ScalarFloat>
    parametrization(const Storage &, dr::WrapMode mode) {
        VolumeParametrization<ScalarFloat> param;
        param.wrap      = true;
        param.wrap_mode = mode;
        return param;
    }
};

template <typename Float, typename Spectrum, typename VolumeBackend,
          typename... MapTypes>
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
     * Mappings are constructed one after the other in the order of occurence
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
        : Base(props), m_mappings(MapTypes(props)...),
          m_backend(BackendDispatch::construct(props)),
          m_wrap_mode(parse_wrap_mode(
              props.get<std::string_view>("wrap_mode", "clamp"))) {
        m_resolution = props.get<ScalarVector3i>(
            "resolution", BackendDispatch::dimensions(*m_backend));
        if (dr::any(m_resolution < 1))
            Throw("mapcoordsvolume: \"resolution\" must be at least 1 along "
                  "each axis (got %s)",
                  m_resolution);
    }

    typename Base::PinGuard pin() const override {
        return BackendDispatch::pin(m_backend, this);
    }

    VolumeParametrization<ScalarFloat> parametrization() const override {
        const auto &first_mapping = std::get<0>(m_mappings);
        VolumeParametrization<ScalarFloat> param =
            BackendDispatch::parametrization(m_backend, m_wrap_mode);
        param.to_world = m_to_local.inverse();
        param.flag     = first_mapping.coord_flag();
        param.uv_range =
            first_mapping.coord_range(); // The first mapping decides the
                                         // chain's external coordinates system.
        return param;
    }

    void traverse(TraversalCallback *cb) override {
        BackendDispatch::traverse(m_backend, cb);
        Base::traverse(cb);
    }

    // Returns the point in backend coordinates, and a value with a mask for
    // points that have no data in the backend (e.g. outside [rmin, rmax] in
    // SphericalMapping). eval() returns that value where the mask is set.
    // The first mapping that sets the mask decides the value.
    std::tuple<Point<Float, 3>, Float, Mask> apply(const Point<Float, 3> &p,
                                                   Mask active) const {
        Point<Float, 3> result = m_to_local * p;
        Float fill_value       = 0.f;
        Mask filled            = false;
        std::apply(
            [&](const auto &...mappings) {
                (([&] {
                     auto [next, v, f] = mappings.apply(result, active);
                     result            = next;
                     fill_value        = dr::select(filled, fill_value, v);
                     filled |= f;
                 }()),
                 ...);
            },
            m_mappings);
        result = BackendDispatch::tex_wrap(m_wrap_mode, result);
        return { result, fill_value, filled };
    }

    inline std::tuple<Interaction3f, Float, Mask>
    map_interaction(const Interaction3f &it, Mask active) const {
        auto [p, fill_value, filled] = apply(it.p, active);
        return { Interaction3f(it.t, it.time, it.wavelengths, p, it.n),
                 fill_value, filled };
    }

    UnpolarizedSpectrum eval(const Interaction3f &it,
                             Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::TextureEvaluate, active);
        auto [mapped_it, fill_value, filled] = map_interaction(it, active);
        UnpolarizedSpectrum value =
            BackendDispatch::eval(*m_backend, mapped_it, active && !filled);
        return dr::select(filled, fill_value, value);
    }

    Float eval_1(const Interaction3f &it, Mask active = true) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::TextureEvaluate, active);
        auto [mapped_it, fill_value, filled] = map_interaction(it, active);
        Float value =
            BackendDispatch::eval_1(*m_backend, mapped_it, active && !filled);
        return dr::select(filled, fill_value, value);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "MappingChain[" << std::endl
            << "  to_local = " << string::indent(m_to_local, 13) << ","
            << std::endl
            << "  bbox = " << string::indent(m_bbox) << "," << std::endl
            << "  mappings = [" << std::endl;
        std::apply(
            [&](const auto &...mappings) {
                ((oss << "    " << string::indent(mappings.to_string(), 4)
                      << "," << std::endl),
                 ...);
            },
            m_mappings);
        oss << "  ]," << std::endl
            << "  wrap_mode = " << wrap_mode_name(m_wrap_mode) << ","
            << std::endl
            << "  volume = ";
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
    std::pair<Float, Float> extremum(BoundingBox3f bbox,
                                     Mask local) const override {
        if (dr::any(!local))
            NotImplementedError("MappingChain only supports local bounds");

        bbox.clip(BoundingBox3f(Point3f(0.f), Point3f(1.f)));

        Float extra_min = dr::Infinity<ScalarFloat>,
              extra_max = -dr::Infinity<ScalarFloat>;
        std::apply(
            [&](const auto &...mappings) {
                (([&] {
                     auto [region, e_min, e_max] = mappings.extremum(bbox);
                     bbox                        = region;
                     extra_min = dr::minimum(extra_min, e_min);
                     extra_max = dr::maximum(extra_max, e_max);
                 }()),
                 ...);
            },
            m_mappings);

        auto [min, max] = BackendDispatch::extremum(*m_backend, bbox, local);
        return { dr::minimum(min, extra_min), dr::maximum(max, extra_max) };
    }

    ScalarVector3i resolution() const override { return m_resolution; }

    ScalarFloat max() const override {
        ScalarFloat result = BackendDispatch::minmax(*m_backend).second;
        std::apply(
            [&](const auto &...mappings) {
                ((result = dr::maximum(result, mappings.max())), ...);
            },
            m_mappings);
        return result;
    }

    ScalarFloat min() const override {
        ScalarFloat result = BackendDispatch::minmax(*m_backend).first;
        std::apply(
            [&](const auto &...mappings) {
                ((result = dr::minimum(result, mappings.min())), ...);
            },
            m_mappings);
        return result;
    }

private:
    static dr::WrapMode parse_wrap_mode(std::string_view name) {
        if (name == "clamp")
            return dr::WrapMode::Clamp;
        if (name == "repeat")
            return dr::WrapMode::Repeat;
        if (name == "mirror")
            return dr::WrapMode::Mirror;
        Throw("mapcoordsvolume: unknown wrap_mode \"%s\"", name);
    }

    static const char *wrap_mode_name(dr::WrapMode mode) {
        switch (mode) {
            case dr::WrapMode::Repeat:
                return "repeat";
            case dr::WrapMode::Mirror:
                return "mirror";
            default:
                return "clamp";
        }
    }

    std::tuple<MapTypes...> m_mappings;
    BackendStorage m_backend;
    dr::WrapMode m_wrap_mode;
    ScalarVector3i m_resolution;
};

/**
 * \brief Given a cell already located along a strictly-increasing boundary
 * array, finish mapping a query x into this axis's output [0,1] fraction.
 *
 */
template <typename Float>
Float irregular_fraction_from_bounds(Float inv_n,
                                     dr::uint32_array_t<Float> cell, Float x,
                                     Float lo, Float hi_b) {
    Float frac = dr::clip((x - lo) / dr::maximum(hi_b - lo, 1e-8f), 0.f, 1.f);

    Float r       = (Float(cell) + frac) * inv_n;
    Float cell_hi = (Float(cell) + 1.f) * inv_n;
    Float margin =
        dr::minimum(8.f * dr::Epsilon<Float> * cell_hi, 0.5f * inv_n);
    return dr::minimum(r, cell_hi - margin);
}

template <typename Float>
Float irregular_cell_to_fraction(const DynamicBuffer<Float> &boundaries,
                                 uint32_t n_cells,
                                 dr::uint32_array_t<Float> cell, Float x,
                                 dr::mask_t<Float> active) {
    Float lo   = dr::gather<Float>(boundaries, cell, active);
    Float hi_b = dr::gather<Float>(boundaries, cell + 1u, active);
    return irregular_fraction_from_bounds<Float>(1.f / (Float) n_cells, cell, x,
                                                 lo, hi_b);
}

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
        static constexpr const char *key = (Axis == 0)   ? "boundaries_x"
                                           : (Axis == 1) ? "boundaries_y"
                                                         : "boundaries_z";

        const TensorXf &bt = props.get_any<TensorXf>(key);
        auto n_points      = (uint32_t) bt.size();
        if (n_points < 2)
            Throw("IrregularSearchMapping: '%s' must have at least 2 "
                  "points (got %d)",
                  key, (int) n_points);

        std::vector<ScalarFloat> b(n_points);
        dr::store(b.data(), bt.array());

        for (uint32_t i = 1; i < n_points; ++i)
            if (!(b[i] > b[i - 1]))
                Throw("IrregularSearchMapping: '%s' must be strictly "
                      "increasing (entry %d: %f <= entry %d: %f)",
                      key, i, b[i], i - 1, b[i - 1]);
        if (b.front() != 0.f || b.back() != 1.f)
            Throw("IrregularSearchMapping: '%s' must start at 0 and end "
                  "at 1 (got %f .. %f)",
                  key, b.front(), b.back());

        m_n_cells    = n_points - 1;
        m_boundaries = dr::load<FloatStorage>(b.data(), b.size());
    }

    template <typename Index, typename ArrT, typename BoundT, typename Value,
              typename MaskT>
    Index searchsorted(const ArrT &arr, BoundT start, BoundT len, Value query,
                       MaskT active) const {
        Index lo = dr::binary_search<Index>(
            start, start + len,
            [&arr, active, query](Index i) DRJIT_INLINE_LAMBDA {
                return active && (dr::gather<Value>(arr, i, active) < query);
            });
        return dr::clip(lo, start + 1u, start + len - 1u) - 1u;
    }

    std::tuple<Point3f, Float, Mask> apply(const Point3f &p,
                                           Mask active) const {
        Point3f result = p;
        result[Axis]   = locate_scalar(p[Axis], active);
        return { result, 0.f, false };
    }

    Float locate_scalar(Float x, Mask active = true) const {
        UInt32 cell =
            searchsorted<UInt32>(m_boundaries, 0u, m_n_cells + 1u, x, active);
        return irregular_cell_to_fraction(m_boundaries, m_n_cells, cell, x,
                                          active);
    }

    std::tuple<BoundingBox3f, Float, Float>
    extremum(BoundingBox3f region) const {
        Float lo         = locate_scalar(region.min[Axis]);
        Float hi         = locate_scalar(region.max[Axis]);
        region.min[Axis] = dr::minimum(lo, hi);
        region.max[Axis] = dr::maximum(lo, hi);
        return { region, dr::Infinity<ScalarFloat>,
                 -dr::Infinity<ScalarFloat> };
    }

    ScalarFloat max() const { return -dr::Infinity<ScalarFloat>; }
    ScalarFloat min() const { return dr::Infinity<ScalarFloat>; }

    VolumeCoordFlag coord_flag() const { return VolumeCoordFlag::Grid; }
    ScalarBoundingBox3f coord_range() const {
        return ScalarBoundingBox3f(ScalarPoint3f(0.f), ScalarPoint3f(1.f));
    }

    std::string to_string() const {
        static constexpr const char *axis_name = (Axis == 0)   ? "x"
                                                 : (Axis == 1) ? "y"
                                                               : "z";
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
 * The ``uv_to_spherical`` transform maps the data's [0, 1]^3 cube to these
 * spherical coordinates. Its inverse is applied to the remapped point. It
 * places data that covers part of the sphere.
 *
 * Points with r outside [rmin, rmax] fall outside this mapping's domain.
 * Rather than a per-sample fallback, their contribution is folded into
 * the fillmin/fillmax "extra" min/max pair extremum() already supports
 * for any mapping that bypasses the backend over part of its domain.
 */

template <typename Float, typename Spectrum> struct SphericalMapping {
public:
    MI_IMPORT_TYPES()

    SphericalMapping(const Properties &props) {
        m_rmin = props.get<ScalarFloat>("rmin", 0.f);
        m_rmax = props.get<ScalarFloat>("rmax", 1.f);
        if (m_rmin >= m_rmax)
            Throw("SphericalMapping: 'rmin' must be lower than 'rmax'");
        m_inv_dr = 1.f / (m_rmax - m_rmin);

        m_fillmin = props.get<ScalarFloat>("fillmin", 0.f);
        m_fillmax = props.get<ScalarFloat>("fillmax", 0.f);

        m_spherical_to_uv =
            props
                .get<ScalarAffineTransform4f>("uv_to_spherical",
                                              ScalarAffineTransform4f())
                .inverse();
    }

    std::tuple<Point3f, Float, Mask> apply(const Point3f &p,
                                           Mask /*active*/) const {
        Float r = dr::norm(p);
        Point3f result =
            m_spherical_to_uv *
            Point3f((r - m_rmin) * m_inv_dr,
                    dr::atan2(dr::sqrt(dr::square(p.x()) + dr::square(p.y())),
                              p.z()) *
                        dr::InvPi<ScalarFloat>,
                    dr::atan2(p.y(), p.x()) * dr::InvTwoPi<ScalarFloat> + .5f);

        Mask below_rmin = r<m_rmin, above_rmax = r> m_rmax;
        Float fill_value =
            dr::select(below_rmin, Float(m_fillmin), Float(m_fillmax));

        return { result, fill_value, below_rmin || above_rmax };
    }

    // Paired with extremum_spherical, the region is already expressed in
    // this mapping's (r, theta, phi) output space. No forward mapping needed
    // here.
    std::tuple<BoundingBox3f, Float, Float>
    extremum(BoundingBox3f region) const {
        Float extra_min = dr::Infinity<ScalarFloat>,
              extra_max = -dr::Infinity<ScalarFloat>;
        Mask below_rmin = region.min.x() < 0.f,
             above_rmax = region.max.x() > 1.f;
        extra_min       = dr::select(below_rmin, Float(m_fillmin), extra_min);
        extra_max       = dr::select(below_rmin, Float(m_fillmin), extra_max);
        extra_min       = dr::select(above_rmax, Float(m_fillmax), extra_min);
        extra_max       = dr::select(above_rmax, Float(m_fillmax), extra_max);

        BoundingBox3f uv;
        for (size_t i = 0; i < 8; ++i)
            uv.expand(m_spherical_to_uv * region.corner(i));
        uv.min = dr::clip(uv.min, 0.f, 1.f);
        uv.max = dr::clip(uv.max, 0.f, 1.f);

        return { uv, extra_min, extra_max };
    }

    ScalarFloat max() const { return dr::maximum(m_fillmin, m_fillmax); }
    ScalarFloat min() const { return dr::minimum(m_fillmin, m_fillmax); }

    VolumeCoordFlag coord_flag() const { return VolumeCoordFlag::Spherical; }
    ScalarBoundingBox3f coord_range() const {
        return ScalarBoundingBox3f(ScalarPoint3f(m_rmin, 0.f, 0.f),
                                   ScalarPoint3f(m_rmax, 1.f, 1.f));
    }

    std::string to_string() const {
        std::ostringstream oss;
        oss << "SphericalMapping[rmin=" << m_rmin << ", rmax=" << m_rmax
            << ", fillmin=" << m_fillmin << ", fillmax=" << m_fillmax << "]";
        return oss.str();
    }

private:
    ScalarFloat m_rmin, m_rmax, m_inv_dr, m_fillmin, m_fillmax;
    ScalarAffineTransform4f m_spherical_to_uv;
};

/**!
.. _volume-mapcoordsvolume:

Irregular coordinate mapping (:monosp:`mapcoordsvolume`)
--------------------------------------------------------

.. pluginparameters:: Common parameters

 * - mapping
   - |string|
   - Selects the remapping strategy. Two values exist: ``search`` and
     ``spherical``. ``search`` locates the Z cell with a binary search
     over the Z boundaries. ``spherical`` remaps Cartesian points to
     normalized spherical coordinates first, then locates the radial
     cell with a binary search over the radial boundaries.
   - —

 * - backend
   - |string|
   - Selects where the remapped density data lives. Two values
     exist: ``volume`` and ``flat``. ``volume`` wraps a nested
     volume plugin and is the default. ``flat`` owns a flat tensor
     directly, with no interpolation. ``mapping = "spherical"``
     only supports ``flat``. The flat backend is single-channel only.
   - —

 * - volume
   - |volume|
   - Nested volume plugin supplying the density data, addressed in
     the remapped coordinates. Required for the ``volume`` backend.
     Unused for the ``flat`` backend.
   - |exposed|, |differentiable|

 * - data
   - |tensor|
   - Flat 3D tensor of density data, addressed in the remapped
     coordinates with one value per cell and no interpolation.
     Required for the ``flat`` backend. Unused for the ``volume``
     backend.
   - |exposed|, |differentiable|

 * - wrap_mode
   - |string|
   - Boundary handling applied after remapping. Valid values are
     ``clamp``, ``repeat`` and ``mirror``, with ``clamp`` as
     default. Only valid for the ``flat`` backend. Setting it for
     the ``volume`` backend raises an error, since the nested volume
     plugin handles its own wrapping in that case.
   - —

 * - resolution
   - |vector|
   - Resolution reported to extremum structures, in the coordinates
     before remapping. ``gridextremum`` caps its requested resolution
     at this value. Defaults to the dimensions of the density data.
   - —

This plugin implements several coordinate mapping strategies from world
frame to local data frame, selected at construction time using the
``mapping`` parameter. The data can be stored as a 3D tensor or a more
generic Mitsuba volume, selected using the ``backend`` parameter.
Each mapping uses a different set of parameters, described below.

.. pluginparameters:: Search mapping-specific parameters

 * - boundaries_z
   - |tensor|
   - Strictly increasing sequence of Z cell boundaries in [0, 1],
     with at least 2 points. Must start at 0 and end at 1. Defines
     the irregular Z grid. Required for ``search``. Unused for
     ``spherical``.
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

 * - uv_to_spherical
   - |transform|
   - Transform from the data's [0, 1]³ cube to the normalized spherical
     coordinates (r, θ / π, φ / 2π + 1/2). Places data that covers part
     of the sphere. Only used by ``spherical``. Default: identity
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

        if (m_mapping == "search") {
            if (m_backend == "volume") {
                if (!props.has_property("volume"))
                    Throw("MapCoordsVolume: backend=\"volume\" requires a "
                          "\"volume\" property");
                props.mark_queried("volume");
            } else if (m_backend == "flat") {
                if (!props.has_property("data"))
                    Throw("MapCoordsVolume: backend=\"flat\" requires a "
                          "\"data\" property");
                props.mark_queried("data");
            } else {
                Throw("MapCoordsVolume: unknown \"backend\" = \"%s\" (expected "
                      "\"volume\" or \"flat\")",
                      m_backend);
            }

            if (!props.has_property("boundaries_z"))
                Throw("MapCoordsVolume: missing required \"boundaries_z\" "
                      "property");
            props.mark_queried("boundaries_z");
        } else if (m_mapping == "spherical") {
            if (m_backend != "flat")
                Throw("MapCoordsVolume: mapping=\"spherical\" only supports "
                      "backend=\"flat\"");

            if (!props.has_property("data"))
                Throw("MapCoordsVolume: backend=\"flat\" requires a \"data\" "
                      "property");
            props.mark_queried("data");

            if (!props.has_property("boundaries_x"))
                Throw("MapCoordsVolume: mapping=\"spherical\" missing required "
                      "\"boundaries_x\" property");
            props.mark_queried("boundaries_x");

            for (const char *key :
                 { "rmin", "rmax", "fillmin", "fillmax", "uv_to_spherical" })
                if (props.has_property(key))
                    props.mark_queried(key);
        } else {
            Throw("MapCoordsVolume: unknown \"mapping\" = \"%s\" (expected "
                  "\"search\" or \"spherical\")",
                  m_mapping);
        }

        if (props.has_property("wrap_mode")) {
            if (m_backend != "flat")
                Throw("MapCoordsVolume: \"wrap_mode\" only applies to "
                      "backend=\"flat\"");
            props.mark_queried("wrap_mode");
        }

        if (props.has_property("resolution"))
            props.mark_queried("resolution");
    }

    std::vector<ref<Object>> expand() const override {
        using VolumeBackend = Volume<Float, Spectrum>;
        using FlatBackend   = dr::Tensor<DynamicBuffer<Float>>;
        using Search        = IrregularSearchMapping<Float, Spectrum, 2>;
        using Spherical     = SphericalMapping<Float, Spectrum>;
        using SearchRadius  = IrregularSearchMapping<Float, Spectrum, 0>;

        ref<Object> result;

        if (m_mapping == "search") {
            if (m_backend == "volume")
                result =
                    (Object *) new MappingChain<Float, Spectrum, VolumeBackend,
                                                Search>(m_props);
            else
                result =
                    (Object *) new MappingChain<Float, Spectrum, FlatBackend,
                                                Search>(m_props);
        } else {
            result =
                (Object *) new MappingChain<Float, Spectrum, FlatBackend,
                                            Spherical, SearchRadius>(m_props);
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
