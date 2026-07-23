#include <mitsuba/core/properties.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/eradiate/extremum.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/volumegrid.h>
#include <nanothread/nanothread.h>

NAMESPACE_BEGIN(mitsuba)

/**!
.. _extremum-extremum_grid:

Extremum grid structure (:monosp:`extremum_grid`)
-------------------------------------------------

.. pluginparameters::

 * - resolution
   - |vector|
   - Grid resolution along the XYZ axis. Does not have to be a multiple of
     the underlying volume. Set to [0,0,0] to trigger an adaptive resolution
     routine. Default: [1,1,1]

This plugin creates a regular grid structure storing local extremum values for
efficient delta tracking in heterogeneous media. All geometric information
(domain, extinction volume, scale) is provided by the owning medium through
\ref ExtremumStructure::build(). The grid is constructed by querying the
volume's extrema over each grid cell, in world space.

The grid covers the intersection of the medium domain and the volume's
bounding box. Queries for border cells are extended outward to the domain
edge, so that the stored bounds remain valid over the whole domain for any
volume edge behavior (clamp, repeat, mirror, none).

At runtime, DDA (Digital Differential Analyzer) traversal through the grid provides
tight-fitting local extrema, dramatically reducing null collisions in media with
high spatial variance (*e.g.* clouds, fog).
*/

template <typename Float, typename Spectrum>
class ExtremumGrid final : public ExtremumStructure<Float, Spectrum> {
public:
    MI_IMPORT_BASE(ExtremumStructure, m_bbox, set_domain)
    MI_IMPORT_TYPES(Volume)

    using TrackingStateType    = TrackingState<Float, Spectrum>;
    using TrackingFunctionType = TrackingFunction<Float, Spectrum>;
    using FloatStorage         = DynamicBuffer<Float>;

    ExtremumGrid(const Properties &props) : Base(props) {
        m_resolution = props.get<ScalarVector3i>("resolution", ScalarVector3i(1,1,1));
        m_adaptive   = dr::all(m_resolution <= ScalarVector3i(0));
    }

    void build(const ScalarBoundingBox3f &domain, const Volume *volume,
               ScalarFloat scale) override {
        Log(Warn, "Extremum grid build, scale: %f", scale);

        set_domain(domain, volume);
        m_scale = scale;

        // The grid only covers the region where the volume defines data;
        // border-cell extension makes its bounds valid over the rest of the
        // domain.
        m_extent = domain;
        m_extent.clip(volume->bbox());
        if (!m_extent.valid())
            Throw("extremum_grid: the medium domain %s does not intersect "
                  "the volume's bounding box %s!", domain,
                  volume->bbox());

        m_to_local = ScalarAffineTransform4f::scale(dr::rcp(m_extent.extents())) *
                     ScalarAffineTransform4f::translate(-ScalarVector3f(m_extent.min));

        if (m_adaptive)
            m_resolution = find_resolution(volume);

        build_grid(volume, m_resolution);
    }

    TrackingStateType traverse_extremum(
        const Ray3f &ray,
        Float mint,
        Float maxt,
        UInt32 channel,
        TrackingStateType state,
        TrackingFunctionType* func,
        Mask active
    ) const override {
        return traverse_dda(
            func,
            state,
            ray,
            mint,
            maxt,
            channel,
            active
        );
    }

    ExtremumSegment next_segment(const Ray3f &ray, Float t,
                                 Mask active) const override {
        auto [hit, d0, d1] = m_bbox.ray_intersect(ray);

        // Forward nudge: select the cell at p(t + eps), but measure exit
        // distances from exact p(t)
        Float eps = dr::maximum(dr::abs(t), 1.f) * math::RayEpsilon<Float>;
        Float tq  = t + eps;

        Mask before = hit && (tq < d0);
        Mask inside = hit && !before && (tq < d1);

        // Grid-coordinate ray in [0, res]^3 over the grid extent
        Vector3f res(m_resolution);
        Point3f  o = (m_to_local * ray.o) * res;
        Vector3f d = (m_to_local * ray.d) * res;
        Vector3i pi = dr::floor2int<Vector3i>(dr::fmadd(d, tq, o));

        // Per-axis next cell boundary. Beyond the grid extent, cells are
        // semi-infinite: the only boundary along such an axis is the grid
        // face (when approaching). Stale boundaries (at or behind tq, e.g.
        // when leaving the grid) resolve to +inf.
        Float t_exit = d1;
        for (size_t i = 0; i < 3; ++i) {
            Int32 b_idx = dr::select(d[i] > 0.f, pi[i] + 1, pi[i]);
            Float b     = Float(dr::clip(b_idx, 0, m_resolution[i]));
            Float ti    = dr::select(d[i] != 0.f, (b - o[i]) / d[i],
                                     dr::Infinity<Float>);
            dr::masked(ti, ti <= tq) = dr::Infinity<Float>;
            t_exit = dr::minimum(t_exit, ti);
        }

        Vector3i piw = dr::clip(pi, 0, m_resolution - 1);
        UInt32 idx = dr::fmadd(dr::fmadd(piw.z(), m_resolution.y(), piw.y()),
                               m_resolution.x(), piw.x());
        Vector2f value =
            dr::gather<Vector2f>(m_extremum_grid, idx, active && inside);

        Float maxt = dr::select(
            inside, dr::maximum(t_exit, tq),
            dr::select(before, d0, dr::Infinity<Float>));

        return ExtremumSegment(t, maxt,
                               dr::select(inside, value, Vector2f(0.f)));
    }

    std::tuple<Float, Float> eval_1(const Interaction3f &it,
                                    Mask active) const override {
        Point3f po = m_to_local * it.p;
        Vector3i piw = dr::clip(dr::floor2int<Vector3i>(po * m_resolution),
                                0, m_resolution - 1);

        UInt32 idx = dr::fmadd(dr::fmadd(piw.z(), m_resolution.y(), piw.y()),
                               m_resolution.x(), piw.x());
        Vector2f extremum = dr::gather<Vector2f>(m_extremum_grid, idx, active);
        return { extremum.x(), extremum.y() };
    }

    void traverse(TraversalCallback *cb) override {
        cb->put("data", m_extremum_grid, ParamFlags::NonDifferentiable);
        cb->put("resolution", m_resolution, ParamFlags::NonDifferentiable);
        Base::traverse(cb);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "ExtremumGrid[" << std::endl
            << "  resolution = " << m_resolution << "," << std::endl
            << "  bbox = " << m_bbox << "," << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS(ExtremumGrid)

private:

    /**
     * \brief Find the best resolution for a given volume.
     *
     * This method performs a ternary search to find the best extremum grid
     * resoluton for a given volume. Note that this is currently a costly
     * operation.
     */
    ScalarVector3i find_resolution(const Volume *volume) {
        using ScalarIndex        = DynamicBuffer<ScalarUInt32>;
        using ScalarFloatStorage = DynamicBuffer<ScalarFloat>;

        ScalarVector3i result   = dr::full<ScalarVector3i>(1);
        ScalarVector3i fine_res = volume->resolution();
        ScalarVector3f dim      = m_extent.extents();

        ScalarFloat denom =
            dim.x() * dim.y() + dim.x() * dim.z() + dim.y() * dim.z();

        // Cost function as defined in Yue et al. 2011
        auto cost_fn = [denom, dim](FloatStorage &grid,
                                    ScalarVector3i resolution) {
            // Select only the majorant from the extremum grid
            ScalarIndex idx =
                dr::arange<ScalarIndex>(dr::prod(resolution)) * 2 + 1;
            ScalarFloatStorage majorant =
                dr::gather<ScalarFloatStorage>(grid, idx);

            ScalarFloat res_prod = dr::prod(resolution);
            ScalarFloat dim_prod = dr::prod(dim);
            ScalarFloat grid_sum = dr::sum(majorant);

            ScalarFloat cost_iter = 2.f * (dim_prod / res_prod) * grid_sum;

            ScalarFloat cost_part = (resolution.x() - 1) * dim.y() * dim.z() +
                                    (resolution.y() - 1) * dim.x() * dim.z() +
                                    (resolution.z() - 1) * dim.x() * dim.y();

            // Timing parameters:
            // t_iter: cost incurred by having to perform a null collision loop
            // t_rewind: cost incurred by having to perform a dda loop
            // Dependent on optimization and target build, left for further
            // optimization.
            ScalarFloat t_iter   = 1.f;
            ScalarFloat t_rewind = 1.f;
            return (t_iter * cost_iter + t_rewind * cost_part) / denom;
        };

        auto to_scalar = [](Float f){
            if constexpr (dr::is_jit_v<Float>)
                return f[0];
            else
                return f;
        };

        // Ternary search over the xyz dimensions to find the best resolution.
        result = dr::maximum(fine_res / 2, 1);
        for (uint8_t i = 0; i < 3; ++i) {
            ScalarUInt32 left  = 1;
            ScalarUInt32 right = fine_res[i];

            ScalarFloat left_cost  = 0.f;
            ScalarFloat right_cost = 0.f;

            while ((right - left) > 3) {
                // Choose new boundaries
                ScalarUInt32 m1 = left + (right - left) / 3;
                ScalarUInt32 m2 = right - (right - left) / 3;

                ScalarVector3i left_res = result, right_res = result;
                left_res[i]  = m1;
                right_res[i] = m2;

                // Evaluate cost at the chosen boundaries
                build_grid(volume, left_res);
                left_cost = to_scalar(cost_fn(m_extremum_grid, left_res));

                build_grid(volume, right_res);
                right_cost = to_scalar(cost_fn(m_extremum_grid, right_res));

                if (left_cost < right_cost)
                    right = m2 - 1;
                else
                    left = m1 + 1;
            }

            result[i] = left_cost <= right_cost ? left : right;
        }

        // safety net, result should never be outside those bounds.
        result = dr::clip(result, 1, m_resolution - 1);
        return result;
    }

    /// World-space query box of a grid cell. Border cells are extended
    /// outward to the domain edge so that their bounds remain valid over
    /// everything outside the grid that a clamped lookup maps to them.
    template <typename Vector3iT, typename PointT>
    auto cell_bounds(const Vector3iT &ijk, const PointT &cell_size,
                     const ScalarVector3i &resolution) const {
        PointT cell_min = PointT(ijk) * cell_size;
        PointT cell_max = cell_min + cell_size;

        // Small inward shrink avoids bleeding into adjacent data texels at
        // exact cell boundaries
        cell_min += math::RayEpsilon<Float>;
        cell_max -= math::RayEpsilon<Float>;

        PointT wmin = PointT(m_extent.min) + cell_min * PointT(m_extent.extents());
        PointT wmax = PointT(m_extent.min) + cell_max * PointT(m_extent.extents());

        for (size_t i = 0; i < 3; ++i) {
            wmin[i] = dr::select(ijk[i] == 0, m_bbox.min[i], wmin[i]);
            wmax[i] = dr::select(ijk[i] == (uint32_t) (resolution[i] - 1),
                                 m_bbox.max[i], wmax[i]);
        }
        return BoundingBox<PointT>(wmin, wmax);
    }

    /**
     * \brief Build the extremum grid from a volume
     *
     * This method constructs a lower-resolution grid where each cell stores
     * the extinction bounds over the corresponding world-space region,
     * queried directly from the volume.
     */
    void build_grid(const Volume *volume, ScalarVector3i resolution) {

        // local space supergrid cell size
        const ScalarVector3f cell_size = dr::rcp(ScalarVector3f(resolution));

        ScalarVector2f safety_factor(1.f - dr::Epsilon<Float>,
                                     1.f + dr::Epsilon<Float>);

        // Allocate extremum grid data
        size_t n = dr::prod(resolution);

        size_t n_threads  = pool_size() + 1;
        size_t grain_size = std::max(n / (4 * n_threads), (size_t) 1);

        m_extremum_grid = dr::empty<FloatStorage>(n * 2);

        // Early return if using the global majorant.
        if (n == 1) {
            ScalarFloat max = volume->max();
            dr::scatter(m_extremum_grid,
                        m_scale * Vector2f(0.f, max) * safety_factor,
                        UInt32(0));
            return;
        }

        if constexpr (!dr::is_jit_v<Float>) {
            auto guard = volume->pin();

            dr::parallel_for(
                dr::blocked_range<size_t>(0, n, grain_size),
                [&](const dr::blocked_range<size_t> &range) {
                    for (auto idx = range.begin(); idx != range.end(); ++idx) {
                        // Store in linear array (Z-slowest, X-fastest)
                        uint32_t x = (uint32_t) (idx % resolution.x());
                        uint32_t y = (uint32_t) ((idx / resolution.x()) % resolution.y());
                        uint32_t z = (uint32_t) (idx / (resolution.x() * resolution.y()));

                        auto bounds = cell_bounds(ScalarVector3u(x, y, z),
                                                  ScalarPoint3f(cell_size),
                                                  resolution);

                        auto [min, maj] = volume->extremum(bounds);

                        dr::scatter(m_extremum_grid,
                                    m_scale * Vector2f(min, maj) *
                                        safety_factor,
                                    UInt32(idx));
                    }
                });
        } else {

            UInt32 idx = dr::arange<UInt32>((uint32_t) n);

            UInt32 x = idx % resolution.x() ;
            UInt32 y = (idx / resolution.x())  % resolution.y();
            UInt32 z =  idx / (resolution.x() * resolution.y());

            auto bounds = cell_bounds(Vector3u(x, y, z), Point3f(cell_size),
                                      resolution);

            auto [min, maj] = volume->extremum(bounds);

            dr::scatter(m_extremum_grid, m_scale * min * safety_factor.x(), idx*2);
            dr::scatter(m_extremum_grid, m_scale * maj * safety_factor.y(), idx*2+1);
            dr::sync_thread();
        }

        Log(Info, "Extremum grid constructed successfully");
    }

    /** \brief General regular grid DDA traversal algorithm.
     *
     * This method traverses the regular grid along the provided ray using the
     * DDA algorithm. At each traversed cell, it calls the function ``func``,
     * that can perform actions on the passed ``segment``, ``state``,
     * ``advance``, and ``active``. This function can be used for sampling
     * segments and perform various tracking algorithms.
     *
     * \param func  Function to be called at each step of the traversal. Must have
     *              the signature:
     *              (ExtremumSegment& segment,
     *               StateD& state,
     *               Mask& advance,
     *               Mask active) -> StateD
     *
     *              Changes to the segment, state, and condition can be done in place.
     * \param state The payload passed to ``func``.
     * \param ray   The ray along which the structure is traversed.
     * \param mint  The minimum distance along the ray.
     * \param maxt  The maximum distance along the ray.
     * \param active
     *
     * \return
     *      Returns the final state at the end of the traversal.
     */
    template<typename FuncT, typename StateT>
    std::decay_t<StateT> traverse_dda(
        FuncT&& func,
        StateT&& state,
        const Ray3f& ray,
        Float mint,
        Float maxt,
        UInt32 channel,
        Mask active
    ) const {
        using StateD = std::decay_t<StateT>;

        ExtremumSegment segment = dr::zeros<ExtremumSegment>();

        // Transform ray to local grid coordinates [0,res]³. Positions beyond
        // the grid extent look up the (extended) border cells via clamping.
        Vector3f res = Vector3f(m_resolution);
        Ray3f local_ray((m_to_local * ray.o) * res, // Normalize origin
                        (m_to_local * ray.d) * res, // Normalize direction
                        ray.time, ray.wavelengths);
        Vector3f rcp_d = dr::rcp(local_ray.d);
        auto inf_t     = local_ray.d == 0.f;
        auto d_pos     = local_ray.d >= 0.f;

        Float t_min = mint;
        Float t_max = maxt;

        active &= t_max > t_min && dr::isfinite(t_max);

        // Advance the ray to the start of the interval
        local_ray.o = dr::fmadd(local_ray.d, t_min, local_ray.o);
        t_max       = t_max - t_min;
        t_min       = 0.f;

        // Compute the integer step direction
        Vector3i step      = dr::select(d_pos, 1, -1);
        Vector3f abs_rcp_d = abs(rcp_d);

        // #NOFIX: Possibly overflow in case the number of cells in the domain
        // is larger than the range of int32. This check prevents this from
        // happening, future iterations will probably change the way periodic bounds
        // will be implemented.

        // Integer grid coordinates
        Vector3i pi = dr::floor2int<Vector3i>(local_ray.o);

        // Fractional entry position
        Vector3f p0 = local_ray.o - Vector3f(pi);
        // Step size to next interaction
        Vector3f dt_v =
            dr::select(d_pos, dr::fmadd(-p0, rcp_d, rcp_d), -p0 * rcp_d);
        dr::masked(dt_v, inf_t) = dr::Infinity<Float>;

        struct LoopState {
            ExtremumSegment segment;
            StateD state;
            Mask advance;
            Mask active;
            Vector3f dt_v;
            Vector3i pi;
            Float t_rem;

            DRJIT_STRUCT(LoopState, segment, state, advance, active, dt_v, \
                pi, t_rem)
        } ls = {
            segment,
            state,
            /*advance=*/active,
            active,
            dt_v,
            pi,
            t_max
        };

        dr::tie(ls) = dr::while_loop(
            dr::make_tuple(ls),
            [](const LoopState& ls) { return ls.active; },
            [this, func, step, abs_rcp_d, t_max, mint, channel](LoopState& ls) {

            ExtremumSegment& segment = ls.segment;
            StateD& state = ls.state;
            Mask& advance    = ls.advance;
            Mask& active    = ls.active;
            Vector3f& dt_v  = ls.dt_v;
            Vector3i& pi    = ls.pi;
            Float& t_rem = ls.t_rem;

            // Select the smallest step. It's possible that dt == 0 when
            // starting directly on a grid line.
            Float dt  = dr::minimum(dr::min(dt_v), t_rem);
            auto mask = dt_v == dt;

            // Note: not multiplying the index by 2 because we gather using
            // Vector2f.
            Vector3i piw = dr::clip(pi, 0, m_resolution - 1);
            UInt32 idx = dr::fmadd(dr::fmadd(piw.z(), m_resolution.y(), piw.y()),
                                   m_resolution.x(), piw.x());

            const Vector2f extremum    = dr::gather<Vector2f>(m_extremum_grid, idx);

            // Store segment for lanes that reached target
            Float t_curr = mint + t_max - t_rem;
            dr::masked(segment, active) = ExtremumSegment(
                t_curr,
                t_curr + dt,
                extremum
            );

            std::tie(advance, active) = func( segment, state, channel, active);

            // Advance
            dr::masked(dt_v, advance) = dr::select(mask, abs_rcp_d, dt_v - dt);
            dr::masked(pi, mask && advance) += step;
            dr::masked(t_rem, advance) -= dt;

            active &= t_rem > 0.f;
        },
        "DDA Traversal");

        return ls.state;
    };

private:
    /// Grid storing pre-computed local majorants
    FloatStorage m_extremum_grid;

    /// World-space region covered by the grid (domain ∩ volume bbox)
    ScalarBoundingBox3f m_extent;
    ScalarFloat m_scale = 1.f;

    bool m_adaptive;
    ScalarVector3i m_resolution;

    ScalarAffineTransform4f m_to_local;
};

MI_EXPORT_PLUGIN(ExtremumGrid)
NAMESPACE_END(mitsuba)
