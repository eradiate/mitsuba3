#include <mitsuba/core/properties.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/eradiate/extremum.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/volumegrid.h>
#include <nanothread/nanothread.h>

NAMESPACE_BEGIN(mitsuba)

enum class SphericalTraversalType { RadialOnly, Full3D };

template <typename Float, typename Spectrum, SphericalTraversalType TraversalType>
class ExtremumSphericalImpl;

/**!
.. _extremum-extremum_spherical:

Extremum spherical structure (:monosp:`extremum_spherical`)
-----------------------------------------------------------

.. pluginparameters::

 * - resolution
   - |vector|
   - Grid resolution as :math:`(r, \theta, \phi)`. Grids with variations on only
     the radial resolution have optimized traversal. Default: [1,1,1]

This plugin creates a spherical extremum structure storing local extremum values
for efficient delta tracking in spherical media. All geometric information
(domain, extinction volume, scale) is provided by the owning medium through
\ref ExtremumStructure::build(); the shell center and radial extent are derived
from the volume's spherical frame. The grid is constructed by querying the
underlying volume's extrema over each spherical cell.

At runtime, concentric shell traversal provides tight-fitting local extremum for
radially-varying media such as planetary atmospheres.

.. warning:: This extremum only implements radial extremum variations.
*/

template <typename Float, typename Spectrum>
class ExtremumSpherical final : public ExtremumStructure<Float, Spectrum> {
public:
    MI_IMPORT_BASE(ExtremumStructure, m_bbox)
    MI_IMPORT_TYPES(Volume)

    using TrackingStateType    = TrackingState<Float, Spectrum>;
    using TrackingFunctionType = TrackingFunction<Float, Spectrum>;

    ExtremumSpherical(const Properties &props) : Base(props), m_props(props) {
        ScalarVector3i resolution = props.get<ScalarVector3i>("resolution", ScalarVector3i(1, 1, 1));

        if (resolution.x() < 1 || resolution.y() < 1 || resolution.z() < 1)
            Throw("All resolution components must be >= 1!");

        // Determine traversal type from resolution
        if (resolution.y() == 1 && resolution.z() == 1) {
            m_traversal_type = SphericalTraversalType::RadialOnly;
        } else {
            m_traversal_type = SphericalTraversalType::Full3D;
        }
    }

    template <SphericalTraversalType TT>
    using Impl = ExtremumSphericalImpl<Float, Spectrum, TT>;

    std::vector<ref<Object>> expand() const override {
        ref<Object> result;
        switch (m_traversal_type) {
            case SphericalTraversalType::RadialOnly:
                result = (Object *) new Impl<SphericalTraversalType::RadialOnly>(m_props);
                break;
            case SphericalTraversalType::Full3D:
                Throw("Full3D spherical traversal is not yet implemented!");
                break;
            default:
                Throw("Unsupported spherical traversal type!");
        }
        return { result };
    }

    // Stub overrides — never called, expand() replaces this object
    TrackingStateType traverse_extremum(const Ray3f &, Float, Float, UInt32,
                    TrackingStateType, TrackingFunctionType*, Mask) const override {
        NotImplementedError("traverse_extremum");
    }

    std::tuple<Float, Float> eval_1(const Interaction3f &,
                                    Mask) const override {
        NotImplementedError("eval_1");
        return { 0.f, 0.f };
    }

    MI_DECLARE_CLASS(ExtremumSpherical)

protected:
    Properties m_props;
    SphericalTraversalType m_traversal_type;
};


// ---------------------------------------------------------------------------
// Implementation class
// ---------------------------------------------------------------------------

template <typename Float, typename Spectrum, SphericalTraversalType TraversalType>
class ExtremumSphericalImpl final : public ExtremumStructure<Float, Spectrum> {
public:
    MI_IMPORT_BASE(ExtremumStructure, m_bbox, set_domain)
    MI_IMPORT_TYPES(Volume)

    using TrackingStateType    = TrackingState<Float, Spectrum>;
    using TrackingFunctionType = TrackingFunction<Float, Spectrum>;
    using FloatStorage         = DynamicBuffer<Float>;

    ExtremumSphericalImpl(const Properties &props) : Base(props) {
        m_resolution =
            props.get<ScalarVector3i>("resolution", ScalarVector3i(1, 1, 1));
    }

    void build(const ScalarBoundingBox3f &domain, const Volume *volume,
               ScalarFloat scale) override {
        set_domain(domain, volume);
        m_scale = scale;

        // Correctness gate, not a tightness hint: the axis meanings of the
        // extremum_local queries below depend on the volume's native
        // parameterization being (r, theta, phi).
        SphericalParameters<ScalarFloat> frame = volume->spherical_frame();
        if (!frame.valid())
            Throw("extremum_spherical: the volume volume does not expose a "
                  "spherical frame (e.g. sphericalcoordsvolume); its "
                  "extremum_local axes would be misinterpreted!");

        m_center = frame.center;
        m_rmin   = frame.rmin;
        m_rmax   = frame.rmax;

        if (m_rmin >= m_rmax)
            Throw("rmin must be less than rmax!");

        m_dr  = (m_rmax - m_rmin) / m_resolution.x();
        m_idr = dr::rcp(m_dr);

        build_grid(volume);
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
        if constexpr (TraversalType == SphericalTraversalType::RadialOnly) {
            return traverse_radial(
                func,
                state,
                ray,
                mint,
                maxt,
                channel,
                active
            );
        } else {
            Throw("Full3D spherical traversal is not yet implemented!");
            return TrackingStateType();
        }
    }

    ExtremumSegment next_segment(const Ray3f &ray, Float t,
                                 Mask active) const override {
        auto [hit, d0, d1] = m_bbox.ray_intersect(ray);

        // Forward nudge: select the radial layer at p(t + eps), but measure
        // exit distances from exact p(t)
        Float eps = dr::maximum(dr::abs(t), 1.f) * math::RayEpsilon<Float>;
        Float tq  = t + eps;

        Mask before = hit && (tq < d0);
        Mask inside = hit && !before && (tq < d1);

        Float r = dr::norm(ray(tq) - m_center);
        Int32 layer = dr::clip(dr::floor2int<Int32>((r - m_rmin) * m_idr),
                               -1, m_resolution.x());
        Mask below = layer < 0, above = layer >= m_resolution.x();

        // Nearest crossing (> tq) of the two shells bounding the current
        // radial region; the domain exit caps the segment
        Vector3f oc  = ray.o - m_center;
        Float a      = dr::squared_norm(ray.d);
        Float b_half = dr::dot(oc, ray.d);
        Float c0     = dr::squared_norm(oc);
        Float inv_a  = dr::rcp(a);

        Float t_exit = d1;
        auto consider = [&](const Float &r_test, const Mask &valid)
            DRJIT_INLINE_LAMBDA {
            Float disc = dr::fmsub(b_half, b_half, a * (c0 - dr::square(r_test)));
            Mask v = valid && (disc >= 0.f);
            Float sq = dr::sqrt(dr::maximum(disc, 0.f));
            Float t0 = (-b_half - sq) * inv_a;
            Float t1 = (-b_half + sq) * inv_a;
            dr::masked(t_exit, v && (t0 > tq) && (t0 < t_exit)) = t0;
            dr::masked(t_exit, v && (t1 > tq) && (t1 < t_exit)) = t1;
        };
        consider(m_rmin + Float(layer) * m_dr, !below);      // inner shell
        consider(m_rmin + Float(layer + 1) * m_dr, !above);  // outer shell

        UInt32 gidx = UInt32(dr::clip(layer, 0, m_resolution.x() - 1));
        Vector2f value = dr::gather<Vector2f>(
            m_extremum_grid, gidx, active && inside && !below && !above);
        dr::masked(value, below) = Vector2f(m_fill_inner);
        dr::masked(value, above) = Vector2f(m_fill_outer);

        Float maxt = dr::select(
            inside, dr::maximum(t_exit, tq),
            dr::select(before, d0, dr::Infinity<Float>));

        return ExtremumSegment(t, maxt,
                               dr::select(inside, value, Vector2f(0.f)));
    }

    std::tuple<Float, Float> eval_1(
        const Interaction3f &it,
        Mask active
    ) const override {
        Float r = dr::norm(it.p - m_center);
        Mask below = r < m_rmin, above = r > m_rmax;

        Vector2f extremum = dr::zeros<Vector2f>();

        if constexpr (TraversalType == SphericalTraversalType::RadialOnly) {
            // Note that this is only valid for radial only for now
            UInt32 ir = dr::clip(
                dr::floor2int<UInt32>((r - m_rmin) * m_idr),
                0u, (uint32_t)(m_resolution.x() - 1)
            );
            extremum = dr::gather<Vector2f>(m_extremum_grid, ir,
                                            active && !below && !above);
        } else if constexpr (TraversalType == SphericalTraversalType::Full3D) {
            Throw("Full3D spherical evaluation is not yet implemented!");
        }

        dr::masked(extremum, below) = Vector2f(m_fill_inner);
        dr::masked(extremum, above) = Vector2f(m_fill_outer);

        return { extremum.x(), extremum.y() };
    }

    void traverse(TraversalCallback *cb) override {
        cb->put("data", m_extremum_grid, ParamFlags::NonDifferentiable);
        cb->put("resolution", m_resolution, ParamFlags::NonDifferentiable);
        Base::traverse(cb);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "ExtremumSpherical[" << std::endl
            << "  traversal = "
            << (TraversalType == SphericalTraversalType::RadialOnly
                ? "RadialOnly" : "Full3D")
            << "," << std::endl
            << "  resolution = " << m_resolution << "," << std::endl
            << "  center = "     << m_center << "," << std::endl
            << "  rmin = "       << m_rmin << "," << std::endl
            << "  rmax = "       << m_rmax << "," << std::endl
            << "  fill_inner = " << m_fill_inner << "," << std::endl
            << "  fill_outer = " << m_fill_outer << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS(ExtremumSphericalImpl)

private:

    // ------------------------------------------------------------------
    // Grid construction
    // ------------------------------------------------------------------

    void build_grid(const Volume *volume) {
        // Cell size in normalized [0,1]^3 space
        const ScalarVector3f cell_size = dr::rcp(ScalarVector3f(m_resolution));
        size_t n = dr::prod(m_resolution);

        ScalarVector2f safety_factor(
            1.f - dr::Epsilon<Float>,
            1.f + dr::Epsilon<Float>
        );

        size_t n_threads = pool_size() + 1;
        size_t grain_size = std::max(n / (4 * n_threads), (size_t) 1);

        m_extremum_grid = dr::empty<FloatStorage>(n * 2);

        if constexpr (!dr::is_jit_v<Float>) {
            auto guard = volume->pin();

            dr::parallel_for(
                dr::blocked_range<size_t>(0, n, grain_size),
                [&](const dr::blocked_range<size_t> &range) {
                    for (auto idx = range.begin(); idx != range.end(); ++idx) {
                        // r-fastest indexing: idx = ir + itheta * res_r + iphi * res_r * res_theta
                        int32_t ir     = idx % m_resolution.x();
                        int32_t itheta = (idx / m_resolution.x()) % m_resolution.y();
                        int32_t iphi   = idx / (m_resolution.x() * m_resolution.y());

                        ScalarPoint3f cell_min =
                            ScalarVector3f(ir, itheta, iphi) * cell_size;
                        ScalarPoint3f cell_max = cell_min + cell_size;
                        ScalarBoundingBox3f cell_bounds(
                            cell_min + math::RayEpsilon<Float>,
                            cell_max - math::RayEpsilon<Float>);

                        auto [min, maj] = volume->extremum_local(cell_bounds);

                        dr::scatter(m_extremum_grid,
                                    m_scale * Vector2f(min, maj) * safety_factor,
                                    UInt32(idx));
                    }
                }
            );
        } else {
            UInt32 idx = dr::arange<UInt32>((uint32_t) n);

            UInt32 ir     = idx % m_resolution.x();
            UInt32 itheta = (idx / m_resolution.x()) % m_resolution.y();
            UInt32 iphi   = idx / (m_resolution.x() * m_resolution.y());

            Point3f cell_min = Vector3f(ir, itheta, iphi) * cell_size;
            Point3f cell_max = cell_min + cell_size;
            BoundingBox3f cell_bounds(
                cell_min + math::RayEpsilon<Float>,
                cell_max - math::RayEpsilon<Float>
            );

            auto [min, maj] = volume->extremum_local(cell_bounds);

            dr::scatter(m_extremum_grid, m_scale * min * safety_factor.x(), idx * 2);
            dr::scatter(m_extremum_grid, m_scale * maj * safety_factor.y(), idx * 2 + 1);
            dr::sync_thread();
        }

        // Fill bounds: the radial analog of border-cell extension. The
        // regions r < rmin and r > rmax are covered by bounds queries over
        // the corresponding (out-of-range) radial intervals, which the
        // volume folds into its fill values.
        ScalarFloat rext = m_rmax - m_rmin;

        auto query_fill = [&](ScalarFloat x0, ScalarFloat x1) -> ScalarVector2f {
            auto [mn, mx] = volume->extremum_local(BoundingBox3f(
                Point3f(x0, 0.f, 0.f), Point3f(x1, 1.f, 1.f)));
            ScalarFloat mn_s, mx_s;
            if constexpr (dr::is_jit_v<Float>) {
                mn_s = mn[0]; mx_s = mx[0];
            } else {
                mn_s = mn; mx_s = mx;
            }
            return m_scale * ScalarVector2f(mn_s, mx_s) * safety_factor;
        };

        m_fill_inner = m_rmin > 0.f
            ? query_fill(-m_rmin / rext, 0.f)
            : ScalarVector2f(0.f);

        // Furthest domain point from the center bounds the outer region
        ScalarFloat r_domain = 0.f;
        for (uint32_t i = 0; i < 8; ++i)
            r_domain = dr::maximum(
                r_domain, dr::norm(m_bbox.corner(i) - m_center));
        m_fill_outer =
            query_fill(1.f, dr::maximum((r_domain - m_rmin) / rext, 1.f));

        Log(Info, "Extremum spherical grid constructed successfully");
    }

    // ------------------------------------------------------------------
    // RadialOnly shell traversal
    // ------------------------------------------------------------------
    /** \brief General radial traversal algorithm.
     *
     * This method traverses the regular concentric-shell structure along the
     * provided ray. At each traversed cell, it calls the function ``func``,
     * that can perform actions on the passed ``segment``, ``state``,``advance``,
     * and ``active``. This function can be used for sampling segments and
     * perform various tracking algorithms.
     *
     * \param func  Function to be called at each step of the traversal. Must have
     *              the signature:
     *              (ExtremumSegment& segment,
     *               StateD& state,
     *               Mask& condition,
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
    std::decay_t<StateT> traverse_radial(
        FuncT&& func,
        StateT&& state,
        const Ray3f ray,
        Float mint,
        Float maxt,
        UInt32 channel,
        Mask active
    ) const {
        using StateD = std::decay_t<StateT>;

        ExtremumSegment segment  = dr::zeros<ExtremumSegment>();
        Mask reached    = false;
        Float current_t = mint;

        // ray-sphere intersection info
        Vector3f o      = ray.o - m_center;
        Float o_squared = dr::squared_norm(o);
        Float a         = dr::squared_norm(ray.d);
        Float b_half    = dr::dot(o, ray.d);

        // Intersection value precomputation
        Float disc_base = b_half * b_half - a * o_squared;
        Float inv_a     = dr::rcp(a);

        // Find the current/next intersection (use this to calculate the
        // midpoint too)
        Point3f pos = ray(mint + dr::Epsilon<Float> * 10.f);
        Vector3f oc = pos - m_center;
        Float r     = dr::norm(oc);

        // Calculate the initial layer index from which we will step through
        // layers.
        Int32 layer_idx = dr::clip(dr::floor2int<Int32>((r - m_rmin) * m_idr),
                                   -1, m_resolution.x());
        Mask passed_midpoint = dr::dot((m_center - pos), ray.d) < 0;
        Int32 shell_padding  = dr::select(passed_midpoint, 1, 0);
        Int32 step           = dr::select(passed_midpoint, 1, -1);

        struct LoopState {
            ExtremumSegment segment;
            StateD state;
            Mask advance;
            Mask active;
            Mask reached;
            Float current_t;
            Int32 layer_idx;
            Int32 step;
            Int32 padding;

            DRJIT_STRUCT(                                                      \
                LoopState, segment, state, advance, active, reached, current_t,\
                 layer_idx, step, padding                                      \
            )
        } ls = {
            segment,
            state,
            /*advance=*/active,
            active,
            reached,
            current_t,
            layer_idx,
            step,
            shell_padding,
        };

        dr::tie(ls) = dr::while_loop(
            dr::make_tuple(ls),
            [](const LoopState &ls) { return ls.active; },
            [this, func, maxt, a, inv_a, disc_base, b_half, channel](LoopState &ls) {

            // Compute radius at current position
            const Float eps = dr::Epsilon<Float> * 2.f;

            // Passed midpoint == exiting the concentric spheres
            const Int32 shell_idx = dr::clip(ls.layer_idx + ls.padding, 0, m_resolution.x());

            // Boundary condition of rmin and rmax
            const Mask below = ls.layer_idx < 0;
            const Mask above = ls.layer_idx >= m_resolution.x();
            const Mask fill  = below || above;

            // Test intersection with the shell
            const Float r_test = m_rmin + Float(shell_idx) * m_dr;

            const Float disc = disc_base + a * dr::square(r_test);
            const Mask valid_test = disc >= 0.f;
            const Float sqrt_disc = dr::sqrt(disc);
            const Float t_test_near = (-b_half - sqrt_disc) * inv_a;
            const Float t_test_far  = (-b_half + sqrt_disc) * inv_a;

            // Update if there is valid intersection that is not tangent.
            Mask update = ls.active
                          && valid_test
                          && dr::abs(t_test_far - t_test_near) > eps;
            const Mask pass_midpoint = !update || ls.layer_idx == -1;

            // Special case at midpoint where we miss the intersection with the
            // shell, bump padding to test outer shell, set step to increase
            // layer index and continue to next iteration.
            dr::masked(ls.padding, pass_midpoint) = 1;
            dr::masked(ls.step, pass_midpoint)    = 1;

            if( dr::any_or<true>(update) ) {

                // Find smallest t > current_t + epsilon among the 4 candidates
                const Float threshold = ls.current_t + eps;
                Float t_next = maxt;

                // Helper: update t_next if candidate > threshold and < t_next
                auto consider = [&](Float t_cand, Mask valid_cand) DRJIT_INLINE_LAMBDA {
                    Mask use = valid_cand && (t_cand > threshold) && (t_cand < t_next);
                    dr::masked(t_next, use) = t_cand;
                };

                consider(t_test_near, valid_test);
                consider(t_test_far,  valid_test);

                // Look up extremum values for this shell
                Vector2f extremum = dr::gather<Vector2f>(
                    m_extremum_grid, ls.layer_idx, ls.active && !fill);
                dr::masked(extremum, below) = Vector2f(m_fill_inner);
                dr::masked(extremum, above) = Vector2f(m_fill_outer);

                dr::masked(ls.segment, update) = ExtremumSegment(
                    ls.current_t, t_next, extremum
                );

                Mask active_update = ls.active && update;
                auto result =
                    func( ls.segment, ls.state, channel, active_update);
                ls.advance    = result.first;
                active_update &= result.second;

                // Advance state for lanes that haven't reached target
                dr::masked(ls.current_t, ls.advance && update) = t_next;
                dr::masked(ls.layer_idx, ls.advance && update) += ls.step;

                // Continue only if not reached and still in bounds
                dr::masked(ls.active, update) &= active_update && (t_next <= maxt);
            }
            ls.active &= ls.layer_idx < m_resolution.x();
        },
        "Spherical Shell Traversal");

        return ls.state;
    }


private:
    FloatStorage m_extremum_grid;
    ScalarVector3i m_resolution;
    ScalarFloat m_rmin, m_rmax;
    ScalarVector2f m_fill_inner, m_fill_outer;
    ScalarPoint3f m_center;
    ScalarFloat m_dr, m_idr;

    ScalarFloat m_scale = 1.f;
};


// ---------------------------------------------------------------------------
// Class name helpers (for expand pattern)
// ---------------------------------------------------------------------------

MI_EXPORT_PLUGIN(ExtremumSpherical)

NAMESPACE_BEGIN(detail)
template <SphericalTraversalType TT>
constexpr const char *extremum_spherical_class_name() {
    if constexpr (TT == SphericalTraversalType::RadialOnly) {
        return "ExtremumSpherical_RadialOnly";
    } else if constexpr (TT == SphericalTraversalType::Full3D) {
        return "ExtremumSpherical_Full3D";
    }
}
NAMESPACE_END(detail)

NAMESPACE_END(mitsuba)
