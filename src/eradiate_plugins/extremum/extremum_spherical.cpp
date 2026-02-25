#include <mitsuba/core/properties.h>
#include <mitsuba/core/plugin.h>
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

 * - volume
   - |volume|
   - Spherical-coordinates volume to build extremum from
   - |exposed|

 * - to_world
   - |transform|
   - Specifies an optional 4x4 transformation matrix that will be applied to volume coordinates.

 * - rmin
   - |float|
   - Inner shell radius (Default: 0)

 * - rmax
   - |float|
   - Outer shell radius (Default: 1)

 * - resolution
   - |vector|
   - Grid resolution as (res_r, res_theta, res_phi) (Default: (1, 1, 1))

 * - scale
   - |float|
   - Scale factor for extinction coefficients (Default: 1.0)

This plugin creates a spherical extremum structure storing local majorant (and
minorant) values for efficient delta tracking in spherical media. The grid is
constructed by querying the extinction volume's extrema over each spherical cell.

At runtime, concentric shell traversal provides tight-fitting local majorants for
radially-varying media such as planetary atmospheres.
*/

template <typename Float, typename Spectrum>
class ExtremumSpherical final : public ExtremumStructure<Float, Spectrum> {
public:
    MI_IMPORT_BASE(ExtremumStructure, m_bbox)
    MI_IMPORT_TYPES(Volume)

    using Segment = ExtremumSegment;

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

        // Mark all properties as queried so they don't warn in expand()
        props.mark_queried("volume");
        props.mark_queried("to_world");
        props.mark_queried("rmin");
        props.mark_queried("rmax");
        props.mark_queried("resolution");
        props.mark_queried("scale");
        props.mark_queried("sample_method");
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
    Segment sample_segment(const Ray3f &, Float, Float, Float,
                           Mask) const override {
        NotImplementedError("sample_segment");
        return dr::zeros<Segment>();
    }

    std::tuple<Float, Float> eval_1(const Interaction3f &,
                                    Mask) const override {
        NotImplementedError("eval_1");
        return { 0.f, 0.f };
    }

    MI_DECLARE_CLASS(ExtremumSpherical)

protected:
    Properties m_props;
    ScalarPoint3f m_center;
    SphericalTraversalType m_traversal_type;
};


// ---------------------------------------------------------------------------
// Implementation class
// ---------------------------------------------------------------------------

template <typename Float, typename Spectrum, SphericalTraversalType TraversalType>
class ExtremumSphericalImpl final : public ExtremumStructure<Float, Spectrum> {
public:
    MI_IMPORT_BASE(ExtremumStructure, m_bbox)
    MI_IMPORT_TYPES(Volume)

    using Segment = ExtremumSegment;
    using FloatStorage = DynamicBuffer<Float>;

    ExtremumSphericalImpl(const Properties &props) : Base(props) {
        // Get volume
        ref<Volume> volume = nullptr;
        for (auto &prop : props.objects()) {
            if (auto *vol = prop.try_get<Volume>()) {
                volume = vol;
                break;
            }
        }

        if (!volume)
            Throw("ExtremumSpherical requires a volume");

        m_volume = volume;
        m_bbox = m_volume->bbox();
        // Register the extremum structure to the volume.
        volume->add_extremum_structure(this);

        ScalarAffineTransform4f to_world = props.get<ScalarAffineTransform4f>("to_world", ScalarAffineTransform4f());
        m_center = to_world.translation();
        m_to_local = to_world.inverse(); 

        m_rmin = props.get<ScalarFloat>("rmin", 0.f);
        m_rmax = props.get<ScalarFloat>("rmax", 1.f);
        m_resolution = props.get<ScalarVector3i>("resolution", ScalarVector3i(1, 1, 1));
        m_scale = props.get<ScalarFloat>("scale", 1.0f);
        
        if (m_rmin >= m_rmax)
            Throw("rmin must be less than rmax!");

        m_dr = (m_rmax - m_rmin) / m_resolution.x();
        m_idr = dr::rcp(m_dr);

        build_grid(m_volume.get());

        Log(Info, "ExtremumSpherical created: resolution=%s, center=%s, "
            "rmin=%f, rmax=%f", m_resolution, m_center, m_rmin, m_rmax);
    }

    void parameters_changed(const std::vector<std::string> &/*keys*/ = {}) override {
        build_grid(m_volume.get());
    }

    Segment sample_segment(
        const Ray3f &ray,
        Float mint, Float maxt,
        Float desired_tau,
        Mask active
    ) const override {

        if constexpr (TraversalType == SphericalTraversalType::RadialOnly) {
            return sample_segment_radial(ray, mint, maxt, desired_tau, active);
            
        } else {
            Throw("Full3D spherical traversal is not yet implemented!");
            return dr::zeros<Segment>();
        }
    }

    std::tuple<Float, Float> eval_1(
        const Interaction3f &it,
        Mask active
    ) const override {
        Float r = dr::norm(it.p - m_center);
        UInt32 ir = dr::clip(
            dr::floor2int<UInt32>((r - m_rmin) / m_dr),
            0u, (uint32_t)(m_resolution.x() - 1)
        );
        Vector2f extremum = dr::gather<Vector2f>(m_extremum_grid, ir, active);
        return { extremum.x(), extremum.y() };
    }

    void traverse(TraversalCallback *cb) override {
        cb->put("data", m_extremum_grid, ParamFlags::NonDifferentiable);
        cb->put("resolution", m_resolution, ParamFlags::NonDifferentiable);
        cb->put("scale", m_scale, ParamFlags::NonDifferentiable);
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
            << "  to_world = "   << m_to_local.inverse() << "," << std::endl
            << "  rmin = "       << m_rmin << "," << std::endl
            << "  rmax = "       << m_rmax << std::endl
            << "  fillmin = "    << m_fillmin << "," << std::endl
            << "  fillmax = "    << m_fillmax << std::endl
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

        ScalarVector2f safety_factor(0.99, 1.01);

        std::vector<ScalarFloat> extremums;

        // Retrieve data pointer before parallel_for to avoid ref count issues
        auto data = volume->array();

        size_t n_threads = pool_size() + 1;
        size_t grain_size = std::max(n / (4 * n_threads), (size_t) 1);

        if constexpr (!dr::is_jit_v<Float>) {

            m_extremum_grid = dr::empty<FloatStorage>(n * 2);

            dr::parallel_for(
                dr::blocked_range<size_t>(0, n, grain_size),
                [&](const dr::blocked_range<size_t> &range) {
                    for (auto idx = range.begin(); idx != range.end(); ++idx) {
                        // r-fastest indexing: idx = ir + itheta * res_r + iphi * res_r * res_theta
                        int32_t ir     = idx % m_resolution.x();
                        int32_t itheta = (idx / m_resolution.x()) % m_resolution.y();
                        int32_t iphi   = idx / (m_resolution.x() * m_resolution.y());

                        ScalarPoint3f cell_min = ScalarVector3f(ir, itheta, iphi) * cell_size;
                        ScalarPoint3f cell_max = cell_min + cell_size;
                        ScalarBoundingBox3f cell_bounds(
                            cell_min + math::RayEpsilon<Float>,
                            cell_max - math::RayEpsilon<Float>
                        );

                        auto [maj, min] = volume->extremum(data, cell_bounds);

                        dr::scatter(m_extremum_grid,
                                    m_scale * Vector2f(min, maj) * safety_factor,
                                    UInt32(idx));
                    }
                }
            );
        } else {

            m_extremum_grid = dr::empty<FloatStorage>(n * 2);

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

            auto [maj, min] = volume->extremum(data, cell_bounds);

            dr::scatter(m_extremum_grid, m_scale * min * safety_factor.x(), idx * 2);
            dr::scatter(m_extremum_grid, m_scale * maj * safety_factor.y(), idx * 2 + 1);
            dr::sync_thread();
        }

        // Retrieve fillmin and fillmax
        Interaction3f it = dr::zeros<Interaction3f>();

        it.p = m_center;
        Float fillmin = volume->eval_1(it, true) * m_scale;
        
        it.p = m_center + m_rmax + 1;
        Float fillmax = volume->eval_1(it, true) * m_scale;

        if constexpr (dr::is_jit_v<Float>) {
            m_fillmin = fillmin[0];
            m_fillmax = fillmax[1];
        } else {
            m_fillmin = fillmin;
            m_fillmax = fillmax;
        }

        Log(Info, "Extremum spherical grid constructed successfully");
    }

    // ------------------------------------------------------------------
    // RadialOnly shell traversal
    // ------------------------------------------------------------------
    Segment sample_segment_radial(
        const Ray3f &ray,
        Float mint, Float maxt,
        Float desired_tau,
        Mask active
    ) const {

        // const ScalarVector3f extents = m_bbox.extents();
        Ray3f local_ray(
            m_to_local * ray.o,  // Normalize origin
            m_to_local * ray.d,  // Normalize direction
            ray.time,
            ray.wavelengths
        );

        Segment result = dr::zeros<Segment>();
        Float tau_acc = 0.f;
        Mask reached = false;
        Float current_t = mint;

        // Log(Debug, "m_rmin: %f, m_rmax: %f, m_dr: %f, desired_tau: %f", m_rmin, m_rmax, m_dr, desired_tau);
        // Log(Debug, "Start Loop");

        // ray-sphere intersection info
        Vector3f o = local_ray.o - m_center;
        Float o_squared = dr::squared_norm(o);
        Float a = dr::squared_norm(local_ray.d);
        Float b_half = dr::dot(o, local_ray.d);
        
        // claude addition
        Float disc_base = b_half * b_half - a * o_squared;  // constant across iterations
        Float inv_a = dr::rcp(a); 

        // Find the current/next intersection (use this to calculate the midpoint too)
        Point3f pos = local_ray(mint+dr::Epsilon<Float>*10.f);
        Vector3f oc = pos - m_center;
        Float r = dr::norm(oc);

        // Calculate the initial layer index from which we will step through layers.
        Int32 layer_idx = dr::clip(dr::floor2int<Int32>((r - m_rmin) * m_idr), -1, m_resolution.x());
        Mask passed_midpoint = dr::dot((m_center - pos), local_ray.d) < 0;
        Int32 shell_padding = dr::select(passed_midpoint, 1, 0);
        Int32 step = dr::select(passed_midpoint, 1, -1);

        // Log(Debug, "r: %f, (r-m_rmin)*m_idr: %f, layer_idx: %f", r, dr::floor2int<Int32>((r - m_rmin) * m_idr), layer_idx);

        struct LoopState {
            Segment result;
            Mask active, reached;
            Float current_t, tau_acc;
            Int32 layer_idx;
            Int32 step;
            Int32 padding;

            DRJIT_STRUCT(                                                      \
                LoopState, result, active, reached, current_t,                 \
                tau_acc, layer_idx, step, padding                              \
            )
        } ls = { 
            result, 
            active, 
            reached, 
            current_t, 
            tau_acc, 
            layer_idx, 
            step,
            shell_padding, 
        };

        dr::tie(ls) = dr::while_loop(
            dr::make_tuple(ls),
            [](const LoopState &ls) { return ls.active; },
            // [this, maxt, desired_tau, o_squared, a, b_half](LoopState &ls) {
            [this, maxt, desired_tau, a, inv_a, disc_base, b_half](LoopState &ls) {
            // Log(Debug, "---------");

            // Compute radius at current position
            const Float eps = math::RayEpsilon<Float>;
            
            // Passed midpoint == exiting the concentric spheres
            const Int32 shell_idx = dr::clip(ls.layer_idx + ls.padding, 0, m_resolution.x());
            // Log(Debug, "current_t: %f, layer_idx: %f, shell_idx: %f", current_t, layer_idx, shell_idx);
            // Log(Debug, "padding: %f , step: %f, what: %f", padding, step, layer_idx + padding);

            // Boundary condition of rmin and rmax
            Float fill_value = -1.f;
            dr::masked(fill_value, ls.layer_idx < 0) = m_fillmin;
            dr::masked(fill_value, ls.layer_idx >= m_resolution.x()) = m_fillmax;
            const Mask fill = fill_value >= 0.f;

            // Test intersection with the shell
            const Float r_test = m_rmin + Float(shell_idx) * m_dr;
            
            // Float c = o_squared - dr::square(r_test);
            // auto [valid_test, t_test_near, t_test_far] = math::solve_quadratic(a, 2.f * b_half, c);
            
            const Float disc = disc_base + a * dr::square(r_test);
            const Mask valid_test = disc >= 0.f;
            const Float sqrt_disc = dr::sqrt(disc);
            const Float t_test_near = (-b_half - sqrt_disc) * inv_a;
            const Float t_test_far  = (-b_half + sqrt_disc) * inv_a;

            // Log(Debug, "r_test: %f, valid_test: %f", r_test, valid_test);
            // Log(Debug, "far: %f, near: %f", t_test_far, t_test_near);
            // Log(Debug, "a, b_half, c_lo, c_hi: %f, %f, %f, %f", a, b_half, c);
            
            // Update if there is valid intersection that is not tangent.
            Mask update = valid_test && dr::abs(t_test_far - t_test_near) > math::RayEpsilon<Float>;
            const Mask pass_midpoint = !update || ls.layer_idx == -1;

            // Mask skip_update = !valid_test || dr::abs(t_test_far - t_test_near) < math::RayEpsilon<Float>;
            // Mask pass_midpoint = skip_update || layer_idx == -1;
            
            // Special case at midpoint where we miss the intersection with the 
            // shell, bump padding to test outer shell, set step to increase
            // layer index and continue to next iteration.
            dr::masked(ls.padding, pass_midpoint) = 1;
            dr::masked(ls.step, pass_midpoint) = 1;

            // Mask update = active && valid_test && !tangent; 
            // Log(Debug, "update: %f, pass_midpoint: %f", update, pass_midpoint);

            if( dr::any_or<true>(update) ) {
            
                // Find smallest t > current_t + epsilon among the 4 candidates
                // Float threshold = current_t + eps;
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
                Vector2f local_extremum = dr::gather<Vector2f>(
                    m_extremum_grid, ls.layer_idx, ls.active && !fill);
                Float local_minorant = dr::select(!fill, local_extremum.x(), fill_value);
                Float local_majorant = dr::select(!fill, local_extremum.y(), fill_value);

                // Accumulate optical depth
                const Float segment_length = t_next - ls.current_t;
                const Float tau_next = dr::fmadd(local_majorant, segment_length, ls.tau_acc);
                
                // Check if desired tau is reached in this segment
                const Mask stops_here = ls.active && (local_majorant > 0.f) &&
                                    (tau_next >= desired_tau) && (t_next <= maxt);

                ls.reached |= stops_here;
                // Log(Debug, "majorant: %f, minorant: %f", local_majorant, local_minorant);
                // Log(Debug, "tau_next: %f, t_next: %f", tau_next, t_next);
                // Log(Debug, "stops_here: %f, reached: %f", stops_here, reached);
                // Store result for lanes that reached target
                dr::masked(ls.result, stops_here) = ExtremumSegment(
                    ls.current_t, t_next, local_majorant, local_minorant, ls.tau_acc
                );

                // Advance state for lanes that haven't reached target
                update &= !ls.reached;
                dr::masked(ls.current_t, update) = t_next;
                dr::masked(ls.tau_acc, update) = tau_next;
                dr::masked(ls.layer_idx, update) += ls.step;

                // Log(Debug, "current_t: %f, layer_idx: %f", current_t, layer_idx);
                // Log(Debug, "count: %f", ls.count);

                // Continue only if not reached and still in bounds
                ls.active &= (!ls.reached && (t_next < maxt)) ;
            }
        },
        "Spherical Shell Traversal");

        // For lanes that didn't reach, invalidate the segment
        dr::masked(ls.result.tmin, !ls.reached) = dr::Infinity<Float>;

        return ls.result;
    }

private:
    FloatStorage m_extremum_grid;
    ScalarVector3i m_resolution;
    ScalarFloat m_rmin, m_rmax;
    ScalarFloat m_fillmin, m_fillmax;
    ScalarPoint3f m_center;
    ScalarFloat m_dr, m_idr;

    ref<Volume> m_volume;

    ScalarFloat m_scale;
    ScalarAffineTransform4f m_to_local;
    ScalarUInt32 m_sample_method;
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
