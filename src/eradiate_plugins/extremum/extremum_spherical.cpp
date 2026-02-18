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

 * - sigma_t
   - |volume|
   - Spherical-coordinates volume to build extremum from
   - |exposed|

 * - center
   - |point|
   - Center of spherical coordinate system (Default: (0, 0, 0))

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
        m_center = props.get<ScalarPoint3f>("center", ScalarPoint3f(0.f));
        ScalarFloat rmin = props.get<ScalarFloat>("rmin", 0.f);
        ScalarFloat rmax = props.get<ScalarFloat>("rmax", 1.f);
        ScalarVector3i resolution = props.get<ScalarVector3i>("resolution", ScalarVector3i(1, 1, 1));

        if (rmin >= rmax)
            Throw("rmin must be less than rmax!");

        if (resolution.x() < 1 || resolution.y() < 1 || resolution.z() < 1)
            Throw("All resolution components must be >= 1!");

        // Determine traversal type from resolution
        if (resolution.y() == 1 && resolution.z() == 1) {
            m_traversal_type = SphericalTraversalType::RadialOnly;
        } else {
            m_traversal_type = SphericalTraversalType::Full3D;
        }

        // Mark all properties as queried so they don't warn in expand()
        props.mark_queried("sigma_t");
        props.mark_queried("center");
        props.mark_queried("rmin");
        props.mark_queried("rmax");
        props.mark_queried("resolution");
        props.mark_queried("scale");
        props.mark_queried("to_world");
        props.mark_queried("sample_method");

        // Set bounding box to enclose the outer sphere
        m_bbox = ScalarBoundingBox3f(
            m_center - ScalarVector3f(rmax),
            m_center + ScalarVector3f(rmax)
        );
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
        // m_volume = props.get_volume<Volume>("volume");
        for (auto &prop : props.objects()) {
            if (auto *vol = prop.try_get<Volume>()) {
                volume = vol;
                break;
            }
        }

        if (!volume)
            Throw("ExtremumSpherical requires a sigma_t volume");

        m_volume = volume;
        volume->add_extremum_structure(this);

        m_center = props.get<ScalarPoint3f>("center", ScalarPoint3f(0.f));
        m_rmin = props.get<ScalarFloat>("rmin", 0.f);
        m_rmax = props.get<ScalarFloat>("rmax", 1.f);
        m_resolution = props.get<ScalarVector3i>("resolution", ScalarVector3i(1, 1, 1));
        m_scale = props.get<ScalarFloat>("scale", 1.0f);
        m_to_local = props.get<ScalarAffineTransform4f>("to_world", ScalarAffineTransform4f()).inverse();

        m_sample_method = props.get<ScalarUInt32>("sample_method", 0);

        m_dr = (m_rmax - m_rmin) / m_resolution.x();
        m_idr = dr::rcp(m_dr);

        
        // Set bounding box to enclose the outer sphere
        m_bbox = m_volume->bbox();

        build_grid(m_volume.get(), m_scale);

        Log(Info, "ExtremumSpherical created: resolution=%s, center=%s, "
            "rmin=%f, rmax=%f", m_resolution, m_center, m_rmin, m_rmax);
    }

    void parameters_changed(const std::vector<std::string> &/*keys*/ = {}) override {
        build_grid(m_volume.get(), m_scale);
    }

    Segment sample_segment(
        const Ray3f &ray,
        Float mint, Float maxt,
        Float desired_tau,
        Mask active
    ) const override {

        if constexpr (TraversalType == SphericalTraversalType::RadialOnly) {
            switch (m_sample_method)
            {
            case 0:
                return sample_segment_radial(ray, mint, maxt, desired_tau, active);
            case 1:
                return sample_segment_radial_3(ray, mint, maxt, desired_tau, active);
            default:
                return dr::zeros<Segment>();
            }
            
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
        cb->put("extremum_grid", m_extremum_grid, ParamFlags::NonDifferentiable);
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
            << "  center = " << m_center << "," << std::endl
            << "  rmin = " << m_rmin << "," << std::endl
            << "  rmax = " << m_rmax << std::endl
            << "  fillmin = " << m_fillmin << "," << std::endl
            << "  fillmax = " << m_fillmax << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS(ExtremumSphericalImpl)

private:

    // ------------------------------------------------------------------
    // Grid construction
    // ------------------------------------------------------------------

    void build_grid(const Volume *sigma_t, ScalarFloat scale) {
        // Cell size in normalized [0,1]^3 space
        const ScalarVector3f cell_size = dr::rcp(ScalarVector3f(m_resolution));

        size_t n = dr::prod(m_resolution);

        std::vector<ScalarFloat> extremums;

        // Retrieve data pointer before parallel_for to avoid ref count issues
        auto data = sigma_t->array();

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

                        auto [maj, min] = sigma_t->extremum(data, cell_bounds);

                        dr::scatter(m_extremum_grid,
                                    scale * Vector2f(min, maj),
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

            auto [maj, min] = sigma_t->extremum(data, cell_bounds);

            dr::scatter(m_extremum_grid, min, idx * 2);
            dr::scatter(m_extremum_grid, maj, idx * 2 + 1);
            dr::sync_thread();
        }

        // Retrieve fillmin and fillmax
        Interaction3f it = dr::zeros<Interaction3f>();

        it.p = m_center;
        Float fillmin = sigma_t->eval_1(it, true) * scale;
        
        it.p = m_center + m_rmax + 1;
        Float fillmax = sigma_t->eval_1(it, true) * scale;

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
            m_to_local * ray.d,                  // Normalize direction
            ray.time,
            ray.wavelengths
        );
        // Ray3f local_ray = ray;
        
        Segment result = dr::zeros<Segment>();
        Float tau_acc = 0.f;
        Mask reached = false;
        Float current_t = mint;

        struct LoopState {
            Segment result;
            Mask active, reached;
            Float current_t, tau_acc;

            DRJIT_STRUCT(LoopState, result, active, reached, current_t, tau_acc)
        } ls = { result, active, reached, current_t, tau_acc };

        // Log(Debug, "m_rmin: %f, m_rmax: %f, m_dr: %f, desired_tau: %f", m_rmin, m_rmax, m_dr, desired_tau);
        // Log(Debug, "Start Loop");

        // ray-sphere intersection info
        Vector3f o = local_ray.o - m_center;
        Float o_squared = dr::squared_norm(o);
        Float a = dr::squared_norm(local_ray.d);
        Float b_half = dr::dot(o, local_ray.d);
                

        dr::tie(ls) = dr::while_loop(
            dr::make_tuple(ls),
            [](const LoopState &ls) { return ls.active; },
            [this, &local_ray, maxt, desired_tau, o_squared, a, b_half](LoopState &ls) {
                // Log(Debug, "-------");

                Segment &result   = ls.result;
                Mask    &active   = ls.active;
                Mask    &reached  = ls.reached;
                Float   &current_t = ls.current_t;
                Float   &tau_acc  = ls.tau_acc;

                // Compute radius at current position
                Float eps = math::RayEpsilon<Float>;

                Point3f pos = local_ray(current_t + eps);
                Vector3f oc = pos - m_center;
                Float r = dr::norm(oc);

                Mask fill = r > m_rmax || r < m_rmin;
                Float fill_value = dr::zeros<Float>();
                dr::masked( fill_value, fill && r > m_rmax) = m_fillmax;
                dr::masked( fill_value, fill && r < m_rmin) = m_fillmin;
                
                // Log(Debug, "fill: %f, fill_value: %f", fill, fill_value);
                // Log(Debug, "current_t: %f, pos: %f, r: %f", current_t, pos, r);
                // Log(Debug, "start tau_acc: %f", tau_acc);

                // Determine shell index
                // TODO: precompute dr::rcp(m_dr)
                // TODO: this is probably not working in case we are inside
                //       the fillmin area. 
                //       -> ir   = 0
                //       -> r_lo = rmin
                //       -> r_hi = rmin + dr 
                UInt32 ir = dr::clip(
                    dr::floor2int<UInt32>((r - m_rmin) * m_idr),
                    0u, (uint32_t)(m_resolution.x() - 1)
                );

                // Shell boundary radii
                Float r_lo = m_rmin + Float(ir) * m_dr;
                Float r_hi = r_lo + m_dr;
                
                // Log(Debug, "ir, r_lo, r_hi: %f, %f, %f", ir, r_lo, r_hi);
                // Sphere-ray intersection for both shell boundaries
                // Ray: p(t) = ray.o + t * ray.d
                // Sphere: |p - center|^2 = R^2
                // Vector3f o = ray.o - m_center;
                // Float a = dr::squared_norm(ray.d);
                // Float b_half = dr::dot(o, ray.d);
                
                // Float c_lo = dr::squared_norm(o) - dr::square(r_lo);
                // Float c_hi = dr::squared_norm(o) - dr::square(r_hi);
                Float c_lo = o_squared - dr::square(r_lo);
                Float c_hi = o_squared - dr::square(r_hi);
                
                // Log(Debug, "r_lo: %f, r_hi: %f", r_lo, r_hi);
                // Log(Debug, "a, b_half, c_lo, c_hi: %f, %f, %f, %f", a, b_half, c_lo, c_hi);
                auto [valid_lo, t_lo_near, t_lo_far] =
                    math::solve_quadratic(a, 2.f * b_half, c_lo);
                auto [valid_hi, t_hi_near, t_hi_far] =
                    math::solve_quadratic(a, 2.f * b_half, c_hi);

                // Log(Debug, "lo; valid: %f, near: %f, far: %f", valid_lo, t_lo_near, t_lo_far);
                // Log(Debug, "hi; valid: %f, near: %f, far: %f", valid_hi, t_hi_near, t_hi_far);

                // Find smallest t > current_t + epsilon among the 4 candidates
                Float threshold = current_t + eps;
                Float t_next = maxt;

                // Helper: update t_next if candidate > threshold and < t_next
                auto consider = [&](Float t_cand, Mask valid_cand) DRJIT_INLINE_LAMBDA {
                    Mask use = valid_cand && (t_cand > threshold) && (t_cand < t_next);
                    dr::masked(t_next, use) = t_cand;
                };

                consider(t_lo_near, valid_lo);
                consider(t_lo_far,  valid_lo);
                consider(t_hi_near, valid_hi);
                consider(t_hi_far,  valid_hi);

                // Clamp to maxt
                t_next = dr::minimum(t_next, maxt);

                // Log(Debug, "t_next: %f", t_next);

                // Look up extremum values for this shell
                Vector2f local_extremum = dr::gather<Vector2f>(
                    m_extremum_grid, ir, active && !fill);
                Float local_minorant = dr::select(!fill, local_extremum.x(), fill_value);
                Float local_majorant = dr::select(!fill, local_extremum.y(), fill_value);

                // Log(Debug, "maj: %f min: %f", local_majorant, local_minorant);

                // Accumulate optical depth
                Float segment_length = t_next - current_t;
                Float tau_next = dr::fmadd(local_majorant, segment_length, tau_acc);
                
                // Check if desired tau is reached in this segment
                Mask stops_here = active && (local_majorant > 0.f) &&
                                  (tau_next >= desired_tau) && (t_next <= maxt);

                reached |= stops_here;

                // Log(Debug, "segment_length: %f tau_next: %f, stops_here: %f", segment_length, tau_next, stops_here);

                // Store result for lanes that reached target
                dr::masked(result, stops_here) = ExtremumSegment(
                    current_t, t_next, local_majorant, local_minorant, tau_acc
                );

                // Advance state for lanes that haven't reached target
                dr::masked(current_t, active && !reached) = t_next;
                dr::masked(tau_acc, active && !reached) = tau_next;

                // Continue only if not reached and still in bounds
                active &= !reached && (t_next < maxt);
            },
            "Spherical Shell Traversal"
        );

        // For lanes that didn't reach, invalidate the segment
        dr::masked(ls.result.tmin, !ls.reached) = dr::Infinity<Float>;

        return ls.result;
    }

    Segment sample_segment_radial_3(
        const Ray3f &ray,
        Float mint, Float maxt,
        Float desired_tau,
        Mask active
    ) const {

        // const ScalarVector3f extents = m_bbox.extents();
        Ray3f local_ray(
            m_to_local * ray.o,  // Normalize origin
            m_to_local * ray.d,                  // Normalize direction
            ray.time,
            ray.wavelengths
        );
        // Ray3f local_ray = ray;

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

        // Find the current/next intersection (use this to calculate the midpoint too)
        Point3f pos = local_ray(mint+dr::Epsilon<Float>*10.f);
        Vector3f oc = pos - m_center;
        Float r = dr::norm(oc);

        // Calculate the initial layer index from which we will step through layers.
        Int32 layer_idx = dr::clip(dr::floor2int<Int32>((r - m_rmin) * m_idr), -1, m_resolution.x());
        Mask passed_midpoint = dr::dot((m_center - pos), local_ray.d) < 0;
        UInt32 shell_padding = dr::select(passed_midpoint, 1, 0);
        Int32 step = dr::select(passed_midpoint, 1, -1);

        // Log(Debug, "r: %f, (r-m_rmin)*m_idr: %f, layer_idx: %f", r, dr::floor2int<UInt32>((r - m_rmin) * m_idr), layer_idx);
        // Calculate distance to midpoint, norm of r cancels out.
        // Float t_midpoint = dr::dot(dr::normalize(ray.d), oc);
        // Float r_midpoint = dr::sqrt(t_midpoint * t_midpoint + r * r) ;
        // Int32 midpoint_idx = dr::clip(dr::floor2int<UInt32>((r_midpoint - m_rmin) * m_idr), -1, m_resolution.x());

        struct LoopState {
            Segment result;
            Mask active, reached;
            Float current_t, tau_acc;
            Int32 layer_idx;
            Int32 step;
            UInt32 padding;

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
            [this, maxt, desired_tau, o_squared, a, b_half](LoopState &ls) {
            // Log(Debug, "---------");
            Segment &result    = ls.result;
            Mask    &active    = ls.active;
            Mask    &reached   = ls.reached;
            Float   &current_t = ls.current_t;
            Float   &tau_acc   = ls.tau_acc;
            Int32   &layer_idx = ls.layer_idx;
            Int32   &step      = ls.step;
            UInt32  &padding   = ls.padding;

            // Compute radius at current position
            Float eps = math::RayEpsilon<Float>;
            
            // Passed midpoint == exiting the concentric spheres
            Int32 shell_idx = dr::clip(layer_idx + padding, 0, m_resolution.x());
            // Log(Debug, "current_t: %f, layer_idx: %f, shell_idx: %f", current_t, layer_idx, shell_idx);
            // Log(Debug, "passed_mp: %f", passed_midpoint);

            // Boundary condition of rmin and rmax
            Float fill_value = -1.f;
            dr::masked(fill_value, layer_idx < 0) = m_fillmin;
            dr::masked(fill_value, layer_idx >= m_resolution.x()) = m_fillmax;
            Mask fill = fill_value >= 0.f;

            // Test intersection with the shell
            Float r_test = m_rmin + Float(shell_idx) * m_dr;
            Float c = o_squared - dr::square(r_test);
            auto [valid_test, t_test_near, t_test_far] = math::solve_quadratic(a, 2.f * b_half, c);
            // Log(Debug, "r_test: %f, valid_test: %f", r_test, valid_test);
            // Log(Debug, "a, b_half, c_lo, c_hi: %f, %f, %f, %f", a, b_half, c);

            // Find smallest t > current_t + epsilon among the 4 candidates
            Float threshold = current_t + eps;
            Float t_next = maxt;

            // Special case at midpoint where we miss the intersection with the 
            // shell, bump padding to test outer shell, set step to increase
            // layer index and continue to next iteration.
            dr::masked(padding, !valid_test) = 1;
            dr::masked(step, !valid_test) = 1;

            Mask update = active && valid_test; 
            // Log(Debug, "update: %f", update);

            if( dr::any_or<true>(update) ) {
                
                // Helper: update t_next if candidate > threshold and < t_next
                auto consider = [&](Float t_cand, Mask valid_cand) DRJIT_INLINE_LAMBDA {
                    Mask use = valid_cand && (t_cand > threshold) && (t_cand < t_next);
                    dr::masked(t_next, use) = t_cand;
                };

                consider(t_test_near, valid_test);
                consider(t_test_far,  valid_test);

                // Clamp to maxt
                t_next = dr::minimum(t_next, maxt);
                // Log(Debug, "t_next: %f", t_next);
                
                // Look up extremum values for this shell
                Vector2f local_extremum = dr::gather<Vector2f>(
                    m_extremum_grid, layer_idx, active && !fill);
                Float local_minorant = dr::select(!fill, local_extremum.x(), fill_value);
                Float local_majorant = dr::select(!fill, local_extremum.y(), fill_value);

                // Accumulate optical depth
                Float segment_length = t_next - current_t;
                Float tau_next = dr::fmadd(local_majorant, segment_length, tau_acc);
                
                // Check if desired tau is reached in this segment
                Mask stops_here = active && (local_majorant > 0.f) &&
                                    (tau_next >= desired_tau) && (t_next <= maxt);

                reached |= stops_here;
                // Log(Debug, "stops_here: %f, reached: %f", stops_here, reached);
                // Store result for lanes that reached target
                dr::masked(result, stops_here) = ExtremumSegment(
                    current_t, t_next, local_majorant, local_minorant, tau_acc
                );

                // Advance state for lanes that haven't reached target
                dr::masked(current_t, update && !reached) = t_next;
                dr::masked(tau_acc, update && !reached) = tau_next;
                dr::masked(layer_idx, update && !reached) += step;

                // Log(Debug, "current_t: %f, layer_idx: %f", current_t, layer_idx);
                // Log(Debug, "count: %f", ls.count);

                // Continue only if not reached and still in bounds
                active &= (!reached && (t_next < maxt)) ;
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
