#include <mitsuba/core/properties.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/eradiate/extremum.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/volumegrid.h>
#include <nanothread/nanothread.h>

NAMESPACE_BEGIN(mitsuba)

// @TODO:
// - Currently only works when perfectly aligned to volume, see if we can extend it :
//   For now, we will assume the extremum grid to be perfectly aligned to the volume 
//   as it simplifies many issues related to building grid outside of the volume 
//   extent (wraping mode).
//   Wraping mode is not available to the user and therefore not to the majorant.
// - Currently works in 1 channel, see how we can extend to multiple channels.
// - Grid building limited to one volume, consider extending to multiple.
//   This is contigent to the axis-aligned volume problem.
// - For build_grid -> detect when volume's res == supergrid res to do a memcpy
//                  -> detect when supergrid res == 1x1x1 to get max value.
// - Is the m_bbox really needed for extremum grids?
// - Add min to the volume functions.

/**!
.. _extremum-extremum_grid:

Extremum grid structure (:monosp:`extremum_grid`)
--------------------------------------------------

.. pluginparameters::

 * - volume
   - |volume|
   - Extinction coefficient volume to build extremum grid from
   - |exposed|

 * - to_world
   - |transform|
   - Specifies a 4x4 transformation matrix of the underlying volume.

 * - scale
   - |float|
   - Scale factor for the extremum values (Default: 1.0).

 * - resolution
   - |vector|
   - Grid resolution along the XYZ axis. Does not have to be a multiple of 
     the underlying volume. `adaptive` is mutually exclusive with `resolution`. 
     (Default: [1,1,1]).

 * - adaptive
   - |bool|
   - Flags the use of an adaptive method to find a suitable majorant grid resolution.
     Note that this search is costly and needs to run on rebuilds. `adaptive` is 
     mutually exclusive with `resolution`. (Default: false)

This plugin creates a regular grid structure storing local extremum values for 
efficient delta tracking in heterogeneous media. The grid is constructed
by querying the extinction volume's extrema over each grid cell.

At runtime, DDA (Digital Differential Analyzer) traversal through the grid provides
tight-fitting local extrema, dramatically reducing null collisions in media with
high spatial variance (e.g., clouds, fog).
*/

template <typename Float, typename Spectrum>
class ExtremumGrid final : public ExtremumStructure<Float, Spectrum> {
public:
    MI_IMPORT_BASE(ExtremumStructure, m_bbox)
    MI_IMPORT_TYPES(Volume)

    using Segment = ExtremumSegment;
    using FloatStorage = DynamicBuffer<Float>;

    ExtremumGrid(const Properties &props) : Base(props) {
        // Volume Parameters
        m_volume = nullptr;
        for (auto &prop : props.objects()) {
            if (auto *vol = prop.try_get<Volume>()) {
                m_volume = vol;
                break;
            }
        }

        if (!m_volume)
            Throw("ExtremumGrid requires at least one volume");
        
        // Register the extremum structure to the volume to trigger 
        // parameter_changed when the volume is modified.
        m_volume->add_extremum_structure(this);
        m_bbox = m_volume->bbox(); 

        m_to_local = props.get<ScalarAffineTransform4f>("to_world", ScalarAffineTransform4f()).inverse();
        m_scale = props.get<ScalarFloat>("scale", 1.0f);

        // Resolution Parameters
        if (props.has_property("adpative") &&
            props.has_property("resolution")) {
            Throw("`adaptive_resolution` and `resolution` are mutually "
                  "exclusive.");
        } else if (props.has_property("adaptive_resolution")) {
            if (props.get<bool>("adaptive_resolution")) {
                m_adaptive   = true;
                m_resolution = find_resolution(m_volume);
            }
        } else if (props.has_property("resolution")) {
            m_adaptive   = false;
            m_resolution = props.get<ScalarVector3i>("resolution");
        }
        m_cell_size = 1.f / ScalarVector3f(m_resolution);

        build_grid(m_volume.get(), m_resolution);

        Log(Info, "ExtremumGrid created: resolution=%s, bbox=%s",
            m_resolution, m_bbox);
    }

    void parameters_changed(const std::vector<std::string> &/*keys*/ = {}) override {
        if (m_adaptive) {
            find_resolution(m_volume.get());
        }
        build_grid(m_volume.get(), m_resolution);
    }

    std::tuple<Segment, Float> sample_segment(
        const Ray3f &ray,
        Float mint, 
        Float maxt,
        Float target_od,
        Mask active
    ) const override {
        // MI_MASKED_FUNCTION(ProfilerPhase::MediumSample, active);

        Segment result = dr::zeros<Segment>();

        // Currently assuming that the majorant aligns perfectly with the 
        // volume and that values outside the bbox cannot be evaluated.
        // Transform ray to local grid coordinates [0,res]³
        Vector3f res = Vector3f(m_resolution);
        Ray3f local_ray(
            (m_to_local * ray.o)*res,  // Normalize origin
            (m_to_local * ray.d)*res,  // Normalize direction
            ray.time,
            ray.wavelengths
        );
        Vector3f rcp_d = dr::rcp(local_ray.d);
        auto inf_t = local_ray.d == 0.f;
        auto d_pos = local_ray.d >= 0.f;

        // Per-axis intersection of the ray with the grid bounds
        Vector3f t_min_v = -local_ray.o * rcp_d;
        Vector3f t_max_v = (res - local_ray.o) * rcp_d;
        Vector3f t_min_v2 = dr::minimum(t_min_v, t_max_v);
        Vector3f t_max_v2 = dr::maximum(t_min_v, t_max_v);

        // Disable extent computation for dims where the ray direction is zero
        dr::masked(t_max_v2, inf_t) = dr::Infinity<Float>;
        dr::masked(t_min_v2, inf_t) = -dr::Infinity<Float>;

        // Reduce constraints to a single ray interval
        Float t_min = dr::maximum(dr::max(t_min_v2), mint);
        Float t_max = dr::minimum(dr::min(t_max_v2), maxt);
        mint=t_min;
        maxt=t_max;

        // Only run the DDA algorithm if the interval is nonempty
        active &= (t_max > t_min) & dr::isfinite(t_max);

        // Deactivate rays that have zero direction along any axis
        // and whose origin along that axis is outside the grid bounds
        active &= dr::all(!inf_t || ((0.f <= local_ray.o) && (local_ray.o <= res)));

        // Advance the ray to the start of the interval
        local_ray.o = dr::fmadd(local_ray.d, t_min, local_ray.o);
        t_max = t_max - t_min;
        t_min = 0.f; // type: ignore

        // Compute the integer step direction
        Vector3i step = dr::select(d_pos, 1, -1);
        Vector3f offset = dr::select(d_pos, 0.f, 1.f);
        Vector3f abs_rcp_d = abs(rcp_d);

        // Integer grid coordinates
        // Vector3i pi = dr::clip(local_ray.o * m_resolution, 0, m_resolution - 1);
        Vector3i pi = dr::clip(Vector3i(local_ray.o), 0, m_resolution - 1);

        // Fractional entry position
        Vector3f p0 = local_ray.o - Vector3f(pi);
        // Log(Debug, "m_resolution: %f", m_resolution);
        // Log(Debug, "l_ray.o: %f, l_ray.d: %f", local_ray.o, local_ray.d);
        // Log(Debug, "pi: %f, p0: %f, cell_size: %f", pi, p0);
        // Step size to next interaction
        Vector3f dt_v = dr::select(d_pos, dr::fmadd(-p0, rcp_d, rcp_d), -p0 * rcp_d);
        dr::masked(dt_v, inf_t) = dr::Infinity<Float>;

        struct LoopState {
            Segment result;
            Mask active;
            Vector3f dt_v;
            Vector3f p0;
            Vector3i pi;
            Float t_rem;
            Float tau_acc;

            DRJIT_STRUCT(LoopState, result, active, dt_v, p0, pi, t_rem, tau_acc)
        } ls = {
            result,
            active,
            dt_v,
            p0,
            pi,
            t_max,
            /* tau_acc = */ 0.f
        };
        
        dr::tie(ls) = dr::while_loop(
            dr::make_tuple(ls),
            [](const LoopState& ls) { return ls.active; },
            [this, local_ray, step, abs_rcp_d, t_max, target_od, mint, offset](LoopState& ls) {
            // Log(Debug, "-----");
            
            Segment& result = ls.result;
            Mask& active    = ls.active;
            Vector3f& dt_v  = ls.dt_v;
            Vector3f& p0    = ls.p0;
            Vector3i& pi    = ls.pi;
            Float& t_rem = ls.t_rem; 
            Float& tau_acc = ls.tau_acc; 
            
            // Select the smallest step. It's possible that dt == 0 when starting
            // directly on a grid line.
            Float dt = dr::minimum(dr::min(dt_v), t_rem); // what type is dt and t_rem?
            auto mask = dt_v == dt;
            // Log(Debug, "t_rem: %f, dt_v: %f", t_rem, dt_v);

            // Compute an updated position
            Vector3f p1 = dr::fmadd(local_ray.d, dt, p0);

            // Note: not multiplying the index by 2 because we gather using Vector2f. 
            UInt32 idx = dr::fmadd( 
                        dr::fmadd( pi.z(), m_resolution.y(), pi.y() ), 
                        m_resolution.x(), pi.x()
                    );
            // Log(Debug, "pi: %f, idx: %f", pi, idx);

            Vector2f extremum = dr::gather<Vector2f>(m_extremum_grid, idx);
            const Float minorant = extremum.x();
            const Float majorant = extremum.y();                    
            // Log(Debug, "local_minorant: %f, local_majorant: %f", minorant, majorant);

            // Accumulate optical depth
            Float tau_next = dr::fmadd(majorant, dt, tau_acc);

            // Check if desired tau reached in this segment
            Mask exit = active && (tau_next >= target_od);

            if (dr::any_or<true>(exit)){
                // Store result for lanes that reached target
                Float t_curr = mint + t_max - t_rem;
                dr::masked(result, exit) = ExtremumSegment(
                    t_curr, 
                    t_curr + dt, 
                    majorant, 
                    minorant
                );
            }
            // Log(Debug, "mint: %f, t_max: %f, t_rem: %f, dt: %f", mint, t_max, t_rem, dt);
            // Log(Debug, "tau_next: %f, t_curr: %f, t_next: %f", tau_next, t_curr, t_curr + dt);
            // Advance
            dt_v = dr::select(mask, abs_rcp_d, dt_v - dt);
            // dr::masked(p1, mask) = dr::select(local_ray.d >= 0, 0.f, 1.f);
            dr::masked(p1, mask) = offset;
            dr::masked(pi, mask) += step;
            t_rem -= dt;
            dr::masked(tau_acc, !exit) = tau_next;
            
            active &= dr::all((pi >= 0) && (pi < m_resolution)) && (t_rem > 0.f) && !exit;
            // Log(Debug, "t_rem: %f, pi: %f, exit: %f, active: %f", t_rem, pi, exit, active );
        },
        "DDA Travesal");


        // For lanes that didn't reach, set tmin to infinity to invalidate the segment
        // dr::masked(ls.result.tmin, !ls.reached) = dr::Infinity<Float>;

        return {ls.result, ls.tau_acc};
    }

    std::tuple<Float, Float> eval_1(
        const Interaction3f &it, 
        Mask active) const override {
        Point3f po = m_to_local*it.p;
        Vector3i pi = dr::clip(Vector3i(po * m_resolution), 0, m_resolution - 1);
        UInt32 idx = dr::fmadd( 
                        dr::fmadd( pi.z(), m_resolution.y(), pi.y() ), 
                        m_resolution.x(), pi.x()
                    );
        Vector2f extremum = dr::gather<Vector2f>(m_extremum_grid, idx, active);
        return {extremum.x(), extremum.y()};
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

    MI_DECLARE_CLASS(ExtremumStructure)

private:

    /**
     * \brief Find the best resolution for a given volume.
     *
     * This method performs a ternary search to find the best extremum grid 
     * resoluton for a given volume. Note that this is currently a costly
     * operation.
     */
    ScalarVector3i find_resolution(const Volume *volume){
        using ScalarIndex =  DynamicBuffer<ScalarUInt32>;
        using ScalarFloatStorage =  DynamicBuffer<ScalarFloat>;
        ScalarVector3i result = dr::full<ScalarVector3i>(1);
        ScalarVector3i fine_res = volume->resolution();
        
        ScalarAffineTransform4f to_world = m_to_local.inverse();
        ScalarVector3f dim;

        // Extract scale from to_world transform. Note this does not work with shear
        for (size_t i = 0; i < 3; ++i) {
            ScalarVector3f basis;
            for (size_t j = 0; j < 3; ++j) {
                basis[j] = to_world.matrix.entry(i, j);
            }
            dim[i] = dr::norm(basis);
        }
 
        Float denom = dim.x()*dim.y() + dim.x()*dim.z() + dim.y()*dim.z();  

        // Cost function as defined in Yue et al. 2011
        auto cost_fn = [denom, dim](FloatStorage& grid, ScalarVector3i resolution){
            // Select only the majorant from the extremum grid
            ScalarIndex idx = dr::arange<ScalarIndex>(dr::prod(resolution))*2+1;
            ScalarFloatStorage majorant = dr::gather<ScalarFloatStorage>(grid, idx);
            
            Float res_prod = dr::prod(resolution);
            Float dim_prod = dr::prod(dim);
            Float grid_sum = dr::sum(majorant);

            Float cost_iter = 2.f * (dim_prod / res_prod) * grid_sum;

            Float cost_part = (resolution.x() - 1)*dim.y()*dim.z()
                            + (resolution.y() - 1)*dim.x()*dim.z()
                            + (resolution.z() - 1)*dim.x()*dim.y(); 

            // Timing parameters: 
            // t_iter: cost incurred by having to perform a null collision loop
            // t_rewind: cost incurred by having to perform a dda loop
            // Dependent on optimization and target build, left for further optimization.
            Float t_iter = 1.f;
            Float t_rewind = 1.f;
            return (t_iter * cost_iter + t_rewind * cost_part)/denom;
        };

        auto to_scalar = [](Float f){
            if constexpr (dr::is_jit_v<Float>)
                return f[0];
            else
                return f;
        };

        // Ternary search over the xyz dimensions to find the best resolution.
        result = dr::maximum( fine_res/2, 1);
        for (uint8_t dim = 0; dim < 3; ++dim){
            ScalarUInt32 left = 1;
            ScalarUInt32 right = fine_res[dim];
            ScalarFloat left_cost, right_cost;

            while( (right - left) > 3) {
                // Choose new boundaries
                ScalarUInt32 m1 = left + (right-left)/3;
                ScalarUInt32 m2 = right - (right-left)/3; 

                ScalarVector3i left_res = result, right_res = result;
                left_res[dim]  = m1;
                right_res[dim] = m2;
    
                // Evaluate cost at the chosen boundaries
                build_grid(volume, left_res);
                left_cost = to_scalar(cost_fn(m_extremum_grid, left_res));

                build_grid(volume, right_res);
                right_cost = to_scalar(cost_fn(m_extremum_grid, right_res));

                if(left_cost < right_cost)
                    right = m2 - 1;
                else 
                    left = m1 + 1;
            }

            result[dim] = left_cost < right_cost ? left : right;
        }
        Log(Warn, "optimal resolution: %f", result);
        return result;
    }

    /**
     * \brief Build the extremum grid from a volume
     *
     * This method constructs a lower-resolution grid where each cell stores
     * the majorant (maximum) extinction value over the corresponding region
     * of the high-resolution volume.
     */
    void build_grid(const Volume *volume, ScalarVector3i resolution) {

        // local space supergrid cell size
        const ScalarVector3f cell_size = dr::rcp(ScalarVector3f(resolution));

        ScalarVector2f safety_factor(
            1.f - dr::Epsilon<Float>, 
            1.f + dr::Epsilon<Float>
        );
        
        // Log(Debug, "Building grid =======");
        // Log(Debug, "m_bbox: %f, cell_size: %f", m_bbox, cell_size);

        // Allocate extremum grid data
        size_t n = dr::prod(resolution);
        
        std::vector<ScalarFloat> extremums;

        // IMPORTANT: need to retrieve the data before to avoid ref count slow-down
        auto data = volume->array();

        size_t n_threads = pool_size() + 1;
        size_t grain_size = std::max( n / (4 * n_threads), (size_t) 1 );

        m_extremum_grid = dr::empty<FloatStorage>(n*2);

        // Early return if using the global majorant.
        if (n == 1) {
            ScalarFloat max = volume->max();
            dr::scatter(
                m_extremum_grid, 
                m_scale * Vector2f(0.f, max) * safety_factor, 
                UInt32(0));
            return;
        }

        if constexpr (!dr::is_jit_v<Float>) {
            
            dr::parallel_for(
                dr::blocked_range<size_t>(0, n, grain_size),
                [&](const dr::blocked_range<size_t> &range) {
                    // Recover x, y, z from block start (one-time div/mod per block)

                    // Log(Info, "begin: %f, end: %f", range.begin(), range.end());
                    for (auto idx = range.begin(); idx != range.end(); ++idx) {
                        // Store in linear array (Z-slowest, X-fastest)
                        int32_t x = idx % resolution.x() ; 
                        int32_t y = (idx / resolution.x())  % resolution.y(); 
                        int32_t z =  idx / (resolution.x() * resolution.y()); 

                        
                        ScalarPoint3f cell_min = ScalarVector3f(x, y, z) * cell_size;
                        ScalarPoint3f cell_max = cell_min + cell_size;
                        ScalarBoundingBox3f cell_bounds(
                            cell_min + math::RayEpsilon<Float>, 
                            cell_max - math::RayEpsilon<Float>
                        );

                        // Query volume for local extremum, currently assume local bounds.
                        auto [maj, min] = volume->extremum(data, cell_bounds);

                        dr::scatter(
                            m_extremum_grid, 
                            m_scale * Vector2f(min, maj) * safety_factor, 
                            UInt32(idx));
                    }
                }
            );
        } else {
            
            UInt32 idx = dr::arange<UInt32>((uint32_t) n);

            UInt32 x = idx % resolution.x() ; 
            UInt32 y = (idx / resolution.x())  % resolution.y(); 
            UInt32 z =  idx / (resolution.x() * resolution.y()); 

            Point3f cell_min = Vector3f(x, y, z) * cell_size;
            Point3f cell_max = cell_min + cell_size;
            BoundingBox3f cell_bounds(
                            cell_min + math::RayEpsilon<Float>, 
                            cell_max - math::RayEpsilon<Float>
                        );

            auto [maj, min] = volume->extremum(data, cell_bounds);

            dr::scatter(m_extremum_grid, m_scale * min * safety_factor.x(), idx*2);
            dr::scatter(m_extremum_grid, m_scale * maj * safety_factor.y(), idx*2+1);
            dr::sync_thread();
        }

        
        Log(Info, "Extremum grid constructed successfully");
    }

private:
    /// Grid storing pre-computed local majorants
    FloatStorage m_extremum_grid;
    
    ref<Volume> m_volume;
    ScalarFloat m_scale;

    bool m_adaptive;
    ScalarVector3i m_resolution;
    ScalarVector3f m_cell_size;
    ScalarAffineTransform4f m_to_local;
};

MI_EXPORT_PLUGIN(ExtremumGrid)
NAMESPACE_END(mitsuba)
