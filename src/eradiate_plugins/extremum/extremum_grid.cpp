#include <mitsuba/core/properties.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/eradiate/extremum.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/volumegrid.h>
#include <nanothread/nanothread.h>

NAMESPACE_BEGIN(mitsuba)


// @TODO:
// - Currently only works on axis aligned volumes, see if we can extend it :
//   For now, we will assume the volume to be axis aligned as it simplifies many
//   issues related to building grid outside of the volume extent (wraping mode).
//   Wraping mode is not available to the user and therefore not to the majorant.
// - currently works in 1D, see how we can extend to multiple channels.
// - Grid building limited to one volume, consider extending to multiple.
//   This is contigent to the axis-aligned volume problem.
// - allow for volumes with accel=True
// - work on adaptive resolution
// - reintroduce the safety factor
// - for build_grid -> detect that the volume's res is the same as the supergrid 
//   and do a memcpy.

/**!
.. _extremum-extremum_grid:

Extremum grid structure (:monosp:`extremum_grid`)
--------------------------------------------------

.. pluginparameters::

 * - sigma_t
   - |volume|
   - Extinction coefficient volume to build extremum grid from
   - |exposed|

 * - scale
   - |float|
   - Scale factor for extinction coefficients (Default: 1.0)

 * - resolution_factor
   - |int|
   - Grid resolution divisor. The extremum grid resolution will be
     sigma_t_resolution / resolution_factor. If 0, automatically determined
     via heuristic. (Default: 0)

 * - safety_factor
   - |float|
   - Safety margin for majorants. All majorants are multiplied by this factor
     to ensure they remain conservative despite floating-point rounding.
     (Default: 1.01)

This plugin creates a regular grid structure storing local majorant (and minorant)
values for efficient delta tracking in heterogeneous media. The grid is constructed
by querying the extinction volume's extrema over each grid cell.

At runtime, DDA (Digital Differential Analyzer) traversal through the grid provides
tight-fitting local majorants, dramatically reducing null collisions in media with
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
        // Get volume
        // TODO: Currently limited to one volume, will need to be extended to multiple.
        ref<Volume> volume = nullptr;
        for (auto &prop : props.objects()) {
            if (auto *vol = prop.try_get<Volume>()) {
                volume = vol;
                break;
            }
        }

        if (!volume)
            Throw("ExtremumGrid requires at least one volume");
        
        m_bbox = volume->bbox(); 
        m_to_local = (ScalarAffineTransform4f::scale(m_bbox.extents()) * ScalarAffineTransform4f::translate(m_bbox.min)).inverse();

        // Get parameters
        ScalarFloat scale = props.get<ScalarFloat>("scale", 1.0f);

        if (props.has_property("resolution_factor") && props.has_property("resolution")) {
            
            Throw("`resolution_factor` and `resolution` are mutually exclusive.");

        } else if (props.has_property("resolution_factor")) {

            const ScalarVector3i full_res = volume->resolution();

            size_t resolution_factor = props.get<size_t>("resolution_factor", 0);

            // Determine grid resolution
            if (resolution_factor == 0) {
                // Heuristic: use cube root of total voxels / 1M, clamped to [1, 16]
                size_t total_voxels = dr::prod(full_res);
                size_t heuristic_factor = std::max(
                    1,
                    static_cast<int>(std::pow(total_voxels / 1000000.0f, 1.0f/3.0f))
                );
                resolution_factor = dr::clip(heuristic_factor, 1, 16);

                Log(Info, "Auto-determined extremum grid resolution factor: %d", resolution_factor);
            }

            m_resolution = full_res / ScalarVector3i(resolution_factor);
            m_resolution = dr::maximum(m_resolution, 1);

            Log(Info, "Building extremum grid %s from sigma_t grid %s (factor=%d)",
            m_resolution, full_res, resolution_factor);

        } else if (props.has_property("resolution")) {
            m_resolution = props.get<ScalarVector3i>("resolution");
            Log(Debug, "resolution: %f", m_resolution);
        }

        build_grid(volume.get(), scale);

        // TODO: remove useless logs
        Log(Info, "ExtremumGrid created: resolution=%s, bbox=%s",
            m_resolution, m_bbox);
    }

    Segment sample_segment(
        const Ray3f &ray,
        Float mint, Float maxt,
        Float desired_tau,
        Mask active
    ) const override {
        // MI_MASKED_FUNCTION(ProfilerPhase::MediumSample, active);

        Segment result = dr::zeros<Segment>();

        // Intersect AABB
        auto [aabb_its, local_mint, local_maxt] = m_bbox.ray_intersect(ray);
        aabb_its &= (dr::isfinite(local_mint) || dr::isfinite(local_maxt));
        active &= aabb_its;
        dr::masked(local_mint, !active) = 0.f;
        dr::masked(local_maxt, !active) = dr::Infinity<Float>;

        mint = dr::maximum(local_mint, mint);
        maxt = dr::minimum(local_maxt, maxt);

        // Prepare DDA traversal
        Float dda_t; 
        Vector3f dda_tmax, dda_tdelta;
        Mask valid;
        std::tie(dda_t, dda_tmax, dda_tdelta, valid) =
            prepare_dda_traversal(ray, mint, maxt, active);
            
        active &= valid;

        Log(Debug, "dda_t: %f, dda_tmax: %f, ddat_tdelta: %f, valid: %f", dda_t, dda_tmax, dda_tdelta, valid);

        // Traverse grid with DDA until desired_tau is reached
        Float tau_acc = 0.f;
        Mask reached = false;

        Log(Debug, "Start Loop");

        struct LoopState {
            Segment result;
            Mask active;
            Mask reached;
            Ray3f ray;
            Float tau_acc;
            Float dda_t; 
            Vector3f dda_tmax;

            DRJIT_STRUCT(LoopState, result, active, reached, ray, tau_acc, dda_t, dda_tmax)
        } ls = {
            result,
            active,
            reached,
            ray,
            tau_acc,
            dda_t,
            dda_tmax
        };
        
        dr::tie(ls) = dr::while_loop(
            dr::make_tuple(ls),
            [](const LoopState& ls) { return ls.active; },
            [this, dda_tdelta, maxt, desired_tau](LoopState& ls) {
            
            Segment& result = ls.result;
            Mask& active = ls.active;
            Mask& reached = ls.reached;
            Ray3f& ray = ls.ray;
            Float& tau_acc = ls.tau_acc;
            Float& dda_t = ls.dda_t;
            Vector3f& dda_tmax = ls.dda_tmax; 
            Log(Debug, "============================");
            Log(Debug, "reached: %f, tau_acc: %f, dda_t: %f", reached, tau_acc, dda_t);
            // Find next voxel boundary (minimum of dda_tmax components)
            Float t_next = dr::minimum(dr::minimum(dda_tmax.x(), dda_tmax.y()), dda_tmax.z());
            t_next = dr::minimum(t_next, maxt);

            Log(Debug, "t_next: %f", t_next);

            // Determine which axis advances
            Vector3f tmax_update = dr::zeros<Vector3f>();
            Mask assigned = false;
            for (size_t k = 0; k < 3; ++k) {
                Mask hit_axis = (dda_tmax[k] == t_next) && !assigned;
                tmax_update[k] = dr::select(hit_axis, dda_tdelta[k], 0.f);
                assigned |= hit_axis;
            }

            
            // Lookup local majorant at cell center
            Float t_mid = 0.5f * (dda_t + t_next);
            Log(Debug, "assigned: %f, tmax_update: %f", assigned, tmax_update);
            Log(Debug, "t_mid: %f, ray(t_mid): %f", t_mid, ray(t_mid));

            // Retrieve index from cell center position.
            const Vector3f pos_i = dr::floor2int<Vector3i>((m_to_local * ray(t_mid)) * Vector3f(m_resolution));
            
            // Note: not multiplying the index by 2 because we gather using Vector2f. 
            UInt32 index( dr::fmadd( 
                        dr::fmadd( UInt32(pos_i.z()), m_resolution.y(), UInt32(pos_i.y()) ), 
                        m_resolution.x(), 
                        UInt32(pos_i.x())
                    ));
            Log(Debug, "pos_i: %f, index: %f", pos_i, index);

            Vector2f local_extremum = dr::gather<Vector2f>(m_extremum_grid, index);
            const Float local_minorant = local_extremum.x();
            const Float local_majorant = local_extremum.y();                    
            Log(Debug, "local_minorant: %f, local_majorant: %f", local_minorant, local_majorant);

            // Accumulate optical depth
            Float segment_length = t_next - dda_t;
            Float tau_next = tau_acc + local_majorant * segment_length;

            // Check if desired tau reached in this segment
            Mask stops_here = active && (local_majorant > 0.f) &&
                             (tau_next >= desired_tau) && (t_next <= maxt);

            Log(Debug, "tau_next: %f", tau_next);
            // Update reached status
            reached |= stops_here;

            // Store result for lanes that reached target
            dr::masked(result, stops_here) = ExtremumSegment(
                dda_t, t_next, local_majorant, local_minorant, tau_acc
            );

            // Advance DDA state
            dr::masked(dda_t, active && !reached) = t_next;
            dr::masked(dda_tmax, active && !reached) = dda_tmax + tmax_update;
            dr::masked(tau_acc, active && !reached) = tau_next;

            // Continue only if not reached and still in bounds
            active &= !reached && (t_next < maxt);
        },
        "DDA Travesal");

        // Here we choose to not compute the precise interaction point  
        // in this function and defer it to the medium. This is for to
        // allow for more flexibility on the sampling method. 
        // Float t_precise = ls.dda_t +
        //         (desired_tau - ls.tau_acc) / dr::maximum(local_majorant, dr::Epsilon<Float>);

        // For lanes that didn't reach, set tmin to infinity to invalidate the segment
        dr::masked(ls.result.tmin, !ls.reached) = dr::Infinity<Float>;

        return ls.result;
    }

    void traverse(TraversalCallback *cb) override {
        cb->put("extremum_grid", m_extremum_grid, ParamFlags::NonDifferentiable);
        cb->put("resolution", m_resolution, ParamFlags::NonDifferentiable);
        Base::traverse(cb);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "ExtremumGrid[" << std::endl
            << "  resolution = " << m_resolution << "," << std::endl
            << "  bbox = " << m_bbox << "," << std::endl
            << "  safety_factor = " << m_safety_factor << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS(ExtremumStructure)

private:

    /**
     * \brief Build the extremum grid from a sigma_t volume
     *
     * This method constructs a lower-resolution grid where each cell stores
     * the majorant (maximum) extinction value over the corresponding region
     * of the high-resolution sigma_t volume.
     */
    // void build_grid(const Volume *sigma_t, ScalarFloat scale) {
    //     using FloatX = DynamicBuffer<ScalarFloat>;

    //     m_bbox = sigma_t->bbox();
    //     // local space supergrid cell size
    //     const ScalarVector3f cell_size = dr::rcp(ScalarVector3f(m_resolution));
        
    //     // Log(Debug, "Building grid =======");
    //     // Log(Debug, "m_bbox: %f, cell_size: %f", m_bbox, cell_size);

    //     // Allocate extremum grid data
    //     size_t n = dr::prod(m_resolution);
    //     std::vector<ScalarFloat> extremums;
    //     // std::unique_ptr<ScalarFloat[]> extremums(new ScalarFloat[n*2]);

    //     // avoid unecessary allocation in scalar mode
    //     if constexpr (!dr::is_jit_v<Float>) {
    //         m_extremum_grid = dr::empty<FloatStorage>(n*2);
    //     } else {
    //         extremums.resize(n*2);
    //     }

    //     // Build grid (consider parallelizing)
    //     for (int32_t z = 0; z < m_resolution.z(); ++z) {
    //         for (int32_t y = 0; y < m_resolution.y(); ++y) {
    //             for (int32_t x = 0; x < m_resolution.x(); ++x) {

    //     // for (int32_t x = 0; x < m_resolution.x(); ++x) {
    //     //     for (int32_t y = 0; y < m_resolution.y(); ++y) {
    //     //         for (int32_t z = 0; z < m_resolution.z(); ++z) {
    //                 // Compute cell bounding box in local space
    //                 ScalarPoint3f cell_min = ScalarVector3f(x, y, z) * cell_size;
    //                 ScalarPoint3f cell_max = cell_min + cell_size;
    //                 ScalarBoundingBox3f cell_bounds(
    //                     cell_min + math::RayEpsilon<Float>, 
    //                     cell_max - math::RayEpsilon<Float>
    //                 );
                    
    //                 // Log(Debug, "[SuperGrid] x,y,z: %f, %f, %f -------", x, y, z);
    //                 // Log(Debug, "[SuperGrid] cell_bounds: %f, %f",cell_bounds.min, cell_bounds.max);

    //                 // Query volume for local extremum, currently assume local bounds.
    //                 auto [maj, min] = sigma_t->extremum(nullptr, cell_bounds);

    //                 // Store in linear array (X-slowest, Z-fastest) <-- current
    //                 // size_t idx = z
    //                 //              + y * m_resolution.z() 
    //                 //              + x * m_resolution.z() * m_resolution.y();

    //                 // test: align all to texture layout
    //                 size_t idx = x
    //                              + y * m_resolution.x() 
    //                              + z * m_resolution.x() * m_resolution.y();

    //                 // Log(Debug, "[SuperGrid] idx: %f",idx);
    //                 // Log(Debug, "[SuperGrid] min, maj: %f, %f",min,maj);

    //                 if constexpr (!dr::is_jit_v<Float>) {
    //                     dr::scatter(m_extremum_grid, scale * Vector2f(min, maj), UInt32(idx));
    //                 } else {
    //                     extremums[idx*2] = min;
    //                     extremums[idx*2+1] = maj;
    //                 }
    //             }
    //         }
    //     }
    //     // Log(Debug, "exetrumum_grid: %f", m_extremum_grid);

    //     if constexpr (dr::is_jit_v<Float>) {
    //         m_extremum_grid = dr::load<FloatStorage>(extremums.data(), n*2);
    //     }
    //     Log(Info, "Extremum grid constructed successfully");
    // }

    void build_grid(const Volume *sigma_t, ScalarFloat scale) {

        m_bbox = sigma_t->bbox();
        // local space supergrid cell size
        const ScalarVector3f cell_size = dr::rcp(ScalarVector3f(m_resolution));
        
        // Log(Debug, "Building grid =======");
        // Log(Debug, "m_bbox: %f, cell_size: %f", m_bbox, cell_size);

        // Allocate extremum grid data
        size_t n = dr::prod(m_resolution);
        
        std::vector<ScalarFloat> extremums;
        // std::unique_ptr<ScalarFloat[]> extremums(new ScalarFloat[n*2]);

        // IMPORTANT: need to retrieve the data before to avoid ref count slow-down
        auto data = sigma_t->array();

        size_t n_threads = pool_size() + 1;
        size_t grain_size = std::max( n / (4 * n_threads), (size_t) 1 );

        if (!dr::is_jit_v<Float> || n < n_threads ) {
        // if constexpr (!dr::is_jit_v<Float> ) {

            // avoid unecessary allocation in scalar mode
            if constexpr (!dr::is_jit_v<Float>) {
                m_extremum_grid = dr::empty<FloatStorage>(n*2);
            } else {
                extremums.resize(n*2);
            }
            
            dr::parallel_for(
                dr::blocked_range<size_t>(0, n, grain_size),
                [&](const dr::blocked_range<size_t> &range) {
                    // Recover x, y, z from block start (one-time div/mod per block)

                    // Log(Info, "begin: %f, end: %f", range.begin(), range.end());
                    for (auto idx = range.begin(); idx != range.end(); ++idx) {
                        // Store in linear array (Z-slowest, X-fastest)
                        int32_t x = idx % m_resolution.x() ; 
                        int32_t y = (idx / m_resolution.x())  % m_resolution.y(); 
                        int32_t z =  idx / (m_resolution.x() * m_resolution.y()); 

                        
                        ScalarPoint3f cell_min = ScalarVector3f(x, y, z) * cell_size;
                        ScalarPoint3f cell_max = cell_min + cell_size;
                        ScalarBoundingBox3f cell_bounds(
                            cell_min + math::RayEpsilon<Float>, 
                            cell_max - math::RayEpsilon<Float>
                        );

                        // // Query volume for local extremum, currently assume local bounds.
                        auto [maj, min] = sigma_t->extremum(data, cell_bounds);

                        if constexpr (!dr::is_jit_v<Float>) {
                            dr::scatter(m_extremum_grid, scale * Vector2f(min, maj), UInt32(idx));
                        } else {
                            extremums[idx*2] = min[0];
                            extremums[idx*2+1] = maj[0];
                        }
                    }
                }
            );
        
            if constexpr (dr::is_jit_v<Float>) {
                m_extremum_grid = dr::load<FloatStorage>(extremums.data(), n*2);
            }
        } else {
            
            m_extremum_grid = dr::empty<FloatStorage>(n*2);

            UInt32 idx = dr::arange<UInt32>((uint32_t) n);

            UInt32 x = idx % m_resolution.x() ; 
            UInt32 y = (idx / m_resolution.x())  % m_resolution.y(); 
            UInt32 z =  idx / (m_resolution.x() * m_resolution.y()); 

            Point3f cell_min = Vector3f(x, y, z) * cell_size;
            Point3f cell_max = cell_min + cell_size;
            BoundingBox3f cell_bounds(
                            cell_min + math::RayEpsilon<Float>, 
                            cell_max - math::RayEpsilon<Float>
                        );

            auto [maj, min] = sigma_t->extremum(data, cell_bounds);

            dr::scatter(m_extremum_grid, min, idx*2);
            dr::scatter(m_extremum_grid, maj, idx*2+1);
            // Log(Debug, "min, maj: %f, %f", min, maj);
            // Log(Debug, "cell_bounds: %f", cell_bounds);
            dr::sync_thread();
        }

        
        Log(Info, "Extremum grid constructed successfully");
    }

    /**
     * \brief Prepare DDA traversal state for a ray
     *
     * This method initializes the DDA (Digital Differential Analyzer) algorithm
     * state for traversing the extremum grid along a ray.
     *
     * \return Tuple of (dda_t, dda_tmax, dda_tdelta, valid)
     *         where valid indicates if the ray intersects the grid
     */
    std::tuple<Float, Vector3f, Vector3f, Mask>
    prepare_dda_traversal(const Ray3f &ray, Float mint, Float maxt, Mask active) const {
        
        Log(Debug, "prepare ray traversals");
        //TODO: account for intersection with volume geometry

        // Transform ray to local grid coordinates [0,1]³
        const ScalarVector3f extents = m_bbox.extents();
        Ray3f local_ray(
            (ray.o - m_bbox.min) / extents,  // Normalize origin
            ray.d / extents,                  // Normalize direction
            ray.time,
            ray.wavelengths
        );
        const Vector3f local_voxel_size = 1.f / Vector3f(m_resolution);

        // Compute current and last voxel indices
        Vector3i current_voxel = dr::floor(local_ray(mint) / local_voxel_size);
        Vector3i last_voxel = dr::floor(local_ray(maxt) / local_voxel_size);

        current_voxel = dr::clip(current_voxel, 0, m_resolution - 1);
        last_voxel = dr::clip(last_voxel, 0, m_resolution - 1);

        Log(Debug, "local_ray: %f, local_voxel_size: %f", local_ray, local_voxel_size);
        Log(Debug, "current_voxel: %f, last_voxel: %f", current_voxel, last_voxel);

        // Traversal direction
        Vector3i step = dr::select(local_ray.d >= 0.f, 1, -1);

        Vector3f next_voxel_boundary =
            Vector3f(current_voxel + step) * local_voxel_size;

        // Handle negative directions
        next_voxel_boundary += dr::select(
            (current_voxel != last_voxel) && (local_ray.d < 0.f),
            local_voxel_size,
            0.f
        );

        Log(Debug, "step: %f, next_voxel_boundary: %f", step, next_voxel_boundary);

        // Compute DDA parameters
        auto ray_nonzero = local_ray.d != 0.f;
        // Value of ray parameter until next intersection with voxel-border along each axis
        Vector3f dda_tmax = dr::select(
            ray_nonzero,
            (next_voxel_boundary - local_ray.o) / local_ray.d,
            dr::Infinity<Float>
        );

        // How far along each component of the ray we must move to move by one voxel
        Vector3f dda_tdelta = dr::select(
            ray_nonzero,
            Vector3f(step) * local_voxel_size / local_ray.d,
            dr::Infinity<Float>
        );

        Log(Debug, "dr::isfinite(dda_tmax): %f", dda_tmax);
        Log(Debug, "dr::isfinite(dda_tdelta): %f", dda_tdelta);
        Log(Debug, "valid %f", dr::all(dr::isfinite(dda_tmax) || dr::isfinite(dda_tdelta)));

        // Current ray parameter throughout DDA traversal
        Float dda_t = mint;
        // Mask valid = active && dr::all(dr::isfinite(dda_tmax) || dr::isfinite(dda_tdelta));
        Mask valid = active;

        return { dda_t, dda_tmax, dda_tdelta, valid };
    }

private:
    /// Grid storing pre-computed local majorants
    // ref<Volume> m_extremum_grid;
    FloatStorage m_extremum_grid;
    /// Safety factor multiplied with all majorants
    ScalarFloat m_safety_factor;

    ScalarVector3i m_resolution;
    ScalarAffineTransform4f m_to_local;
};

MI_EXPORT_PLUGIN(ExtremumGrid)
NAMESPACE_END(mitsuba)
