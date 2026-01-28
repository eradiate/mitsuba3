#include <mitsuba/core/properties.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/eradiate/extremum.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/volumegrid.h>

NAMESPACE_BEGIN(mitsuba)

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
        }

        // Build extremum grid
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

        // Prepare DDA traversal
        Float dda_t; 
        Vector3f dda_tmax, dda_tdelta;
        Mask valid;
        std::tie(dda_t, dda_tmax, dda_tdelta, valid) =
            prepare_dda_traversal(ray, mint, maxt, active);

        active &= valid;

        // Traverse grid with DDA until desired_tau is reached
        Float tau_acc = 0.f;
        Mask reached = false;

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

            // Find next voxel boundary (minimum of dda_tmax components)
            Float t_next = dr::minimum(dr::minimum(dda_tmax.x(), dda_tmax.y()), dda_tmax.z());
            t_next = dr::minimum(t_next, maxt);

            // Determine which axis advances
            Vector3f tmax_update = dr::zeros<Vector3f>();
            Mask assigned = false;
            // TODO: what about if we start in the middle of the medium?
            for (size_t k = 0; k < 3; ++k) {
                Mask hit_axis = (dda_tmax[k] == t_next) && !assigned;
                tmax_update[k] = dr::select(hit_axis, dda_tdelta[k], 0.f);
                assigned |= hit_axis;
            }

            // Lookup local majorant at cell center
            Float t_mid = 0.5f * (dda_t + t_next);
            // MediumInteraction3f lookup_mei = dr::zeros<MediumInteraction3f>();
            // lookup_mei.t = t_mid;
            // lookup_mei.p = ray(t_mid);

            // Retrieve index from cell center position.
            const Vector3f pos = (m_to_local * ray(t_mid)) * Vector3f(m_resolution);
            const Vector3i pos_i = dr::floor2int<Vector3i>(pos);
            UInt32 index( 2 * dr::fmadd( 
                        dr::fmadd( UInt32(pos_i.z()), m_resolution.y(), UInt32(pos_i.y()) ), 
                        m_resolution.x(), 
                        UInt32(pos_i.x())
                    ));

            dr::Array<Float, 2> local_extremum = dr::gather<dr::Array<Float, 2>>(m_extremum_grid, index);
            const Float local_minorant = local_extremum[0];
            const Float local_majorant = local_extremum[1];                    
            // const Float local_majorant = m_extremum_grid->eval_1(lookup_mei, active);

            // Accumulate optical depth
            Float segment_length = t_next - dda_t;
            Float tau_next = tau_acc + local_majorant * segment_length;

            // Check if desired tau reached in this segment
            Mask stops_here = active && (local_majorant > 0.f) &&
                             (tau_next >= desired_tau) && (t_next <= maxt);

            // Compute precise intersection point
            Float t_precise = dda_t +
                (desired_tau - tau_acc) / dr::maximum(local_majorant, dr::Epsilon<Float>);

            // Update reached status
            reached |= stops_here;

            // Store result for lanes that reached target
            dr::masked(result.tmin, stops_here) = t_precise;
            dr::masked(result.tmax, stops_here) = t_next;
            dr::masked(result.sigma_min, stops_here) = local_minorant;
            dr::masked(result.sigma_maj, stops_here) = local_majorant;

            // Advance DDA state
            dr::masked(dda_t, active && !reached) = t_next;
            dr::masked(dda_tmax, active && !reached) = dda_tmax + tmax_update;
            dr::masked(tau_acc, active && !reached) = tau_next;

            // Continue only if not reached and still in bounds
            active &= !reached && (t_next < maxt);
        },
        "DDA Travesal");

        // For lanes that didn't reach, set tmin to infinity
        dr::masked(ls.result.tmin, !reached) = dr::Infinity<Float>;

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
    void build_grid(const Volume *sigma_t, ScalarFloat scale) {
        m_bbox = sigma_t->bbox();

        

        const ScalarVector3f cell_size = m_bbox.extents() / ScalarVector3f(m_resolution);

        // Allocate extremum grid data
        size_t n = dr::prod(m_resolution);
        std::vector<ScalarFloat> extremums(n*2);

        // Build grid (consider parallelizing)
        for (int32_t z = 0; z < m_resolution.z(); ++z) {
            for (int32_t y = 0; y < m_resolution.y(); ++y) {
                for (int32_t x = 0; x < m_resolution.x(); ++x) {
                    // Compute cell bounding box in world space
                    ScalarPoint3f cell_min = m_bbox.min +
                        ScalarVector3f(x, y, z) * cell_size;
                    ScalarPoint3f cell_max = cell_min + cell_size;
                    ScalarBoundingBox3f cell_bounds(cell_min, cell_max);

                    // Query volume for local extremum
                    // !! The extremum function transforms the cell bound to 
                    // local reference system. It is currently limited to 
                    // axis-aligned queries. Need to find some way to include local 
                    // bound queries
                    auto [maj, min] = sigma_t->extremum(cell_bounds);
                    min = scale * min;
                    maj = scale * maj;

                    // Store in linear array (Z-slowest, X-fastest)
                    size_t idx = x * 2
                                 + y * 2 * m_resolution.x() 
                                 + z * 2 * m_resolution.x() * m_resolution.y();
                    extremums[idx] = min;
                    extremums[idx+1] = maj;
                }
            }
        }


        m_extremum_grid = dr::load<FloatStorage>(extremums.data(), n*2);
        // // Create grid volume from computed majorants
        // Properties grid_props("gridvolume");
        // grid_props.set("raw", true);  // No color conversion

        // // Create tensor from majorants
        // size_t shape[4] = {
        //     (size_t) m_resolution.z(),
        //     (size_t) m_resolution.y(),
        //     (size_t) m_resolution.x(),
        //     2  // Minorant and Majorant
        // };

        // TensorXf tensor(majorants.data(), 4, shape);
        // grid_props.set_any("data", tensor);

        // // Set transform to match sigma_t bbox
        // ScalarAffineTransform4f to_world = ScalarAffineTransform4f::scale(m_bbox.extents());
        // to_world = to_world * ScalarAffineTransform4f::translate(m_bbox.min);
        // grid_props.set("to_world", to_world);

        // // Disable filtering (use nearest neighbor for grid cells)
        // grid_props.set("filter_type", "nearest");

        // m_extremum_grid = PluginManager::instance()->create_object<Volume>(grid_props);

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

        // Clamp to valid range
        current_voxel = dr::clip(current_voxel, 0, m_resolution - 1);
        last_voxel = dr::clip(last_voxel, 0, m_resolution - 1);

        // Traversal direction
        Vector3i step = dr::select(local_ray.d >= 0.f, 1, -1);

        // Next voxel boundaries
        Vector3f next_voxel_boundary =
            Vector3f(current_voxel + step) * local_voxel_size;

        // Handle negative directions
        next_voxel_boundary += dr::select(
            (current_voxel != last_voxel) && (local_ray.d < 0.f),
            local_voxel_size,
            0.f
        );

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

        // Current ray parameter throughout DDA traversal
        Float dda_t = mint;
        Mask valid = active && dr::all(dr::isfinite(dda_tmax) || dr::isfinite(dda_tdelta));

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
