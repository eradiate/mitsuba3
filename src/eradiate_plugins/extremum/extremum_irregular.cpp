#include <mitsuba/core/properties.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/eradiate/extremum.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/volumegrid.h>
#include <nanothread/nanothread.h>

#include "extremum_build_common.h"

NAMESPACE_BEGIN(mitsuba)

/**!
.. _extremum-extremum_irregular:

Irregular extremum grid structure (:monosp:`extremum_irregular`)
----------------------------------------------------------------

.. pluginparameters::

 * - volume
   - |volume|
   - Extinction coefficient volume to build the irregular extremum grid from.
   - |exposed|

 * - to_world
   - |transform|
   - Specifies the 4×4 transformation matrix of the underlying volume.

 * - scale
   - |float|
   - Scale factor for the extremum values. Default: 1.0

 * - resolution
   - |vector|
   - Virtual fine grid resolution along the XYZ axis. Set to [0,0,0] to use
     the resolution of the underlying volume. Default: [0,0,0]

 * - wrap_mode
   - |string|
   - Wrapping policy applied outside of the volume bounds. One of
     :monosp:`clamp`, :monosp:`repeat`, or :monosp:`mirror`.
     Default: :monosp:`clamp`

 * - merge_threshold
   - |float|
   - Relative extremum-span tolerance used when merging fine cells into
     boxes. A value of 0 only merges cells with identical extrema.
     Default: 0.05

 * - rebuild_threshold
   - |float|
   - Lazy topology rebuild policy applied on volume updates. Negative values
     disable topology rebuilds (per-box extrema are still refit). A positive
     value ``r`` triggers a rebuild when the traversal quality degrades by
     more than a factor ``1 + r`` relative to the last build. Default: -1.0

This plugin creates an irregular grid structure storing local extremum values
for efficient delta tracking in heterogeneous media, following the approach
of Pérard-Gayot et al., *GPU Ray Tracing using Irregular Grids* (2017).

A virtual fine grid is built over the volume, and neighboring cells with
similar extrema are merged into axis-aligned boxes. Each fine cell stores the
index of the box covering it, forming an exact partition of the volume. At
runtime, traversal proceeds like a DDA but advances an entire box per step:
homogeneous regions (*e.g.* clear sky) are crossed in large steps while
heterogeneous regions (*e.g.* clouds) retain fine-grained extrema.

When the underlying volume is updated (*e.g.* at each spectral iteration),
the box extrema are refit in place without changing the topology, which is
considerably cheaper than a full rebuild. The exposed parameters are
recomputed by this refit, so externally written values are overwritten on
the next volume update.
*/

template <typename Float, typename Spectrum>
class ExtremumIrregularGrid final : public ExtremumStructure<Float, Spectrum> {
public:
    MI_IMPORT_BASE(ExtremumStructure, m_bbox)
    MI_IMPORT_TYPES(Volume)

    using TrackingStateType    = TrackingState<Float, Spectrum>;
    using TrackingFunctionType = TrackingFunction<Float, Spectrum>;
    using FloatStorage         = DynamicBuffer<Float>;
    using UInt32Storage        = DynamicBuffer<UInt32>;
    using Int32Storage         = DynamicBuffer<Int32>;
    using Mask3                = dr::mask_t<Vector3i>;

    ExtremumIrregularGrid(const Properties &props) : Base(props) {
        // Volume Parameters
        m_volume = nullptr;
        for (auto &prop : props.objects()) {
            if (auto *vol = prop.try_get<Volume>()) {
                m_volume = vol;
                break;
            }
        }

        if (!m_volume)
            Throw("ExtremumIrregularGrid requires at least one volume");

        // Register the extremum structure to the volume to trigger
        // parameter_changed when the volume is modified.
        m_volume->add_extremum_structure(this);
        m_bbox = m_volume->bbox();

        m_to_local = props.get<ScalarAffineTransform4f>("to_world", ScalarAffineTransform4f()).inverse();
        m_scale = props.get<ScalarFloat>("scale", 1.0f);

        std::string_view wrap_mode_str = props.get<std::string_view>("wrap_mode", "clamp");
        if (wrap_mode_str == "repeat")
            m_wrap_mode = dr::WrapMode::Repeat;
        else if (wrap_mode_str == "mirror")
            m_wrap_mode = dr::WrapMode::Mirror;
        else if (wrap_mode_str == "clamp")
            m_wrap_mode = dr::WrapMode::Clamp;
        else
            Throw("Invalid wrap_mode \"%s\", must be one of: \"repeat\", "
                  "\"mirror\", or \"clamp\"!", wrap_mode_str);

        // Resolution Parameters
        m_resolution = props.get<ScalarVector3i>("resolution", ScalarVector3i(0));
        if (dr::any(m_resolution <= ScalarVector3i(0)))
            m_resolution = m_volume->resolution();

        size_t n = (size_t) m_resolution.x() * (size_t) m_resolution.y() *
                   (size_t) m_resolution.z();
        if (n > (size_t) 256 * 256 * 256)
            Log(Warn, "ExtremumIrregularGrid: virtual fine grid resolution %s "
                      "requires a %zu MiB cell index, consider specifying a "
                      "coarser \"resolution\"",
                m_resolution, (n * sizeof(uint32_t)) >> 20);

        for (size_t i = 0; i < 3; ++i)
            m_inv_resolution[i] = dr::divisor<int32_t>((int32_t) m_resolution[i]);

        m_merge_threshold   = props.get<ScalarFloat>("merge_threshold", 0.05f);
        if (m_merge_threshold < 0.f)
            Throw("\"merge_threshold\" must be non-negative");
        m_rebuild_threshold = props.get<ScalarFloat>("rebuild_threshold", -1.f);

        m_cn = props.get<ScalarFloat>("cn", 1.f);
        m_ct = props.get<ScalarFloat>("ct", 1.f);

        build(m_volume.get());
    }

    void parameters_changed(const std::vector<std::string> &/*keys*/ = {}) override {
        refit(m_volume.get());
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
        return traverse_irregular(
            func,
            state,
            ray,
            mint,
            maxt,
            channel,
            active
        );
    }

    std::tuple<Float, Float> eval_1(const Interaction3f &it,
                                    Mask active) const override {
        Point3f po = m_to_local * it.p;
        auto [piw, tile, flip] =
            wrap_with_tile(dr::floor2int<Vector3i>(po * m_resolution));

        UInt32 idx = dr::fmadd(dr::fmadd(piw.z(), m_resolution.y(), piw.y()),
                               m_resolution.x(), piw.x());
        UInt32 box = dr::gather<UInt32>(m_index_grid, idx, active);
        Vector2f extremum = dr::gather<Vector2f>(m_box_extrema, box, active);
        return { extremum.x(), extremum.y() };
    }

    void traverse(TraversalCallback *cb) override {
        cb->put("index_grid", m_index_grid, ParamFlags::NonDifferentiable);
        cb->put("box_lo", m_box_lo, ParamFlags::NonDifferentiable);
        cb->put("box_hi", m_box_hi, ParamFlags::NonDifferentiable);
        cb->put("box_extrema", m_box_extrema, ParamFlags::NonDifferentiable);
        cb->put("resolution", m_resolution, ParamFlags::NonDifferentiable);
        Base::traverse(cb);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "ExtremumIrregularGrid[" << std::endl
            << "  resolution = " << m_resolution << "," << std::endl
            << "  box_count = " << m_box_count << "," << std::endl
            << "  quality_at_build = " << m_quality_at_build << "," << std::endl
            << "  bbox = " << m_bbox << "," << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS(ExtremumIrregularGrid)

private:

    /**
     * \brief Build the irregular grid from a volume.
     *
     * Computes a fine extremum grid and merges its cells into boxes.
     */
    void build(const Volume *volume) {
        FloatStorage fine = build_fine_extremum_grid<Float, Spectrum>(
            volume, m_resolution, m_scale);
        build_topology(fine);

        Log(Info, "Irregular extremum grid constructed successfully (%u boxes)",
            m_box_count);
    }

    /**
     * \brief Refit the box extrema to an updated volume, keeping the
     * topology.
     *
     * Recomputes the fine extremum grid and reduces it per box through the
     * cell index. This is the spectral update fast path: it runs at every
     * volume update, so unlike the topology build it is fully vectorized
     * in JIT modes. When a rebuild threshold is set and the traversal
     * quality has degraded past it, the topology is rebuilt from the
     * already-computed fine grid.
     */
    void refit(const Volume *volume) {
        FloatStorage fine = build_fine_extremum_grid<Float, Spectrum>(
            volume, m_resolution, m_scale);

        const size_t n = (size_t) m_resolution.x() * (size_t) m_resolution.y() *
                         (size_t) m_resolution.z();
        ScalarFloat quality = 1.f;

        if constexpr (dr::is_jit_v<Float>) {
            FloatStorage extrema = dr::zeros<FloatStorage>(2 * m_box_count);
            dr::scatter(extrema,
                        dr::full<Float>(dr::Infinity<ScalarFloat>, m_box_count),
                        dr::arange<UInt32>(m_box_count) * 2);

            UInt32 cell = dr::arange<UInt32>((uint32_t) n);
            UInt32 box  = dr::gather<UInt32>(m_index_grid, cell);
            Float fmin  = dr::gather<Float>(fine, cell * 2);
            Float fmaj  = dr::gather<Float>(fine, cell * 2 + 1);
            dr::scatter_reduce(ReduceOp::Min, extrema, fmin, box * 2);
            // Mixing reduction types on an unevaluated target is illegal
            // (dr.ReduceMode.Expand), so resolve the Min reduction first
            dr::eval(extrema);
            dr::scatter_reduce(ReduceOp::Max, extrema, fmaj, box * 2 + 1);
            m_box_extrema = extrema;

            if (m_rebuild_threshold >= 0.f) {
                // One host sync per volume update, only when the lazy
                // rebuild policy is enabled
                Float box_maj = dr::gather<Float>(extrema, box * 2 + 1);
                ScalarFloat sum_box  = dr::sum(box_maj)[0],
                            sum_cell = dr::sum(fmaj)[0];
                quality = sum_cell > 0.f ? sum_box / sum_cell : 1.f;
            }
        } else {
            const ScalarFloat *fine_data = fine.data();
            const uint32_t *index        = m_index_grid.data();

            // Serial reduce: a parallel version would race on the box slots,
            // and this pass is negligible next to the fine grid build
            std::vector<ScalarFloat> extrema(2 * (size_t) m_box_count);
            for (size_t b = 0; b < m_box_count; ++b) {
                extrema[b * 2 + 0] = dr::Infinity<ScalarFloat>;
                extrema[b * 2 + 1] = 0.f;
            }
            for (size_t cell = 0; cell < n; ++cell) {
                size_t b = index[cell];
                extrema[b * 2 + 0] = std::min(extrema[b * 2 + 0],
                                              fine_data[cell * 2 + 0]);
                extrema[b * 2 + 1] = std::max(extrema[b * 2 + 1],
                                              fine_data[cell * 2 + 1]);
            }
            m_box_extrema = dr::load<FloatStorage>(extrema.data(),
                                                   extrema.size());

            if (m_rebuild_threshold >= 0.f)
                quality = compute_quality(fine_data, index, extrema.data(), n);
        }

        if (m_rebuild_threshold >= 0.f &&
            quality > (1.f + m_rebuild_threshold) * m_quality_at_build) {
            Log(Debug, "Irregular extremum grid quality degraded (%f -> %f), "
                       "rebuilding topology", m_quality_at_build, quality);
            build_topology(fine);
        }
    }

    /**
     * \brief Merge the fine grid cells into boxes and build the cell index.
     *
     * Greedy scan-order merge: for each unassigned cell, a box is grown
     * along +x (cell strip), then +y (row slabs), then +z (plane slabs),
     * as long as the extended cells are unassigned and the merged extrema
     * remain within the relative span tolerance. The result is an exact
     * partition of the fine grid.
     *
     * Runs on the host: topology builds are rare (refits do not change the
     * topology), so a scalar implementation is sufficient. Note that reading
     * ``fine.data()`` on the host assumes CPU-resident storage; a CUDA
     * backend would require a ``dr::migrate`` to host memory first.
     */
    void build_topology(const FloatStorage &fine) {
        if constexpr (dr::is_jit_v<Float>) {
            dr::eval(fine);
            dr::sync_thread();
        }
        const auto *fine_data = fine.data();

        const int32_t rx = m_resolution.x(),
                      ry = m_resolution.y(),
                      rz = m_resolution.z();
        const size_t n = (size_t) rx * (size_t) ry * (size_t) rz;

        auto cell_of = [&](int32_t x, int32_t y, int32_t z) {
            return ((size_t) z * ry + y) * rx + x;
        };

        // Aggregates over a candidate box: the box extrema (min_min,
        // maj_max) plus the spans needed by the merge criterion
        struct Aggregate {
            ScalarFloat min_min, min_max, maj_min, maj_max;
            ScalarVector3f lo, hi;
            ScalarBoundingBox3f bbox;

            void merge_cell(const ScalarFloat *data, size_t cell, ScalarVector3f new_hi) {
                ScalarFloat mn = data[cell * 2 + 0],
                            mj = data[cell * 2 + 1];
                min_min = std::min(min_min, mn);
                min_max = std::max(min_max, mn);
                maj_min = std::min(maj_min, mj);
                maj_max = std::max(maj_max, mj);
                hi = new_hi;
                bbox.expand(ScalarPoint3f(new_hi));
            }

            void merge_agg(const Aggregate& agg) {
                min_min = std::min(min_min, agg.min_min);
                min_max = std::max(min_max, agg.min_max);
                maj_min = std::min(maj_min, agg.maj_min);
                maj_max = std::max(maj_max, agg.maj_max);

                lo = dr::minimum(lo, agg.lo);
                hi = dr::maximum(hi, agg.hi);
                bbox.expand(agg.bbox);
            }
        };
        auto aggregate_cell = [&](size_t cell, ScalarVector3f lo, ScalarVector3f hi) {
            ScalarFloat mn = fine_data[cell * 2 + 0],
                        mj = fine_data[cell * 2 + 1];

            return Aggregate{ mn, mn, mj, mj, lo, hi, ScalarBoundingBox3f(lo,hi) };
        };

        // Both spans are normalized by the majorant: the majorant span
        // bounds the wasted null collisions, the minorant span bounds the
        // residual magnitude for residual ratio tracking. Normalizing the
        // minorant span by the minorant itself would diverge in empty
        // regions.
        const ScalarFloat thr = m_merge_threshold;
        // auto criterion = [thr](const Aggregate &a) {
        //     return (a.maj_max - a.maj_min) <= thr * a.maj_max &&
        //            (a.min_max - a.min_min) <= thr * a.maj_max;
        // };

        // auto criterion = [thr](const Aggregate &a) {
        //     return (a.maj_max - a.min_min) <= thr * a.maj_max;
        // };

        auto criterion = [thr](const Aggregate &/*a*/, const Aggregate &/*a*/, const Aggregate &a) {
            return (a.maj_max - a.min_min) * dr::norm(a.hi - a.lo) < thr;
        };

        ScalarVector3f cell_size =  m_bbox.extents() * dr::rcp(ScalarVector3f(m_resolution));
        ScalarFloat L = dr::mean(cell_size);
        ScalarFloat sigma_ref = m_volume->max();

        // auto cost = [this](const Aggregate &a) {
        //     return m_cn * a.bbox.volume() *  (a.maj_max; - a.min_min) /
        //            + m_ct * a.bbox.surface_area();
        // };
        // auto cost = [this, L, sigma_ref](const Aggregate &a) {
        //     ScalarFloat null = a.bbox.volume() * (a.maj_max - a.min_min) / sigma_ref;
        //     ScalarFloat trav = a.bbox.surface_area() * L * m_ct;
        //     return null + trav;
        // };

        // auto criterion = [this, cost](const Aggregate &a0, const Aggregate &a1, const Aggregate &am){
        //     bool veto = (am.maj_max - am.min_min) > m_cn                     // absolute, or
        //                 || am.maj_max > (1.f + dr::Epsilon<ScalarFloat>) * std::max(a0.maj_max, a1.maj_max); // relative

        //     return !veto && (cost(am) < cost(a0) + cost (a1));
        // };


        const uint32_t unassigned = 0xFFFFFFFFu;
        std::vector<uint32_t> index(n, unassigned);
        std::vector<int32_t> lo, hi;
        std::vector<ScalarFloat> extrema;


        for (int32_t z0 = 0; z0 < rz; ++z0)
        for (int32_t y0 = 0; y0 < ry; ++y0)
        for (int32_t x0 = 0; x0 < rx; ++x0) {
            if (index[cell_of(x0, y0, z0)] != unassigned)
                continue;

            ScalarVector3f b_lo = ScalarVector3f(x0, y0, z0) * cell_size;
            ScalarVector3f b_hi = b_lo + cell_size;

            Aggregate agg = aggregate_cell(cell_of(x0, y0, z0), b_lo, b_hi);

            // Grow a strip along +x
            int32_t x1 = x0 + 1;
            for (; x1 < rx; ++x1) {
                size_t cell = cell_of(x1, y0, z0);
                ScalarVector3f b_lo = ScalarVector3f(x1, y0, z0) * cell_size;
                ScalarVector3f b_hi = b_lo + cell_size;
                Aggregate new_agg = aggregate_cell(cell_of(x1, y0, z0), b_lo, b_hi);

                Aggregate cand = agg;
                ScalarVector3f new_hi = (ScalarVector3f(x1, y0, z0)+1.f) * cell_size;
                cand.merge_cell(fine_data, cell, new_hi);
                if (index[cell] != unassigned || !criterion(agg, new_agg, cand))
                    break;
                agg = cand;
            }

            // Grow row slabs along +y
            int32_t y1 = y0 + 1;
            for (; y1 < ry; ++y1) {
                Aggregate cand = agg;
                ScalarVector3f s_lo = ScalarVector3f(x0, y1, z0) * cell_size;
                ScalarVector3f s_hi = s_lo + cell_size;
                Aggregate new_agg = aggregate_cell(cell_of(x0, y1, z0), s_lo, s_hi);

                bool accept = true;
                for (int32_t x = x0; x < x1 && accept; ++x) {
                    size_t cell = cell_of(x, y1, z0);
                    ScalarVector3f new_hi = (ScalarVector3f(x, y1, z0)+1.f) * cell_size;
                    accept &= index[cell] == unassigned;
                    new_agg.merge_cell(fine_data, cell, new_hi);
                }
                cand.merge_agg(new_agg);
                if (!accept || !criterion(agg, new_agg, cand))
                    break;
                agg = cand;
            }

            // Grow plane slabs along +z
            int32_t z1 = z0 + 1;
            for (; z1 < rz; ++z1) {
                Aggregate cand = agg;
                ScalarVector3f s_lo = ScalarVector3f(x0, y0, z1) * cell_size;
                ScalarVector3f s_hi = s_lo + cell_size;
                Aggregate new_agg = aggregate_cell(cell_of(x0, y0, z1), s_lo, s_hi);

                bool accept = true;
                for (int32_t y = y0; y < y1 && accept; ++y)
                    for (int32_t x = x0; x < x1 && accept; ++x) {
                        size_t cell = cell_of(x, y, z1);
                        accept &= index[cell] == unassigned;
                        ScalarVector3f new_hi = (ScalarVector3f(x, y, z1)+1.f) * cell_size;
                        new_agg.merge_cell(fine_data, cell, new_hi);
                    }
                cand.merge_agg(new_agg);
                if (!accept || !criterion(agg, new_agg, cand))
                    break;
                agg = cand;
            }

            uint32_t box_id = (uint32_t) (lo.size() / 3);
            for (int32_t z = z0; z < z1; ++z)
                for (int32_t y = y0; y < y1; ++y)
                    for (int32_t x = x0; x < x1; ++x)
                        index[cell_of(x, y, z)] = box_id;

            lo.insert(lo.end(), { x0, y0, z0 });
            hi.insert(hi.end(), { x1, y1, z1 });
            extrema.insert(extrema.end(), { agg.min_min, agg.maj_max });
        }

        m_quality_at_build = compute_quality(fine_data, index.data(),
                                             extrema.data(), n);

        m_box_count   = (uint32_t) (lo.size() / 3);
        m_index_grid  = dr::load<UInt32Storage>(index.data(), n);
        m_box_lo      = dr::load<Int32Storage>(lo.data(), lo.size());
        m_box_hi      = dr::load<Int32Storage>(hi.data(), hi.size());
        m_box_extrema = dr::load<FloatStorage>(extrema.data(), extrema.size());
    }

    /**
     * \brief Traversal quality metric of a cell-to-box assignment.
     *
     * Ratio of the per-cell sums of box majorants over fine cell majorants
     * (>= 1; cell volumes are uniform and cancel out). A value of 1 means
     * the boxes are as tight as the fine grid; quality degrades as the
     * underlying volume drifts away from the topology it was built for.
     */
    ScalarFloat compute_quality(const ScalarFloat *fine_data,
                                const uint32_t *index,
                                const ScalarFloat *extrema, size_t n) const {
        double sum_box = 0.0, sum_cell = 0.0;
        for (size_t cell = 0; cell < n; ++cell) {
            sum_cell += (double) fine_data[cell * 2 + 1];
            sum_box  += (double) extrema[(size_t) index[cell] * 2 + 1];
        }
        return sum_cell > 0.0 ? (ScalarFloat) (sum_box / sum_cell) : 1.f;
    }

    /** \brief Irregular grid traversal algorithm.
     *
     * This method traverses the irregular grid along the provided ray. Each
     * step looks up the box covering the current cell and advances directly
     * to the box exit, calling ``func`` with the resulting segment. See
     * \ref ExtremumGrid::traverse_dda for the callback contract; the segment
     * distance conventions are identical.
     */
    template<typename FuncT, typename StateT>
    std::decay_t<StateT> traverse_irregular(
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

        // Transform ray to local grid coordinates [0,res]³
        Vector3f res = Vector3f(m_resolution);
        Ray3f local_ray((m_to_local * ray.o) * res,
                        (m_to_local * ray.d) * res,
                        ray.time, ray.wavelengths);
        Vector3f rcp_d = dr::rcp(local_ray.d);
        auto inf_t     = local_ray.d == 0.f;
        auto d_pos     = local_ray.d >= 0.f;

        Float t_min = mint;
        Float t_max = maxt;

        active &= t_max > t_min;

        // Advance the ray to the start of the interval
        local_ray.o = dr::fmadd(local_ray.d, t_min, local_ray.o);
        t_max       = t_max - t_min;
        t_min       = 0.f;

        // Integer grid coordinates, in unwrapped space (may lie outside
        // [0,res) under repeat/mirror policies)
        Vector3i pi = dr::floor2int<Vector3i>(local_ray.o);

        Vector3f local_o = local_ray.o;
        Vector3f local_d = local_ray.d;

        struct LoopState {
            ExtremumSegment segment;
            StateD state;
            Mask advance;
            Mask active;
            Vector3i pi;
            Float t_rem;

            DRJIT_STRUCT(LoopState, segment, state, advance, active, pi, \
                t_rem)
        } ls = {
            segment,
            state,
            /*advance=*/active,
            active,
            pi,
            t_max
        };

        dr::tie(ls) = dr::while_loop(
            dr::make_tuple(ls),
            [](const LoopState& ls) { return ls.active; },
            [this, func, local_o, local_d, rcp_d, inf_t, d_pos, t_max, mint,
             channel](LoopState& ls) {

            ExtremumSegment& segment = ls.segment;
            StateD& state   = ls.state;
            Mask& advance   = ls.advance;
            Mask& active    = ls.active;
            Vector3i& pi    = ls.pi;
            Float& t_rem    = ls.t_rem;

            Vector3i res_i = Vector3i(m_resolution);

            // Look up the box covering the current cell
            auto [piw, tile, flip] = wrap_with_tile(pi);
            UInt32 cell = dr::fmadd(dr::fmadd(piw.z(), m_resolution.y(), piw.y()),
                                    m_resolution.x(), piw.x());

            UInt32 box        = dr::gather<UInt32>(m_index_grid, cell, active);
            Vector3i lo       = dr::gather<Vector3i>(m_box_lo, box, active);
            Vector3i hi       = dr::gather<Vector3i>(m_box_hi, box, active);
            Vector2f extremum = dr::gather<Vector2f>(m_box_extrema, box, active);

            // Unwrap the canonical box bounds into the ray's current tile.
            // Mirrored tiles map local cells c -> res-1-c, so the half-open
            // cell range [lo,hi) reflects to [res-hi, res-lo).
            Vector3i blo_i, bhi_i;
            if (m_wrap_mode == dr::WrapMode::Mirror) {
                blo_i = dr::select(flip, res_i - hi, lo);
                bhi_i = dr::select(flip, res_i - lo, hi);
            } else {
                blo_i = lo;
                bhi_i = hi;
            }
            blo_i += tile * res_i;
            bhi_i += tile * res_i;

            Vector3f blo = Vector3f(blo_i);
            Vector3f bhi = Vector3f(bhi_i);
            Vector3i clip_lo = blo_i;
            Vector3i clip_hi = bhi_i;

            if (m_wrap_mode == dr::WrapMode::Clamp) {
                // The clamped field is constant outside the volume, so
                // boundary boxes extend to infinity along out-of-domain
                // axes. Sentinel integer bounds keep the cell clip below
                // consistent with the widened slab planes.
                const int32_t sentinel = 0x40000000;
                auto below = pi < 0;
                auto above = pi >= res_i;
                dr::masked(blo, below)     = -dr::Infinity<Float>;
                dr::masked(bhi, above)     = dr::Infinity<Float>;
                dr::masked(clip_lo, below) = -sentinel;
                dr::masked(clip_hi, above) = sentinel;
            }

            // Distance to the box exit, computed from the absolute local
            // origin to avoid accumulating floating point drift
            Float t_acc       = t_max - t_rem;
            Vector3f plane    = dr::select(d_pos, bhi, blo);
            Vector3f t_exit_v = dr::fmadd(plane - local_o, rcp_d, -t_acc);
            dr::masked(t_exit_v, inf_t) = dr::Infinity<Float>;

            Float dt_box = dr::min(t_exit_v);
            Float dt     = dr::maximum(dr::minimum(dt_box, t_rem), 0.f);
            auto exit_axis = t_exit_v == dt_box;

            // Store segment for lanes that reached target
            Float t_curr = mint + t_acc;
            dr::masked(segment, active) = ExtremumSegment(
                t_curr,
                t_curr + dt,
                extremum
            );

            std::tie(advance, active) = func(segment, state, channel, active);

            // Advance to the first cell beyond the box. On exit axes the
            // next cell is derived from the integer box bounds, which
            // guarantees progress; on the remaining axes the floored exit
            // position is clipped into the box extent to remove drift.
            Vector3i pi_exit  = dr::select(d_pos, bhi_i, blo_i - 1);
            Vector3f p_exit   = dr::fmadd(local_d, t_acc + dt, local_o);
            Vector3i pi_floor = dr::clip(dr::floor2int<Vector3i>(p_exit),
                                         clip_lo, clip_hi - 1);

            dr::masked(pi, advance)    = dr::select(exit_axis, pi_exit, pi_floor);
            dr::masked(t_rem, advance) -= dt;

            active &= t_rem > 0.f;
        },
        "Irregular Grid Traversal");

        return ls.state;
    };

    /**
     * \brief Applies the configured texture wrapping mode to an integer
     * position.
     *
     * \return
     *     A tuple ``(wrapped, tile, flip)`` holding the canonical cell
     *     coordinates in ``[0,res)``, the tile index of the position
     *     (floor division by the resolution; 0 in clamp mode), and a
     *     per-axis mirror flip mask (tiles with odd index are reflected).
     */
    std::tuple<Vector3i, Vector3i, Mask3> wrap_with_tile(const Vector3i &pos) const {
        Vector3i res_i = Vector3i(m_resolution);

        if (m_wrap_mode == dr::WrapMode::Clamp) {
            return { dr::clip(pos, 0, res_i - 1), Vector3i(0), Mask3(false) };
        } else {
            Vector3i value_shift_neg = dr::select(pos < 0, pos + 1, pos);

            Vector3i div;
            for (size_t i = 0; i < 3; ++i)
                div[i] = m_inv_resolution[i](value_shift_neg[i]);

            Vector3i mod = pos - div * res_i;
            mod[mod < 0] += res_i;

            // True floor-division tile: pos - mod is an exact multiple of
            // the resolution, so the truncated division is exact
            Vector3i tile;
            for (size_t i = 0; i < 3; ++i)
                tile[i] = m_inv_resolution[i](pos[i] - mod[i]);

            Mask3 flip = (tile & 1) != 0;

            if (m_wrap_mode == dr::WrapMode::Mirror)
                mod = dr::select(flip, res_i - 1 - mod, mod);

            return { mod, tile, flip };
        }
    }

private:
    /// Per-cell box index of the virtual fine grid
    UInt32Storage m_index_grid;
    /// Inclusive integer lower cell bounds of each box
    Int32Storage m_box_lo;
    /// Exclusive integer upper cell bounds of each box
    Int32Storage m_box_hi;
    /// Box extrema stored as interleaved [minorant, majorant] pairs
    FloatStorage m_box_extrema;
    uint32_t m_box_count;

    ref<Volume> m_volume;
    ScalarFloat m_scale;

    ScalarVector3i m_resolution;
    dr::divisor<int32_t> m_inv_resolution[3] { };
    dr::WrapMode m_wrap_mode;

    ScalarFloat m_merge_threshold;
    ScalarFloat m_rebuild_threshold;
    ScalarFloat m_quality_at_build = 1.f;

    ScalarAffineTransform4f m_to_local;

    ScalarFloat m_cn;
    ScalarFloat m_ct;
};

MI_EXPORT_PLUGIN(ExtremumIrregularGrid)
NAMESPACE_END(mitsuba)
