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
.. _extremum-extremum_irregular2:

Irregular extremum grid structure, paper-faithful merge (:monosp:`extremum_irregular2`)
---------------------------------------------------------------------------------------

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

 * - cn
   - |float|
   - Weight of the null-collision term in the merge cost. Default: 1.0

 * - ct
   - |float|
   - Weight of the traversal (surface area) term in the merge cost.
     Default: 1.0

 * - alpha
   - |float|
   - Merge-loop stopping factor in ``(0, 1)``. Axis passes are repeated as
     long as an iteration reduces the box count below ``alpha`` times the
     previous count. Default: 0.9

 * - rebuild_threshold
   - |float|
   - Lazy topology rebuild policy applied on volume updates. Negative values
     disable topology rebuilds (per-box extrema are still refit). A positive
     value ``r`` triggers a rebuild when the traversal quality degrades by
     more than a factor ``1 + r`` relative to the last build. Default: -1.0

This plugin is a second implementation of the irregular grid structure of
Pérard-Gayot et al., *GPU Ray Tracing using Irregular Grids* (2017). It shares
its runtime traversal and storage with :monosp:`extremum_irregular` but builds
its topology with the **iterative, axis-aligned merge** of the paper (and its
reference implementation, HaGrid) instead of a greedy scan-order merge.

Starting from one box per fine cell, the build repeats passes over the X, Y and
Z axes. In each pass every box independently considers merging with its
neighbour along the axis, accepting the merge when it lowers a surface-area
heuristic (SAH) cost adapted to extremum structures:

.. math::

   \mathrm{cost}(\text{box}) = c_n \, V \, (\mu_{\max} - \mu_{\min})
                             + c_t \, \mathrm{SA},

where :math:`V` and :math:`\mathrm{SA}` are the world-space volume and surface
area of the box, :math:`\mu_{\min}/\mu_{\max}` its minorant/majorant, and
:math:`c_n/c_t` the ``cn``/``ct`` weights. The first term is the null-collision
waste (the extremum analogue of the paper's per-primitive intersection cost);
the second is the traversal cost. A chain-breaking step keeps only every other
box along a merge chain, preventing runaway cascades and yielding better-shaped
boxes than the greedy variant. Passes repeat until the box count stops
shrinking (controlled by ``alpha``).

When the underlying volume is updated (*e.g.* at each spectral iteration), the
box extrema are refit in place without changing the topology, exactly as in
:monosp:`extremum_irregular`.
*/

template <typename Float, typename Spectrum>
class ExtremumIrregularGrid2 final : public ExtremumStructure<Float, Spectrum> {
public:
    MI_IMPORT_BASE(ExtremumStructure, m_bbox)
    MI_IMPORT_TYPES(Volume)

    using TrackingStateType    = TrackingState<Float, Spectrum>;
    using TrackingFunctionType = TrackingFunction<Float, Spectrum>;
    using FloatStorage         = DynamicBuffer<Float>;
    using UInt32Storage        = DynamicBuffer<UInt32>;
    using Int32Storage         = DynamicBuffer<Int32>;
    using Mask3                = dr::mask_t<Vector3i>;

    ExtremumIrregularGrid2(const Properties &props) : Base(props) {
        // Volume Parameters
        m_volume = nullptr;
        for (auto &prop : props.objects()) {
            if (auto *vol = prop.try_get<Volume>()) {
                m_volume = vol;
                break;
            }
        }

        if (!m_volume)
            Throw("ExtremumIrregularGrid2 requires at least one volume");

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
            Log(Warn, "ExtremumIrregularGrid2: virtual fine grid resolution %s "
                      "requires a %zu MiB cell index, consider specifying a "
                      "coarser \"resolution\"",
                m_resolution, (n * sizeof(uint32_t)) >> 20);

        for (size_t i = 0; i < 3; ++i)
            m_inv_resolution[i] = dr::divisor<int32_t>((int32_t) m_resolution[i]);

        // Merge cost weights and stopping factor
        m_cn = props.get<ScalarFloat>("cn", 1.f);
        m_ct = props.get<ScalarFloat>("ct", 1.f);
        m_alpha = props.get<ScalarFloat>("alpha", 0.9f);
        if (m_alpha <= 0.f || m_alpha >= 1.f)
            Throw("\"alpha\" must lie in the open interval (0, 1)");

        m_rebuild_threshold = props.get<ScalarFloat>("rebuild_threshold", -1.f);

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
        oss << "ExtremumIrregularGrid2[" << std::endl
            << "  resolution = " << m_resolution << "," << std::endl
            << "  box_count = " << m_box_count << "," << std::endl
            << "  quality_at_build = " << m_quality_at_build << "," << std::endl
            << "  bbox = " << m_bbox << "," << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS(ExtremumIrregularGrid2)

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

        Log(Info, "Irregular extremum grid (paper merge) constructed "
                  "successfully (%u boxes)", m_box_count);
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
     * \brief Axis-aligned box used during the merge build.
     *
     * ``lo``/``hi`` are integer fine-grid coordinates with ``hi`` exclusive
     * (matching the storage layout and HaGrid's exclusive ``Cell.max``);
     * ``mn``/``mj`` are the box minorant/majorant.
     */
    struct Box {
        int32_t lo[3], hi[3];
        ScalarFloat mn, mj;
    };

    /**
     * \brief Merge the fine grid cells into boxes and build the cell index.
     *
     * Paper-faithful merge (Pérard-Gayot et al. 2017 / HaGrid ``merge.cu``):
     * starting from one box per fine cell, iterate passes over the X, Y and Z
     * axes. Each pass lets every box independently merge with its neighbour
     * along the axis when a surface-area heuristic cost improves, then breaks
     * merge chains (keeping every other box) to avoid cascades. Passes repeat
     * until the box count stops shrinking (``alpha``). The result is an exact
     * partition of the fine grid.
     *
     * Runs on the host and is serial for now: topology builds are rare (refits
     * do not change the topology). Reading ``fine.data()`` on the host assumes
     * CPU-resident storage; a CUDA backend would require a ``dr::migrate`` to
     * host memory first.
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

        // X-fastest / Z-slowest linearization, matching the fine grid build
        auto cell_of = [&](int32_t x, int32_t y, int32_t z) -> size_t {
            return ((size_t) z * ry + y) * rx + x;
        };

        // World-space size of a fine cell; drives the SAH cost below
        ScalarVector3f cell_size =
            m_bbox.extents() * dr::rcp(ScalarVector3f(m_resolution));

        // SAH cost adapted to extremum structures: a null-collision waste
        // term (volume x extremum span) plus a traversal term (surface area)
        auto cost = [&](const Box &b) -> ScalarFloat {
            ScalarFloat ex = (b.hi[0] - b.lo[0]) * cell_size.x();
            ScalarFloat ey = (b.hi[1] - b.lo[1]) * cell_size.y();
            ScalarFloat ez = (b.hi[2] - b.lo[2]) * cell_size.z();
            ScalarFloat V  = ex * ey * ez;
            ScalarFloat SA = 2.f * (ex * ey + ey * ez + ex * ez);
            return m_cn * V * (b.mj - b.mn) + m_ct * SA;
        };

        // A box may merge with its axis neighbour only if they are adjacent
        // along the axis and share the exact same extent on the other two
        // axes; this keeps merged boxes rectangular.
        auto aligned = [](int axis, const Box &a, const Box &b) -> bool {
            int a1 = (axis + 1) % 3, a2 = (axis + 2) % 3;
            return a.hi[axis] == b.lo[axis] &&
                   a.lo[a1] == b.lo[a1] && a.lo[a2] == b.lo[a2] &&
                   a.hi[a1] == b.hi[a1] && a.hi[a2] == b.hi[a2];
        };

        auto union_box = [](const Box &a, const Box &b) -> Box {
            Box u;
            for (int k = 0; k < 3; ++k) {
                u.lo[k] = std::min(a.lo[k], b.lo[k]);
                u.hi[k] = std::max(a.hi[k], b.hi[k]);
            }
            u.mn = std::min(a.mn, b.mn);
            u.mj = std::max(a.mj, b.mj);
            return u;
        };

        // Initialize: one box per fine cell, index[c] = c
        std::vector<Box> boxes(n);
        std::vector<uint32_t> index(n);
        for (int32_t z = 0; z < rz; ++z)
            for (int32_t y = 0; y < ry; ++y)
                for (int32_t x = 0; x < rx; ++x) {
                    size_t c = cell_of(x, y, z);
                    Box &b = boxes[c];
                    b.lo[0] = x;     b.lo[1] = y;     b.lo[2] = z;
                    b.hi[0] = x + 1; b.hi[1] = y + 1; b.hi[2] = z + 1;
                    b.mn = fine_data[c * 2 + 0];
                    b.mj = fine_data[c * 2 + 1];
                    index[c] = (uint32_t) c;
                }

        // One merge pass along a single axis. ``empty_mask`` staggers which
        // box boundaries may merge this iteration (see the merge loop).
        auto merge_pass = [&](int axis, int32_t empty_mask) {
            const size_t nb = boxes.size();
            std::vector<int32_t> nexts(nb, -1), prevs(nb, -1);

            // Step 1: each box decides whether to merge with its +axis
            // neighbour when the SAH cost improves.
            for (size_t id = 0; id < nb; ++id) {
                const Box &b1 = boxes[id];

                // merge_allowed: on a flat dense grid HaGrid's octree shift is
                // 0, so is_top_level is always true and the gate reduces to
                // (pos & empty_mask) == 0. It staggers merge boundaries per
                // iteration and is unrestricted once empty_mask == 0.
                if ((b1.lo[axis] & empty_mask) != 0)
                    continue;

                int32_t np[3] = { b1.lo[0], b1.lo[1], b1.lo[2] };
                np[axis] = b1.hi[axis];
                if (np[axis] >= m_resolution[axis])
                    continue;

                uint32_t nid = index[cell_of(np[0], np[1], np[2])];
                const Box &b2 = boxes[nid];
                if (!aligned(axis, b1, b2))
                    continue;

                Box u = union_box(b1, b2);
                if (cost(u) <= cost(b1) + cost(b2)) {
                    nexts[id]  = (int32_t) nid;
                    prevs[nid] = (int32_t) id;
                }
            }

            // Step 2: break merge chains. Keep the chain head and every other
            // box (keep, discard, keep, ...), so each kept box merges with the
            // single discarded box that follows it.
            std::vector<uint8_t> flag(nb, 0);
            for (size_t id = 0; id < nb; ++id) {
                if (prevs[id] >= 0)
                    continue;
                flag[id] = 1;
                int32_t nx = nexts[id];
                int count = 1;
                while (nx >= 0) {
                    flag[nx] = (count % 2) ? 0 : 1;
                    nx = nexts[nx];
                    ++count;
                }
            }

            // Step 3: rebuild the box list and remap the fine cell index. A
            // discarded box is always claimed by the kept box preceding it, so
            // every fine cell maps to a valid new box.
            std::vector<Box> new_boxes;
            new_boxes.reserve(nb);
            std::vector<uint32_t> new_id(nb, 0xFFFFFFFFu);
            for (size_t id = 0; id < nb; ++id) {
                if (!flag[id])
                    continue;
                uint32_t out = (uint32_t) new_boxes.size();
                new_id[id] = out;
                if (nexts[id] >= 0) {
                    uint32_t nid = (uint32_t) nexts[id];
                    new_id[nid] = out;
                    new_boxes.push_back(union_box(boxes[id], boxes[nid]));
                } else {
                    new_boxes.push_back(boxes[id]);
                }
            }

            for (size_t c = 0; c < n; ++c)
                index[c] = new_id[index[c]];

            boxes.swap(new_boxes);
        };

        // Iterate X/Y/Z passes until the box count stops shrinking. The
        // empty_mask restricts merges to progressively coarser-aligned
        // boundaries in early iterations and becomes unrestricted afterwards.
        int iter = 0;
        size_t prev;
        do {
            prev = boxes.size();
            int32_t empty_mask = iter > 3 ? 0 : ((1 << (iter + 1)) - 1);
            merge_pass(0, empty_mask);
            merge_pass(1, empty_mask);
            merge_pass(2, empty_mask);
            ++iter;
        } while ((ScalarFloat) boxes.size() < m_alpha * (ScalarFloat) prev);

        // Flatten into the storage layout shared with extremum_irregular
        m_box_count = (uint32_t) boxes.size();
        std::vector<int32_t> lo(3 * boxes.size()), hi(3 * boxes.size());
        std::vector<ScalarFloat> extrema(2 * boxes.size());
        for (size_t i = 0; i < boxes.size(); ++i) {
            lo[i * 3 + 0] = boxes[i].lo[0];
            lo[i * 3 + 1] = boxes[i].lo[1];
            lo[i * 3 + 2] = boxes[i].lo[2];
            hi[i * 3 + 0] = boxes[i].hi[0];
            hi[i * 3 + 1] = boxes[i].hi[1];
            hi[i * 3 + 2] = boxes[i].hi[2];
            extrema[i * 2 + 0] = boxes[i].mn;
            extrema[i * 2 + 1] = boxes[i].mj;
        }

        m_quality_at_build = compute_quality(fine_data, index.data(),
                                             extrema.data(), n);

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

    ScalarFloat m_rebuild_threshold;
    ScalarFloat m_quality_at_build = 1.f;

    ScalarAffineTransform4f m_to_local;

    /// Merge cost weights (null-collision and traversal terms)
    ScalarFloat m_cn;
    ScalarFloat m_ct;
    /// Merge-loop stopping factor in (0, 1)
    ScalarFloat m_alpha;
};

MI_EXPORT_PLUGIN(ExtremumIrregularGrid2)
NAMESPACE_END(mitsuba)
