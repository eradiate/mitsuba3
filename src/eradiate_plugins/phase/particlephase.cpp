#include <algorithm>
#include <array>
#include <cmath>
#include <vector>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/string.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/eradiate/phase_utils.h>
#include <iomanip>
#include <sstream>


NAMESPACE_BEGIN(mitsuba)

/**!

.. _phase-particlephase:

Particle phase function (:monosp:`particlephase`)
-------------------------------------------------

.. pluginparameters::

 * - r_eff_volume, v_eff_volume
   - |volume|
   - Effective radius and effective variance of the particulate medium in the
     scene space.
   - |exposed|, |differentiable|

 * - r_eff_grid, v_eff_grid
   - |tensor|
   - Effective radius and effective variance increasing irregular grids. Shape
     :math:`(n_r, 1)`, :math:`(n_v, 1)`, respectively.
   - |exposed|

 * - grid_start
   - |tensor| (unsigned integer)
   - Sparse index segmenting ``nodes`` and ``phase_mueller`` into one sub-range
     per ``(r_eff_grid, v_eff_grid)`` entry. Element ``e`` is the start of
     entry ``e``'s sub-range; the last element is the number of
     ``nodes``/``phase_mueller`` elements in use. Shape
     :math:`(n_r \times n_v + 1, 1)`.
   - |exposed|

 * - nodes
   - |tensor|
   - Sparsely indexed scattering angle cosine discretizations, optionally
     padded past ``grid_start``'s last element. Shape :math:`(N, 1)`.
   - |exposed|

 * - phase_mueller
   - |tensor|
   - Sparsely indexed 6-components-per-node Mueller matrix table. Shape
     :math:`(N, 6)`.
   - |exposed|, |differentiable|, |discontinuous|

 * - sigma_s_weight
   - |tensor|
   - Scattering probability weight, one value per ``(r_eff_grid, v_eff_grid)``
     entry, used to rescale the weight of blended phase functions. Shape
     :math:`(n_r, n_v)`.
   - |exposed|, |differentiable|

This plugin implements a phase function tabulated on a 2D grid of effective
radius (``r_eff_grid``) and effective variance (``v_eff_grid``), each grid entry
carrying its own tabulated, polarized phase matrix over an independent irregular
scattering :math:`\cos\theta` discretization (``nodes``).

At a given interaction point, ``r_eff_volume``/``v_eff_volume`` are evaluated
to find the corresponding ``(r_eff, v_eff)`` couple. The 4 grid entries
bracketing that value are located in ``r_eff_grid``/``v_eff_grid``. A bilinear
blend weight is calculated for each corner, and further rescaled by the
corresponding scattering probability in ``sigma_s_weight``.

The irregular grids of :math:`\cos\theta` are all different in length and
values. This plugin uses an index defined by ``grid_start`` with one position
``e`` for each couple in ``(r_eff_grid, v_eff_grid)``. The contiguous portion of
``nodes`` and ``phase_mueller`` located in ``[grid_start[e]; grid_start[e+1][``
defines a separate phase function corresponding to that couple in
``(r_eff_grid, v_eff_grid)``. This configuration allows an efficient memory
representation of the ragged structure of angles typical of Mie phase
function lookup tables.

All phase discretizations encoded in ``nodes`` and indexed via ``grid_start``
must cover the full range :math:`[-1; 1]` of :math:`\cos\theta`.

.. note::
   ``r_eff_grid``, ``v_eff_grid``, ``nodes``, ``phase_mueller``,
   ``grid_start`` and ``sigma_s_weight`` can only be specified when building
   this plugin at runtime from Python or C++ and cannot be specified in the
   XML scene description.

.. note::
   Updating any of ``r_eff_grid``, ``v_eff_grid``, ``nodes``,
   ``phase_mueller``, ``grid_start`` or ``sigma_s_weight`` through the
   differentiable parameter map and calling ``params.update()`` rebuilds
   the internal state that depends on it, including CDF recalculation.
   ``nodes``/``phase_mueller`` can be over-allocated and padded so the
   active length ``grid_start`` indexes into them can vary between
   updates.
*/

template <typename Float, typename Spectrum>
class ParticlePhaseFunction final : public PhaseFunction<Float, Spectrum> {
public:
    MI_IMPORT_BASE(PhaseFunction, m_flags, m_components)
    MI_IMPORT_TYPES(PhaseFunctionContext, Volume)

    using FloatStorage  = DynamicBuffer<Float>;
    using UInt32Storage = DynamicBuffer<UInt32>;
    using Vector6u      = dr::Array<UInt32, 6>;
    using Vector6f      = dr::Array<Float, 6>;

    struct BlendWeight {
        Vector4u idx;
        Vector4f w;
    };


    /// Load property `name` as a `storage_t` tensor.
    template <typename storage_t>
    dr::Tensor<storage_t> load_grid(const Properties &props, const char *name) {
        return props.get_any<dr::Tensor<storage_t>>(name);
    }

    explicit ParticlePhaseFunction(const Properties &props) : Base(props) {
        m_r_eff_volume   = props.get_volume<Volume>("r_eff_volume");
        m_v_eff_volume   = props.get_volume<Volume>("v_eff_volume");
        m_r_eff_grid     = load_grid<FloatStorage>(props, "r_eff_grid");
        m_v_eff_grid     = load_grid<FloatStorage>(props, "v_eff_grid");
        m_nodes          = load_grid<FloatStorage>(props, "nodes");
        m_mueller        = load_grid<FloatStorage>(props, "phase_mueller");
        m_sigma_s_weight = load_grid<FloatStorage>(props, "sigma_s_weight");
        m_grid_start     = load_grid<UInt32Storage>(props, "grid_start");

        check_sizes();
        check_volume_range();
        precompute_cdf();

        m_flags = +PhaseFunctionFlags::Anisotropic;
        m_components.clear();
        m_components.push_back(m_flags);
    }

private:

    // =========================================================================
    //  CDF Calculation
    // =========================================================================

    /// Validate that a r_eff/v_eff grid axis is strictly increasing, so
    /// that it can be located with a binary search (\ref searchsorted).
    /// Validate the [-1.f, 1.f] range
    void validate_grid_axis(const FloatStorage &grid, size_t n, const char *name) {
        for (size_t i = 1; i < n; ++i) {
            ScalarFloat prev = dr::slice(grid, i - 1);
            ScalarFloat cur  = dr::slice(grid, i);
            if (cur <= prev)
                Throw("ParticlePhaseFunction: %s not strictly increasing: "
                      "%s[%zu]=%.17g >= %s[%zu]=%.17g.",
                      name, name, i - 1, prev, name, i, cur);
        }
    }

    /// Refresh m_n_r/m_n_v and validate every array size invariant.
    void check_sizes() {
        const FloatStorage  &r_eff_grid = m_r_eff_grid.array();
        const FloatStorage  &v_eff_grid = m_v_eff_grid.array();
        const UInt32Storage &grid_start = m_grid_start.array();

        m_n_r = dr::width(m_r_eff_grid);
        m_n_v = dr::width(m_v_eff_grid);
        validate_grid_axis(r_eff_grid, m_n_r, "r_eff_grid");
        validate_grid_axis(v_eff_grid, m_n_v, "v_eff_grid");

        size_t n = m_n_r * m_n_v;
        if (dr::width(m_mueller) != dr::width(m_nodes) * 6u)
            Throw("ParticlePhaseFunction: phase_mueller must have 6 * len(nodes) elements");
        if (dr::width(m_sigma_s_weight) != n)
            Throw("ParticlePhaseFunction: sigma_s_weight must have n_r * n_v = %zu "
                  "elements, got %zu.", n, dr::width(m_sigma_s_weight));
        if (dr::width(m_grid_start) != n + 1)
            Throw("ParticlePhaseFunction: grid_start must have n_r * n_v + 1 = %zu "
                  "elements (last one is the sentinel len(nodes)), got %zu.",
                  n + 1, dr::width(m_grid_start));
        if (dr::slice(grid_start, n) > dr::width(m_nodes))
            Throw("ParticlePhaseFunction: grid_start's last element must not exceed "
                  "len(nodes) = %zu, got %zu.",
                  dr::width(m_nodes), (size_t) dr::slice(grid_start, n));
    }

    /// Check r_eff_volume/v_eff_volume's value range against r_eff_grid/v_eff_grid's coverage.
    void check_volume_range() const {
        if (m_n_r > 1) {
            const FloatStorage &r_eff_grid = m_r_eff_grid.array();
            ScalarFloat r_min = dr::slice(r_eff_grid, 0);
            ScalarFloat r_max = dr::slice(r_eff_grid, m_n_r - 1);
            ScalarFloat vol_min = m_r_eff_volume->min();
            ScalarFloat vol_max = m_r_eff_volume->max();
            if (dr::isnan(vol_min) || dr::isnan(vol_max) ||
                vol_min < r_min || vol_max > r_max)
                Throw("ParticlePhaseFunction: r_eff_volume's value range "
                      "[%s, %s] is invalid or exceeds r_eff_grid's [%s, %s].",
                      vol_min, vol_max, r_min, r_max);
        }
        if (m_n_v > 1) {
            const FloatStorage &v_eff_grid = m_v_eff_grid.array();
            ScalarFloat v_min = dr::slice(v_eff_grid, 0);
            ScalarFloat v_max = dr::slice(v_eff_grid, m_n_v - 1);
            ScalarFloat vol_min = m_v_eff_volume->min();
            ScalarFloat vol_max = m_v_eff_volume->max();
            if (dr::isnan(vol_min) || dr::isnan(vol_max) ||
                vol_min < v_min || vol_max > v_max)
                Throw("ParticlePhaseFunction: v_eff_volume's value range "
                      "[%s, %s] is invalid or exceeds v_eff_grid's [%s, %s].",
                      vol_min, vol_max, v_min, v_max);
        }
    }

    /// Validate the node layout, then build the per-entry CDFs and their
    /// normalization factors (\ref m_cdf, \ref m_norm), vectorized over one
    /// lane per (r_eff, v_eff) entry.
    void precompute_cdf() {
        const FloatStorage  &nodes      = m_nodes.array();
        const UInt32Storage &grid_start = m_grid_start.array();

        size_t total_pts = dr::width(m_nodes);
        size_t n = m_n_r * m_n_v;

        size_t expected_start = 0;

        for (size_t i = 0; i < n; ++i) {
            auto [start, len] = entry_range_host(i);

            if (len < 2)
                Throw("ParticlePhaseFunction: entry %zu has grid_len=%zu (need >= 2).", i, len);
            if (start + len > total_pts)
                Throw("ParticlePhaseFunction: entry %zu slice [%zu, %zu) exceeds nodes size %zu.",
                      i, start, start + len, total_pts);
            if (start != expected_start)
                Throw("ParticlePhaseFunction: entry %zu grid_start=%zu, expected %zu "
                      "(start/len must form a contiguous concatenation; is 'nodes' a "
                      "global union instead of the per-entry concatenation?).",
                      i, start, expected_start);
            expected_start = start + len;

            bool invalid_s = dr::slice(nodes, start) != -1.f;
            bool invalid_e = dr::slice(nodes, start + len - 1) != 1.f;
            if (invalid_s || invalid_e)
                Throw(
                    "ParticlePhaseFunction: invalid phase bound for entry %zu. Start is 1.0: %b, end is 1.0: %b",
                    i, invalid_s, invalid_e
                );

            for (size_t k = 1; k < len; ++k) {
                ScalarFloat prev = dr::slice(nodes, start + k - 1);
                ScalarFloat cur  = dr::slice(nodes, start + k);
                if (cur <= prev)
                    Throw("ParticlePhaseFunction: entry %zu not strictly increasing: "
                          "nodes[%zu]=%.17g >= nodes[%zu]=%.17g (%s).",
                          i, start + k - 1, prev, start + k, cur,
                          cur == prev ? "duplicate" : "decreasing - descending cos(theta)?");
            }
        }

        if (expected_start > total_pts)
            Throw("ParticlePhaseFunction: entries cover %zu nodes, 'nodes' has %zu.",
                  expected_start, total_pts);

        size_t active_pts = expected_start;

        // Vectorized CDF construction via a segmented prefix sum, LLVM
        // friendly
        using MaskStorage = dr::mask_t<FloatStorage>;
        m_norm = dr::zeros<FloatStorage>(n);

        // Helper indexing arrays
        UInt32Storage entry = dr::arange<UInt32Storage>((uint32_t) n);
        UInt32Storage idx = dr::arange<UInt32Storage>((uint32_t) active_pts);
        UInt32Storage owning_entry = searchsorted<UInt32Storage>(
            grid_start, 0u, (uint32_t) (n + 1), idx + 1u, MaskStorage(true));
        UInt32Storage owning_start = dr::gather<UInt32Storage>(grid_start, owning_entry);
        MaskStorage is_first = idx == owning_start;
        UInt32Storage idx_prev = dr::select(is_first, idx, idx - 1u);

        auto [start, len] = entry_range(entry);

        // Fill increment with the trapezoid eval of each M11 component of the
        // index
        FloatStorage x1 = dr::gather<FloatStorage>(nodes, idx);
        FloatStorage x0 = dr::gather<FloatStorage>(nodes, idx_prev);
        FloatStorage y0 = gather_mueller_m11<FloatStorage>(idx_prev);
        FloatStorage y1 = gather_mueller_m11<FloatStorage>(idx);
        FloatStorage increment = dr::select(
            is_first, FloatStorage(0.f), 0.5f * (y0 + y1) * (x1 - x0));

        // Compute the global prefix sum of the increment and subtract the
        // cumulated value at the start of the entry, removing what leaked
        // in from the previous entry
        FloatStorage global_incl = dr::prefix_sum(increment, /*exclusive*/ false);
        FloatStorage global_incl_at_start = dr::gather<FloatStorage>(global_incl, owning_start);
        m_cdf = global_incl - global_incl_at_start;

        // Normalize
        UInt32Storage last_idx = start + len - 1u;
        FloatStorage running = dr::gather<FloatStorage>(m_cdf, last_idx);
        MaskStorage has_mass = running > 0.f;
        FloatStorage norm_per_entry = dr::select(has_mass, dr::rcp(running), FloatStorage(0.f));
        dr::scatter(m_norm, norm_per_entry, entry);
        FloatStorage norm_per_node = dr::gather<FloatStorage>(norm_per_entry, owning_entry);
        m_cdf *= norm_per_node;
    }

private:

    // =========================================================================
    //  Weights Calculation
    // =========================================================================

    /// Locate `value` along a `(r_eff, v_eff)` grid axis
    std::pair<UInt32, Float> locate_axis(const FloatStorage &grid, uint32_t n,
                                         Float value, Float value_min,
                                         Float value_max, Mask active) const {
        if (n <= 1)
            return { 0u, Float(0.f) };

        UInt32 i = searchsorted<UInt32>(grid, 0u, n, value, active);
        Vector2f v01 = dr::gather<Vector2f>(grid, Vector2u(i, i + 1u), active);
        Float v0 = v01.x(), v1 = v01.y();
        Float d  = v1 - v0;
        Float scale = dr::maximum(dr::abs(value_max - value_min), Float(1.f));
        Float t = dr::select(dr::abs(d) > dr::Epsilon<Float> * scale,
                             (value - v0) / d, Float(0.f));
        return { i, t };
    }

private:
    // =========================================================================
    //  Helper methods
    // =========================================================================

    /// Gather the 6 Mueller matrix components stored at `point`
    Vector6f gather_mueller(UInt32 point, Mask active) const {
        return dr::gather<Vector6f>(m_mueller.array(), point * 6u + Vector6u(0u, 1u, 2u, 3u, 4u, 5u), active);
    }

    /// Gather the m11 Mueller component stored at node index `idx`. `FloatX`
    /// and `UInt32X` must be the same width (e.g `Float`/`UInt32`, `Vector4f`
    /// /`Vector4u` or `FloatStorage`/`UInt32Storage`).
    template <typename FloatX, typename UInt32X>
    FloatX gather_mueller_m11(const UInt32X &idx,
                              dr::mask_t<UInt32X> active = true) const {
        return dr::gather<FloatX>(m_mueller.array(), idx * 6u, active);
    }

    /// Gather the pair of cosine values stored at node indices `idx_lo` and
    /// `idx_hi` (see \ref gather_mueller_m11).
    template <typename FloatX, typename UInt32X>
    std::pair<FloatX, FloatX> gather_node_pair(const UInt32X &idx_lo, const UInt32X &idx_hi,
                                                dr::mask_t<UInt32X> active = true) const {
        const FloatStorage &nodes = m_nodes.array();
        return { dr::gather<FloatX>(nodes, idx_lo, active),
                 dr::gather<FloatX>(nodes, idx_hi, active) };
    }

    /// Gather the pair of m11 Mueller components stored at node indices
    /// `idx_lo` and `idx_hi` (see \ref gather_mueller_m11).
    template <typename FloatX, typename UInt32X>
    std::pair<FloatX, FloatX> gather_mueller_m11_pair(const UInt32X &idx_lo, const UInt32X &idx_hi,
                                                       dr::mask_t<UInt32X> active = true) const {
        return { gather_mueller_m11<FloatX>(idx_lo, active),
                 gather_mueller_m11<FloatX>(idx_hi, active) };
    }

    /// Gather the (start, len) pair describing a grid entry's own sub-range
    /// of `m_nodes`/`m_mueller`, derived from `m_grid_start[entry]` and
    /// `m_grid_start[entry + 1]`. `entry` may be a scalar-width or
    /// vectorized (e.g. `Vector4u`)
    template <typename UInt32X>
    std::pair<UInt32X, UInt32X> entry_range(const UInt32X &entry,
                                            dr::mask_t<UInt32X> active = true) const {
        const UInt32Storage &grid_start = m_grid_start.array();
        UInt32X start = dr::gather<UInt32X>(grid_start, entry,      active);
        UInt32X next  = dr::gather<UInt32X>(grid_start, entry + 1u, active);
        return { start, next - start };
    }

    /// Host-side counterpart of entry_range(): the (start, len) pair for
    /// grid entry `i`, read via dr::slice(). Used where concrete values
    /// are needed right away (e.g. validation Throw() messages)
    std::pair<size_t, size_t> entry_range_host(size_t i) const {
        const UInt32Storage &grid_start = m_grid_start.array();
        size_t start = (size_t) dr::slice(grid_start, i);
        size_t next  = (size_t) dr::slice(grid_start, i + 1);
        return { start, next - start };
    }

    /// Locate sample `u` in a CDF entry and perform the quadratic
    /// inversion of its cosine.
    Float sample_cos_theta_entry(Float u, UInt32 start, UInt32 len,
                                 Float norm, Mask active) const {
        UInt32 lo = searchsorted<UInt32>(m_cdf, start, len, u, active);
        UInt32 hi = lo + 1u;

        auto [x0, x1] = gather_node_pair<Float>(lo, hi, active);
        auto [y0, y1] = gather_mueller_m11_pair<Float>(lo, hi, active);
        Float c0      = dr::gather<Float>(m_cdf, lo, active);
        y0 *= norm;
        y1 *= norm;

        Float h = x1 - x0;
        Float s = (u - c0) / h;

        // Numerically stable form of the quadratic root. The classical
        // form (y0 - sqrt(y0^2 + 2*s*(y1-y0))) / (y0 - y1) subtracts two
        // nearly-equal values whenever y0 ~ y1.
        Float disc = dr::fmadd(y0, y0, 2.f*s*(y1-y0));
        Float t    = 2.f*s * dr::rcp(y0 + dr::safe_sqrt(disc));
        return dr::fmadd(dr::clip(t, 0.f, 1.f), h, x0);
    }

private:
    // =========================================================================
    //  Search helper
    // =========================================================================

    /// Specialized call to dr::binary_search for readability purposes.
    template <typename Index, typename ArrT, typename BoundT, typename Value, typename MaskT>
    Index searchsorted(const ArrT &arr, BoundT start, BoundT len,
                        Value query, MaskT active) const {
        Index lo = dr::binary_search<Index>(
            start, start + len,
            [&arr, active, query](Index i) DRJIT_INLINE_LAMBDA {
                return active && (dr::gather<Value>(arr, i, active) < query);
            });
        Index result = dr::clip(lo, start + 1u, start + len - 1u) - 1u;
        return result;
    }

private:

    // =========================================================================
    //  Blend method implementation
    // =========================================================================

    /// Select one of the 4 corner's CDFs using \c u_select as a weight and
    /// invert it.
    std::tuple<Vector3f, Spectrum, Float>
    sample_stochastic(const PhaseFunctionContext &ctx,
                         const MediumInteraction3f &mei,
                         Float u_select, const Point2f &sample2,
                         const BlendWeight &bw, Mask active) const {
        Vector4f cumw;
        cumw[0] = bw.w[0];
        for (int c = 1; c < 4; ++c)
            cumw[c] = cumw[c - 1] + bw.w[c];

        // Corner selected by `u_select`, and its weight.
        UInt32 idx   = bw.idx[0];
        Float  w_sel = bw.w[0];
        for (int c = 1; c < 4; ++c) {
            Mask sel = u_select > cumw[c - 1];
            idx      = dr::select(sel, bw.idx[c], idx);
            w_sel    = dr::select(sel, bw.w[c], w_sel);
        }

        auto [s, l] = entry_range(idx, active);
        Float  norm = dr::gather<Float>(m_norm, idx, active);

        Float cos_t = sample_cos_theta_entry(sample2.x(), s, l, norm, active);
        Float sin_t = dr::safe_sqrt(1.f - cos_t * cos_t);
        auto [sin_phi, cos_phi] = dr::sincos(2.f * dr::Pi<ScalarFloat> * sample2.y());
        Vector3f wo = -mei.to_world(Vector3f(sin_t * cos_phi, sin_t * sin_phi, cos_t));
        wo = dr::select(active, wo, Vector3f(0.f));

        auto [phase_val, pdf] = eval_pdf_flat(ctx, mei, wo, bw, active);

        Spectrum weight = phase_val * dr::rcp(pdf);

        // Gradient correction for the discrete corner choice above;
        // does not impact the weight, but injects d(log w_sel)/d(theta)
        // into weight's gradient.
        if constexpr (dr::is_diff_v<Float>) {
            if (dr::grad_enabled(w_sel)) {
                Float score = dr::replace_grad(Float(1.f), w_sel / dr::detach(w_sel));
                weight *= score;
            }
        }

        weight = dr::select(active, weight, Spectrum(0.f));
        pdf    = dr::select(active, pdf, Float(0.f));

        return { wo, weight, pdf };
    }

    /// Locate the 4 entry indices in r_eff_grid and v_eff_grid bracketing the
    /// (r_eff, v_eff) value at the interaction point, and compute the
    /// associated 4 blend weights. The r_eff and v_eff grids need to be
    /// strictly increasing.
    BlendWeight get_blend_weight(const MediumInteraction3f &mei, Mask active) const {
        Float r_eff = m_r_eff_volume->eval_1(mei, active);
        Float v_eff = m_v_eff_volume->eval_1(mei, active);

        // These sanity checks can only be done in scalar variants.
        if constexpr (!dr::is_jit_v<Float>) {
            Mask nan_mask = dr::isnan(r_eff) || dr::isnan(v_eff);
            // Here we assess if the sampled r_eff or v_eff is defined.
            // If it isn't, this typically means the volume which possess
            // this phase has a non null extinction where the phase is 
            // not defined, or is incorrectly sampled there. This is the
            // sign of a severe problem in the scene construction or this
            // plugin's parameters.
            if (dr::any(active && nan_mask)) {
                std::ostringstream oss;
                oss << "ParticlePhaseFunction: NaN r_eff or v_eff at interaction point.\n";
                oss << "  mei.p       = " << mei.p << "\n";
                oss << "  r_eff value = " << r_eff << "\n";
                oss << "  v_eff value = " << v_eff << "\n";
                oss << "  r_eff is NaN: " << dr::isnan(r_eff) << "\n";
                oss << "  v_eff is NaN: " << dr::isnan(v_eff) << "\n";
                oss << "  r_eff_volume to_local:\n" << m_r_eff_volume->bbox() << "\n";
                oss << "  v_eff_volume to_local:\n" << m_v_eff_volume->bbox() << "\n";
                Throw("%s", oss.str());
            }
        }

        const FloatStorage &r_eff_grid = m_r_eff_grid.array();
        const FloatStorage &v_eff_grid = m_v_eff_grid.array();

        Float r_min = 0.f, r_max = 0.f;
        if (m_n_r > 1) {
            Vector2f r_bounds = dr::gather<Vector2f>(
                r_eff_grid, Vector2u(0u, (uint32_t) (m_n_r - 1)), active);
            r_min = r_bounds.x();
            r_max = r_bounds.y();
        }

        Float v_min = 0.f, v_max = 0.f;
        if (m_n_v > 1) {
            Vector2f v_bounds = dr::gather<Vector2f>(
                v_eff_grid, Vector2u(0u, (uint32_t) (m_n_v - 1)), active);
            v_min = v_bounds.x();
            v_max = v_bounds.y();
        }

        auto [ir, tr] = locate_axis(r_eff_grid, (uint32_t) m_n_r, r_eff, r_min, r_max, active);
        auto [iv, tv] = locate_axis(v_eff_grid, (uint32_t) m_n_v, v_eff, v_min, v_max, active);

        UInt32 n_v   = (uint32_t) m_n_v;
        UInt32 ir_hi = m_n_r > 1 ? ir + 1u : ir;
        UInt32 iv_hi = m_n_v > 1 ? iv + 1u : iv;

        // The 4 bracketing corners, in (r_eff, v_eff) order: (lo,lo), (lo,hi), (hi,lo), (hi,hi)
        Vector4u ir4(ir, ir, ir_hi, ir_hi);
        Vector4u iv4(iv, iv_hi, iv, iv_hi);
        Vector4f tr4(1.f - tr, 1.f - tr, tr, tr);
        Vector4f tv4(1.f - tv, tv, 1.f - tv, tv);

        BlendWeight bw;
        bw.idx = ir4 * n_v + iv4;
        bw.w   = tr4 * tv4;

        // Rescale each corner according to its scattering probability
        Vector4f sca  = dr::gather<Vector4f>(m_sigma_s_weight.array(), bw.idx, active);
        bw.w         *= sca;
        Float sca_sum = dr::sum(bw.w);
        bw.w          = dr::select(sca_sum > dr::Epsilon<Float>, bw.w / sca_sum, bw.w);

        return bw;
    }

private:

    // =========================================================================
    //  Post-processing
    // =========================================================================

    /// Locate a cosine in the 4 corner's PDF sub-sections, blend, and
    /// evaluate the PDF given the associated blend weights.
    std::pair<Spectrum, Float>
    eval_pdf_flat(const PhaseFunctionContext &ctx,
                  const MediumInteraction3f &mei,
                  const Vector3f &wo,
                  const BlendWeight &bw,
                  Mask active) const {
        Float cos_t     = -dr::dot(wo, mei.wi);
        const FloatStorage &nodes = m_nodes.array();
        auto [starts, lens] = entry_range(bw.idx, active);
        Vector4u lo4;
        for (int c = 0; c < 4; ++c)
            lo4[c] = searchsorted<UInt32>(nodes, starts[c], lens[c], cos_t, active);

        Vector4u hi4 = lo4 + 1u;
        auto [x0_4, x1_4] = gather_node_pair<Vector4f>(lo4, hi4, active);
        auto [domain_lo, domain_hi] =
            gather_node_pair<Vector4f>(starts, starts + lens - 1u, active);

        Vector4f t4    = (Vector4f(cos_t) - x0_4) / (x1_4 - x0_4);
        Vector4f norms = dr::gather<Vector4f>(m_norm, bw.idx, active);
        Vector4f w     = dr::select(Vector4f(cos_t) >= domain_lo && Vector4f(cos_t) <= domain_hi,
                                    bw.w, Vector4f(0.f));

        Spectrum phase_val(0.f);
        Float pdf;

        if constexpr (is_polarized_v<Spectrum>) {
            Vector4f m11_4;
            for (size_t c = 0; c < 4; ++c) {
                Vector6f v0  = gather_mueller(lo4[c], active);
                Vector6f v1  = gather_mueller(hi4[c], active);
                Vector6f mi  = dr::fmadd(t4[c], v1 - v0, v0);
                m11_4[c]     = mi[0];
                Float scale  = norms[c] * dr::InvTwoPi<ScalarFloat> * w[c];
                phase_val += MuellerMatrix<Float>(
                    mi[0]*scale,  mi[1]*scale, 0.f,            0.f,
                    mi[1]*scale,  mi[2]*scale, 0.f,            0.f,
                    0.f,            0.f,           mi[3]*scale,  mi[4]*scale,
                    0.f,            0.f,          -mi[4]*scale,  mi[5]*scale);
            }
            pdf = dr::dot(w * norms, m11_4) * dr::InvTwoPi<ScalarFloat>;

            Vector3f wo_hat = ctx.mode == TransportMode::Radiance ? wo : mei.wi,
                     wi_hat = ctx.mode == TransportMode::Radiance ? mei.wi : wo;
            Vector3f x_hat      = dr::cross(-wo_hat, wi_hat),
                     p_axis_in  = dr::normalize(dr::cross(x_hat, -wo_hat)),
                     p_axis_out = dr::normalize(dr::cross(x_hat,  wi_hat));
            phase_val = mueller::rotate_mueller_basis(
                phase_val,
                -wo_hat, p_axis_in,  mueller::stokes_basis(-wo_hat),
                 wi_hat, p_axis_out, mueller::stokes_basis( wi_hat));
            dr::masked(phase_val, dr::isnan(phase_val)) = depolarizer<Spectrum>(0.f);
        } else {
            auto [y0_m11, y1_m11] = gather_mueller_m11_pair<Vector4f>(lo4, hi4, active);
            Vector4f m11_4        = dr::fmadd(t4, y1_m11 - y0_m11, y0_m11);
            pdf = dr::dot(w * norms, m11_4) * dr::InvTwoPi<ScalarFloat>;
            phase_val = Spectrum(pdf);
        }

        return { phase_val, pdf };
    }


public:

    std::tuple<Vector3f, Spectrum, Float>
    sample(const PhaseFunctionContext &ctx,
           const MediumInteraction3f &mei,
           Float sample1,
           const Point2f &sample2,
           Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::PhaseFunctionSample, active);

        BlendWeight bw = get_blend_weight(mei, active);

        return sample_stochastic(ctx, mei, sample1, sample2, bw, active);
    }

    std::pair<Spectrum, Float>
    eval_pdf(const PhaseFunctionContext &ctx,
             const MediumInteraction3f &mei,
             const Vector3f &wo,
             Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::PhaseFunctionEvaluate, active);

        BlendWeight bw = get_blend_weight(mei, active);
        auto [val, pdf] = eval_pdf_flat(ctx, mei, wo, bw, active);

        val = dr::select(active, val, Spectrum(0.f));
        pdf = dr::select(active, pdf, Float(0.f));

        return { val, pdf };
    }

    void parameters_changed(const std::vector<std::string> &keys = {}) override {
        bool volume_changed  = keys.empty() || string::contains(keys, "r_eff_volume") ||
                                string::contains(keys, "v_eff_volume");
        bool grid_resized    = keys.empty() || string::contains(keys, "r_eff_grid") ||
                                string::contains(keys, "v_eff_grid");
        bool nodes_changed   = keys.empty() || string::contains(keys, "nodes");
        bool mueller_changed = keys.empty() || string::contains(keys, "phase_mueller");
        bool start_changed   = keys.empty() || string::contains(keys, "grid_start");
        bool sigma_changed   = keys.empty() || string::contains(keys, "sigma_s_weight");

        if (grid_resized || nodes_changed || mueller_changed || sigma_changed || start_changed)
            check_sizes();
        if (volume_changed || grid_resized)
            check_volume_range();
        if (grid_resized || nodes_changed || mueller_changed || start_changed)
            precompute_cdf();
    }

    void traverse(TraversalCallback *cb) override {
        cb->put("r_eff_volume",  m_r_eff_volume.get(), ParamFlags::Differentiable);
        cb->put("v_eff_volume",  m_v_eff_volume.get(), ParamFlags::Differentiable);
        cb->put("r_eff_grid",    m_r_eff_grid,         ParamFlags::NonDifferentiable);
        cb->put("v_eff_grid",    m_v_eff_grid,         ParamFlags::NonDifferentiable);
        cb->put("nodes",           m_nodes,              ParamFlags::NonDifferentiable);
        cb->put("phase_mueller",   m_mueller,            ParamFlags::Differentiable | ParamFlags::Discontinuous);
        cb->put("grid_start",      m_grid_start,         ParamFlags::NonDifferentiable);
        cb->put("sigma_s_weight",  m_sigma_s_weight,     ParamFlags::Differentiable);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "ParticlePhaseFunction[" << std::endl
            << "  n_r = " << m_n_r << "," << std::endl
            << "  n_v = " << m_n_v << "," << std::endl
            << "  total_pts = " << dr::slice(m_grid_start.array(), m_n_r * m_n_v) << std::endl
            << "]";
        return oss.str();
    }

    FloatStorage get_envelope_nodes() const override {
        const UInt32Storage &grid_start = m_grid_start.array();
        const FloatStorage  &nodes      = m_nodes.array();
        size_t n = (size_t) dr::slice(grid_start, m_n_r * m_n_v);
        std::vector<ScalarFloat> tmp(n);
        for (size_t i = 0; i < n; ++i)
            tmp[i] = dr::slice(nodes, i);

        std::sort(tmp.begin(), tmp.end());
        auto last = std::unique(tmp.begin(), tmp.end());
        tmp.erase(last, tmp.end());

        return dr::load<FloatStorage>(tmp.data(), tmp.size());
    }


    void accumulate_envelope(const FloatStorage &nodes, FloatStorage &values) const override {
        const FloatStorage &own_nodes = m_nodes.array();
        const FloatStorage &mueller   = m_mueller.array();

        size_t n = dr::width(nodes);
        std::vector<ScalarFloat> tmp_values(n, ScalarFloat(0));

        for (size_t ir = 0; ir < m_n_r; ++ir) {
            for (size_t iv = 0; iv < m_n_v; ++iv) {
                size_t entry = ir * m_n_v + iv;
                auto [start, len] = entry_range_host(entry);
                ScalarFloat norm = dr::slice(m_norm,      entry);

                size_t lo = start;

                for (size_t i = 0; i < n; ++i) {
                    ScalarFloat cos_t = dr::slice(nodes, i);

                    while (lo < start + len - 2 && dr::slice(own_nodes, lo + 1) < cos_t)
                        ++lo;

                    size_t hi     = lo + 1;
                    ScalarFloat x0  = dr::slice(own_nodes, lo);
                    ScalarFloat x1  = dr::slice(own_nodes, hi);
                    ScalarFloat t   = dr::clip((cos_t - x0) / (x1 - x0), ScalarFloat(0.f), ScalarFloat(1.f));
                    ScalarFloat y0  = dr::slice(mueller, lo * 6u);
                    ScalarFloat y1  = dr::slice(mueller, hi * 6u);
                    ScalarFloat m11 = dr::fmadd(t, y1 - y0, y0) * norm * dr::InvTwoPi<ScalarFloat>;
                    tmp_values[i] = dr::maximum(tmp_values[i], m11);
                }
            }
        }

        values = dr::maximum(
	        dr::load<FloatStorage>(tmp_values.data(), tmp_values.size()),
	        values
        );
    }

    MI_DECLARE_CLASS(ParticlePhaseFunction)

private:
    ref<Volume> m_r_eff_volume;
    ref<Volume> m_v_eff_volume;

    dr::Tensor<FloatStorage>  m_r_eff_grid;
    dr::Tensor<FloatStorage>  m_v_eff_grid;
    dr::Tensor<FloatStorage>  m_nodes;
    dr::Tensor<FloatStorage>  m_mueller;
    dr::Tensor<UInt32Storage> m_grid_start;

    FloatStorage  m_cdf;
    FloatStorage  m_norm;
    dr::Tensor<FloatStorage>  m_sigma_s_weight;

    size_t m_n_r = 0;
    size_t m_n_v = 0;
};

MI_EXPORT_PLUGIN(ParticlePhaseFunction)

NAMESPACE_END(mitsuba)
