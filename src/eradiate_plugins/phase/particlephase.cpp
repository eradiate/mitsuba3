#include <algorithm>
#include <array>
#include <cmath>
#include <vector>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/volumegrid.h>
#include <iomanip>
#include <sstream>


NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
class ParticlePhaseFunction final : public PhaseFunction<Float, Spectrum> {
public:
    MI_IMPORT_BASE(PhaseFunction, m_flags, m_components)
    MI_IMPORT_TYPES(PhaseFunctionContext, Volume, VolumeGrid)

    using FloatStorage  = DynamicBuffer<Float>;
    using UInt32Storage = DynamicBuffer<UInt32>;
    using Vector6u      = dr::Array<UInt32, 6>;
    using Vector6f      = dr::Array<Float, 6>;

    enum class BlendingMethod { Stochastic, Reconstruct, Tabulate, Bisect, Search };

    static BlendingMethod method_from_string(const std::string &method) {
        if (method == "stochastic")  return BlendingMethod::Stochastic;
        if (method == "reconstruct") return BlendingMethod::Reconstruct;
        if (method == "tabulate")    return BlendingMethod::Tabulate;
        if (method == "bisect")      return BlendingMethod::Bisect;
        if (method == "search")      return BlendingMethod::Search;
        Throw("ParticlePhaseFunction: unknown blending_method '%s'", method);
    }

    static std::string method_to_string(BlendingMethod method) {
        switch (method) {
            case BlendingMethod::Stochastic:  return "stochastic";
            case BlendingMethod::Reconstruct: return "reconstruct";
            case BlendingMethod::Tabulate:     return "tabulate";
            case BlendingMethod::Bisect:       return "bisect";
            case BlendingMethod::Search:       return "search";
            default: Throw("ParticlePhaseFunction: unknown blending method");
        }
    }

    struct BilinearWeights {
        Vector4u idx;
        Vector4f w;
    };

    struct BisectionBounds {
        Vector4u starts, lens;
        Vector4f norm;
        Float lo, hi, target, tol;
        Vector4u lo4_at_hi;
    };

    template<typename data_t, typename storage_t>
    storage_t load_grid(const Properties &props, const char *name) {
        try {
            ref<Object> obj = props.get<ref<Object>>(name);
            auto *vg = dynamic_cast<VolumeGrid *>(obj.get());
            if (!vg)
                Throw("ParticlePhaseFunction: property \"%s\" must be a VolumeGrid.", name);
            size_t n = (size_t) vg->size().x();
            const ScalarFloat *buf = vg->data();
            storage_t out = dr::zeros<storage_t>(n);
            for (size_t i = 0; i < n; ++i)
                dr::scatter(out, data_t(buf[i]), UInt32(i));
            return out;
        } catch (...) {}
        return props.get_any<storage_t>(name);
    }

    explicit ParticlePhaseFunction(const Properties &props) : Base(props) {
        m_r_eff_volume = props.get_volume<Volume>("r_eff_volume");
        m_v_eff_volume = props.get_volume<Volume>("v_eff_volume");

        m_n_r = props.get<int>("n_r");
        m_n_v = props.get<int>("n_v");

        m_r_eff_grid    = load_grid<Float, FloatStorage>(props, "r_eff_grid");
        m_v_eff_grid    = load_grid<Float, FloatStorage>(props, "v_eff_grid");
        m_nodes         = props.get_any<FloatStorage>("nodes");
        m_mueller       = props.get_any<FloatStorage>("phase_mueller");
        m_grid_start    = load_grid<UInt32, UInt32Storage>(props, "grid_start");
        m_grid_len      = load_grid<UInt32, UInt32Storage>(props, "grid_len");
        m_sigma_s_weight = load_grid<Float, FloatStorage>(props, "sigma_s_weight");

        m_blending_method = method_from_string(props.get<std::string>("blending_method", "bisect"));

        if (dr::width(m_mueller) != dr::width(m_nodes) * 6u)
            Throw("ParticlePhaseFunction: phase_mueller must have 6 * len(nodes) elements");
        if (dr::width(m_r_eff_grid) != (size_t) m_n_r)
            Throw("ParticlePhaseFunction: r_eff_grid must have n_r elements");
        if (dr::width(m_v_eff_grid) != (size_t) m_n_v)
            Throw("ParticlePhaseFunction: v_eff_grid must have n_v elements");
        if (dr::width(m_sigma_s_weight) < (size_t)(m_n_r * m_n_v))
            Throw("ParticlePhaseFunction: sigma_s_weight must have at least n_r * n_v = %d elements", m_n_r * m_n_v);

        if (props.has_property("cdf") && props.has_property("norm")) {
            m_cdf  = props.get_any<FloatStorage>("cdf");
            m_norm = props.get_any<FloatStorage>("norm");
        } else if (props.has_property("cdf") || props.has_property("norm")) {
            Throw("ParticlePhaseFunction: cdf and norm must be provided together");
        } else {
            precompute_cdf();
        }

        m_flags = +PhaseFunctionFlags::Anisotropic;
        m_components.clear();
        m_components.push_back(m_flags);
    }

    void precompute_cdf() {
        size_t total_pts = dr::width(m_nodes);
        size_t n = (size_t)(m_n_r * m_n_v);

        m_cdf  = dr::zeros<FloatStorage>(total_pts);
        m_norm = dr::zeros<FloatStorage>(n);

        for (size_t i = 0; i < n; ++i) {
            size_t start = dr::slice(m_grid_start, i);
            size_t len   = dr::slice(m_grid_len,   i);

            if (len < 2)
                Throw("ParticlePhaseFunction: entry %zu has fewer than 2 nodes!", i);

            ScalarFloat running = 0.0;
            dr::scatter(m_cdf, Float(0.f), UInt32(start));
            for (size_t k = 0; k < len - 1; ++k) {
                ScalarFloat x0 = dr::slice(m_nodes,   start + k);
                ScalarFloat x1 = dr::slice(m_nodes,   start + k + 1);
                ScalarFloat y0 = dr::slice(m_mueller, (start + k)     * 6);
                ScalarFloat y1 = dr::slice(m_mueller, (start + k + 1) * 6);
                running += 0.5 * (y0 + y1) * (x1 - x0);
                dr::scatter(m_cdf, Float(running), UInt32(start + k + 1));
            }

            ScalarFloat total = running;
            if (total > 0.0) {
                dr::scatter(m_norm, Float(1.0 / total), UInt32(i));
                for (size_t k = 0; k < len; ++k) {
                    ScalarFloat val = dr::slice(m_cdf, start + k);
                    dr::scatter(m_cdf, Float(val / total), UInt32(start + k));
                }
            }
        }

        validate_node_layout();
    }

    void validate_node_layout() const {
        size_t total_pts = dr::width(m_nodes);
        size_t n_entries = (size_t)(m_n_r * m_n_v);

        size_t expected_start = 0;
        for (size_t e = 0; e < n_entries; ++e) {
            size_t start = dr::slice(m_grid_start, e);
            size_t len   = dr::slice(m_grid_len,   e);

            if (len < 2)
                Throw("ParticlePhaseFunction: entry %zu has grid_len=%zu (need >= 2).", e, len);
            if (start + len > total_pts)
                Throw("ParticlePhaseFunction: entry %zu slice [%zu, %zu) exceeds nodes size %zu.",
                      e, start, start + len, total_pts);
            if (start != expected_start)
                Throw("ParticlePhaseFunction: entry %zu grid_start=%zu, expected %zu "
                      "(start/len must form a contiguous concatenation; is 'nodes' a "
                      "global union instead of the per-entry concatenation?).",
                      e, start, expected_start);
            expected_start = start + len;

            ScalarFloat prev = dr::slice(m_nodes, start);
            for (size_t k = 1; k < len; ++k) {
                ScalarFloat cur = dr::slice(m_nodes, start + k);
                if (cur <= prev)
                    Throw("ParticlePhaseFunction: entry %zu not strictly increasing: "
                          "nodes[%zu]=%.17g >= nodes[%zu]=%.17g (%s).",
                          e, start + k - 1, prev, start + k, cur,
                          cur == prev ? "duplicate" : "decreasing - descending cos(theta)?");
                prev = cur;
            }
        }
        if (expected_start != total_pts)
            Throw("ParticlePhaseFunction: entries cover %zu nodes, 'nodes' has %zu.",
                  expected_start, total_pts);
    }

    BilinearWeights get_bilinear_weights(Float r_eff, Float v_eff,
                                         Float r_min, Float r_max,
                                         Float v_min, Float v_max,
                                         Mask active) const {
        UInt32 ir = 0u;
        Float tr = Float(0.f);
        if (m_n_r > 1) {
            Float dr_   = (r_max - r_min) / Float(m_n_r - 1);
            ir  = dr::clip(UInt32((r_eff - r_min) / dr_), 0u, UInt32(m_n_r - 2));
            Float r0 = dr::fmadd(Float(ir),      dr_, r_min);
            Float r1 = dr::fmadd(Float(ir + 1u), dr_, r_min);
            Float dri = r1 - r0;
            tr = dr::select(dr::abs(dri) > dr::Epsilon<Float>, (r_eff - r0) / dri, Float(0.f));
        }

        UInt32 iv = 0u;
        Float tv = Float(0.f);
        if (m_n_v > 1) {
            Float dv_   = (v_max - v_min) / Float(m_n_v - 1);
            iv  = dr::clip(UInt32((v_eff - v_min) / dv_), 0u, UInt32(m_n_v - 2));
            Float v0 = dr::fmadd(Float(iv),      dv_, v_min);
            Float v1 = dr::fmadd(Float(iv + 1u), dv_, v_min);
            Float dvi = v1 - v0;
            tv = dr::select(dr::abs(dvi) > dr::Epsilon<Float>, (v_eff - v0) / dvi, Float(0.f));
        }

        Float _tr = 1.f - tr, _tv = 1.f - tv;
        UInt32 n_v   = m_n_v;
        UInt32 ir_hi = m_n_r > 1 ? ir + 1u : ir;
        UInt32 iv_hi = m_n_v > 1 ? iv + 1u : iv;

        BilinearWeights bw;
        bw.idx = Vector4u(ir    * n_v + iv,
                          ir    * n_v + iv_hi,
                          ir_hi * n_v + iv,
                          ir_hi * n_v + iv_hi);
        bw.w   = Vector4f(_tr * _tv, _tr * tv, tr * _tv, tr * tv);

        Vector4f sca  = dr::gather<Vector4f>(m_sigma_s_weight, bw.idx, active);
        bw.w         *= sca;
        Float sca_sum = dr::sum(bw.w);
        bw.w          = dr::select(sca_sum > dr::Epsilon<Float>, bw.w / sca_sum, bw.w);

        return bw;
    }

    static constexpr Vector6u k_mueller_offsets() { return Vector6u(0u, 1u, 2u, 3u, 4u, 5u); }

    Vector6f gather_mueller(UInt32 point, Mask active) const {
        return dr::gather<Vector6f>(m_mueller, point * 6u + k_mueller_offsets(), active);
    }

    Vector4u ragged_searchsorted4(const FloatStorage &arr,
                                   Vector4u start, Vector4u len,
                                   Float query, Mask active) const {
        Vector4u lo = start;
        Vector4u hi = start + len;
        dr::mask_t<Vector4u> running = dr::mask_t<Vector4u>(active) && (lo < hi);

        Vector4f val;
        Vector4u mid;
        while (dr::any_nested(running)) {
            mid  = (lo + hi) >> 1u;
            val  = dr::gather<Vector4f>(arr, mid, running);
            dr::mask_t<Vector4u> go_right = running && (val < Vector4f(query));
            lo = dr::select(go_right, mid + 1u, lo);
            hi = dr::select(running && !go_right, mid, hi);
            running = running && (lo < hi);
        }

        return dr::clip(lo, start + 1u, start + len - 1u) - 1u;
    }

    UInt32 searchsorted1(const FloatStorage &arr, UInt32 start, UInt32 len,
                         Float query, Mask active) const {
        UInt32 lo = dr::binary_search<UInt32>(
            start, start + len,
            [&](UInt32 i) DRJIT_INLINE_LAMBDA {
                return active && (dr::gather<Float>(arr, i, active) < query);
            });
        UInt32 result = dr::clip(lo, start + 1u, start + len - 1u) - 1u;
        return result;
    }

    Float sample_cos_theta_entry(Float u, UInt32 start, UInt32 len,
                                 Float norm, Mask active) const {
        UInt32 lo = searchsorted1(m_cdf, start, len, u, active);
        UInt32 hi = lo + 1u;

        Float x0 = dr::gather<Float>(m_nodes,   lo,      active);
        Float x1 = dr::gather<Float>(m_nodes,   hi,      active);
        Float c0 = dr::gather<Float>(m_cdf,     lo,      active);
        Float y0 = dr::gather<Float>(m_mueller, lo * 6u, active) * norm;
        Float y1 = dr::gather<Float>(m_mueller, hi * 6u, active) * norm;

        Float h = x1 - x0;
        Float s = (u - c0) / h;
        Float t_linear = (y0 - dr::safe_sqrt(dr::fmadd(y0, y0, 2.f*s*(y1-y0)))) * dr::rcp(y0 - y1);
        Float t_const  = s * dr::rcp(y0);
        Float t        = dr::select(dr::abs(y1 - y0) < dr::Epsilon<Float>, t_const, t_linear);
        return dr::fmadd(dr::clip(t, 0.f, 1.f), h, x0);
    }


    Float eval_blended_cdf_at(Float x, const BilinearWeights &bw, Vector4u lo4,
                               Vector4f norm, Mask active) const {
        Vector4u hi4   = lo4 + 1u;
        Vector4f x0_4  = dr::gather<Vector4f>(m_nodes,   lo4,      active);
        Vector4f x1_4  = dr::gather<Vector4f>(m_nodes,   hi4,      active);
        Vector4f y0_4  = dr::gather<Vector4f>(m_mueller, lo4 * 6u, active);
        Vector4f y1_4  = dr::gather<Vector4f>(m_mueller, hi4 * 6u, active);
        Vector4f cdf0  = dr::gather<Vector4f>(m_cdf,     lo4,      active);
        Vector4f dx    = x1_4 - x0_4;
        Vector4f t4    = dr::clip(dr::select(dx > dr::Epsilon<Float>,
                                             (Vector4f(x) - x0_4) / dx, 0.f), 0.f, 1.f);
        Vector4f cdf_c = cdf0 + norm * dx * t4 * dr::fmadd(0.5f * t4, y1_4 - y0_4, y0_4);
        return dr::dot(bw.w, cdf_c);
    }

    Float eval_blended_cdf(Float x, const BilinearWeights &bw,
                            Vector4u starts, Vector4u lens, Vector4f norm, Mask active) const {
        Vector4u lo4 = ragged_searchsorted4(m_nodes, starts, lens, x, active);
        return eval_blended_cdf_at(x, bw, lo4, norm, active);
    }

    BisectionBounds
    bisection_bounds(Float u, const BilinearWeights &bw, Mask active) const {
        Vector4u starts = dr::gather<Vector4u>(m_grid_start, bw.idx, active);
        Vector4u lens   = dr::gather<Vector4u>(m_grid_len,   bw.idx, active);
        Vector4f norm   = dr::gather<Vector4f>(m_norm,       bw.idx, active);

        Float lo = dr::max(dr::gather<Vector4f>(m_nodes, starts,             active));
        Float hi = dr::min(dr::gather<Vector4f>(m_nodes, starts + lens - 1u, active));

        Vector4u lo4_at_hi = ragged_searchsorted4(m_nodes, starts, lens, hi, active);
        Float target = u * eval_blended_cdf_at(hi, bw, lo4_at_hi, norm, active);
        Float tol    = 16 * dr::Epsilon<Float> * dr::maximum(dr::abs(hi), Float(1.f));

        return { starts, lens, norm, lo, hi, target, tol, lo4_at_hi };
    }

    Float sample_cos_theta_bisect(Float u, const BilinearWeights &bw, Mask active) const {
        BisectionBounds bb = bisection_bounds(u, bw, active);
        Float lo = bb.lo, hi = bb.hi;

        constexpr size_t kMaxIters = 60;
        for (size_t it = 0; it < kMaxIters; ++it) {
            Mask refine = active && ((hi - lo) > bb.tol);
            if (!dr::any_nested(refine))
                break;
            Float mid = 0.5f * (lo + hi);
            Mask go_right = refine && (eval_blended_cdf(mid, bw, bb.starts, bb.lens, bb.norm, active) < bb.target);
            lo = dr::select(go_right,            mid, lo);
            hi = dr::select(refine && !go_right, mid, hi);
        }

        return 0.5f * (lo + hi);
    }

    Float sample_cos_theta_tabulate(Float u, const BilinearWeights &bw, Mask active) const {
        ScalarFloat lo_bound = 0.0, hi_bound = 0.0;
        auto entries = gather_column_entries<1>(bw, active, lo_bound, hi_bound);

        std::vector<ScalarFloat> merged;
        std::array<std::vector<ScalarFloat>, 1> Y_cols;
        merge_and_interpolate<1>(entries, true, lo_bound, hi_bound, merged, Y_cols);

        const std::vector<ScalarFloat> &Y = Y_cols[0];
        size_t n = merged.size();
        std::vector<ScalarFloat> cdf(n, 0.0);
        for (size_t k = 1; k < n; ++k)
            cdf[k] = cdf[k - 1] + 0.5 * (Y[k - 1] + Y[k]) * (merged[k] - merged[k - 1]);

        ScalarFloat total = cdf.back();
        if (n < 2 || total <= 0.0)
            return Float(lo_bound);

        ScalarFloat target = dr::slice(u, 0) * total;

        size_t hi_idx = (size_t)(std::lower_bound(cdf.begin(), cdf.end(), target) - cdf.begin());
        hi_idx = std::min(std::max(hi_idx, (size_t) 1), n - 1);
        size_t lo_idx = hi_idx - 1;

        ScalarFloat x0 = merged[lo_idx], x1 = merged[hi_idx];
        ScalarFloat y0 = Y[lo_idx],      y1 = Y[hi_idx];
        ScalarFloat c0 = cdf[lo_idx];

        ScalarFloat h = x1 - x0;
        ScalarFloat s = h > 0.0 ? (target - c0) / h : 0.0;

        ScalarFloat t;
        if (std::abs(y1 - y0) < 1e-12 * std::max({ std::abs(y0), std::abs(y1), ScalarFloat(1.0) }))
            t = y0 > 0.0 ? s / y0 : ScalarFloat(0.0);
        else
            t = (y0 - std::sqrt(std::max(ScalarFloat(0.0), y0 * y0 + ScalarFloat(2.0) * s * (y1 - y0)))) / (y0 - y1);
        t = std::min(std::max(t, ScalarFloat(0.0)), ScalarFloat(1.0));

        ScalarFloat cos_t = x0 + t * h;
        return Float(cos_t);
    }

    std::tuple<Float, Vector4u, Vector4u, Vector4u, Mask>
    sample_cos_theta_search(Float u, const BilinearWeights &bw, Mask active) const {
        auto [starts, lens, norm, lo, hi, target, tol, lo4_at_hi] = bisection_bounds(u, bw, active);

        Vector4u lo4_at_lo = ragged_searchsorted4(m_nodes, starts, lens, lo, active);

        constexpr size_t kMaxIters = 30;
        for (size_t it = 0; it < kMaxIters; ++it) {
            Mask stable = dr::all(lo4_at_lo == lo4_at_hi);
            Mask refine = active && !stable && ((hi - lo) > tol);
            if (!dr::any_nested(refine))
                break;

            Float mid = 0.5f * (lo + hi);
            Vector4u lo4_at_mid = ragged_searchsorted4(m_nodes, starts, lens, mid, active);
            Mask go_right = refine && (eval_blended_cdf_at(mid, bw, lo4_at_mid, norm, active) < target);

            Mask update_lo = refine && go_right;
            Mask update_hi = refine && !go_right;
            lo         = dr::select(update_lo, mid, lo);
            hi         = dr::select(update_hi, mid, hi);
            lo4_at_lo  = dr::select(update_lo, lo4_at_mid, lo4_at_lo);
            lo4_at_hi  = dr::select(update_hi, lo4_at_mid, lo4_at_hi);
        }

        Mask stable = active && dr::all(lo4_at_lo == lo4_at_hi);

        Vector4u lo4 = lo4_at_lo;
        Vector4u hi4 = lo4 + 1u;
        Vector4f x0_4   = dr::gather<Vector4f>(m_nodes,   lo4,      active);
        Vector4f x1_4   = dr::gather<Vector4f>(m_nodes,   hi4,      active);
        Vector4f y0_4   = dr::gather<Vector4f>(m_mueller, lo4 * 6u, active);
        Vector4f y1_4   = dr::gather<Vector4f>(m_mueller, hi4 * 6u, active);
        Vector4f cdf0_4 = dr::gather<Vector4f>(m_cdf,     lo4,      active);
        Vector4f dx4    = x1_4 - x0_4;

        Vector4f wn = bw.w * norm;
        Vector4f k2 = dr::select(dx4 > dr::Epsilon<Float>, 0.5f * wn * (y1_4 - y0_4) / dx4, Vector4f(0.f));
        Vector4f k1 = wn * y0_4;

        Float A = dr::sum(k2);
        Float B = dr::sum(k1 - 2.f * k2 * x0_4);
        Float C = dr::sum(bw.w * cdf0_4 - k1 * x0_4 + k2 * x0_4 * x0_4) - target;

        Float disc  = dr::maximum(dr::fmadd(B, B, -4.f * A * C), 0.f);
        Float sq    = dr::sqrt(disc);
        Float root1 = dr::select(dr::abs(A) > dr::Epsilon<Float>, (-B + sq) / (2.f * A), 0.f);
        Float root2 = dr::select(dr::abs(A) > dr::Epsilon<Float>, (-B - sq) / (2.f * A), 0.f);
        Float root_lin = dr::select(dr::abs(B) > dr::Epsilon<Float>, -C / B, 0.5f * (lo + hi));

        Mask root1_in = root1 >= lo && root1 <= hi;
        Float x_quad  = dr::select(root1_in, root1, root2);
        Float cos_t   = dr::select(dr::abs(A) > dr::Epsilon<Float>, x_quad, root_lin);
        cos_t = dr::select(stable, dr::clip(cos_t, lo, hi), 0.5f * (lo + hi));

        return { cos_t, starts, lens, lo4, stable };
    }


    std::tuple<Vector3f, Spectrum, Float>
    sample_stochastic(const PhaseFunctionContext &ctx,
                         const MediumInteraction3f &mei,
                         Float u_select, const Point2f &sample2,
                         const BilinearWeights &bw, Mask active) const {
        Vector4f cumw(bw.w[0],
                      bw.w[0] + bw.w[1],
                      bw.w[0] + bw.w[1] + bw.w[2],
                      bw.w[0] + bw.w[1] + bw.w[2] + bw.w[3]);

        UInt32 idx = bw.idx[0];
        idx = dr::select(u_select > cumw[0], bw.idx[1], idx);
        idx = dr::select(u_select > cumw[1], bw.idx[2], idx);
        idx = dr::select(u_select > cumw[2], bw.idx[3], idx);

        UInt32 s = dr::gather<UInt32>(m_grid_start, idx, active);
        UInt32 l = dr::gather<UInt32>(m_grid_len,   idx, active);
        Float  norm = dr::gather<Float>(m_norm, idx, active);

        Float cos_t = sample_cos_theta_entry(sample2.x(), s, l, norm, active);
        Float sin_t = dr::safe_sqrt(1.f - cos_t * cos_t);
        auto [sin_phi, cos_phi] = dr::sincos(2.f * dr::Pi<ScalarFloat> * sample2.y());
        Vector3f wo = -mei.to_world(Vector3f(sin_t * cos_phi, sin_t * sin_phi, cos_t));
        wo = dr::select(active, wo, Vector3f(0.f));

        auto [phase_val, pdf] = eval_pdf_flat(ctx, mei, wo, bw, active);

        Spectrum weight = dr::select(active, phase_val * dr::rcp(pdf), Spectrum(0.f));
        pdf              = dr::select(active, pdf, Float(0.f));

        return { wo, weight, pdf };
    }

    std::pair<BilinearWeights, Mask>
    get_weights(const MediumInteraction3f &mei, Mask active) const {
        Float r_eff = m_r_eff_volume->eval_1(mei, active);
        Float v_eff = m_v_eff_volume->eval_1(mei, active);

        Mask nan_mask = dr::isnan(r_eff) || dr::isnan(v_eff);
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

        Float r_min = 0.f, r_max = 0.f;
        if (m_n_r > 1) {
            r_min = dr::gather<Float>(m_r_eff_grid, UInt32(0u), active);
            r_max = dr::gather<Float>(m_r_eff_grid, UInt32(m_n_r - 1), active);
            if (dr::any(active && (r_eff < r_min || r_eff > r_max)))
                Throw("ParticlePhaseFunction: r_eff out of dataset range. %s out of [%s; %s]", r_eff, r_min, r_max);
        }

        Float v_min = 0.f, v_max = 0.f;
        if (m_n_v > 1) {
            v_min = dr::gather<Float>(m_v_eff_grid, UInt32(0u), active);
            v_max = dr::gather<Float>(m_v_eff_grid, UInt32(m_n_v - 1), active);
            if (dr::any(active && (v_eff < v_min || v_eff > v_max)))
                Throw("ParticlePhaseFunction: v_eff out of dataset range.");
        }

        auto weights = get_bilinear_weights(r_eff, v_eff, r_min, r_max, v_min, v_max, active);

        return { weights, active };
    }

    std::pair<Spectrum, Float>
    eval_pdf_flat_at(const PhaseFunctionContext &ctx,
                      const MediumInteraction3f &mei,
                      const Vector3f &wo,
                      const BilinearWeights &bw,
                      Float cos_t,
                      Vector4u starts, Vector4u lens, Vector4u lo4,
                      Mask active) const {
        Vector4u hi4        = lo4 + 1u;
        Vector4f x0_4       = dr::gather<Vector4f>(m_nodes, lo4, active);
        Vector4f x1_4       = dr::gather<Vector4f>(m_nodes, hi4, active);
        Vector4f t4         = (Vector4f(cos_t) - x0_4) / (x1_4 - x0_4);
        Vector4f norms      = dr::gather<Vector4f>(m_norm, bw.idx, active);

        Vector4f domain_lo  = dr::gather<Vector4f>(m_nodes, starts,             active);
        Vector4f domain_hi  = dr::gather<Vector4f>(m_nodes, starts + lens - 1u, active);
        Vector4f w          = dr::select(Vector4f(cos_t) >= domain_lo && Vector4f(cos_t) <= domain_hi,
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
                    mi[0]*scale,  mi[1]*scale, 0,            0,
                    mi[1]*scale,  mi[2]*scale, 0,            0,
                    0,            0,           mi[3]*scale,  mi[4]*scale,
                    0,            0,          -mi[4]*scale,  mi[5]*scale);
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
            Vector4f y0_m11 = dr::gather<Vector4f>(m_mueller, lo4 * 6u, active);
            Vector4f y1_m11 = dr::gather<Vector4f>(m_mueller, hi4 * 6u, active);
            Vector4f m11_4  = dr::fmadd(t4, y1_m11 - y0_m11, y0_m11);
            pdf = dr::dot(w * norms, m11_4) * dr::InvTwoPi<ScalarFloat>;
            phase_val = Spectrum(pdf);
        }

        return { phase_val, pdf };
    }

    std::pair<Spectrum, Float>
    eval_pdf_flat(const PhaseFunctionContext &ctx,
                  const MediumInteraction3f &mei,
                  const Vector3f &wo,
                  const BilinearWeights &bw,
                  Mask active) const {
        Float cos_t     = -dot(wo, mei.wi);
        Vector4u starts = dr::gather<Vector4u>(m_grid_start, bw.idx, active);
        Vector4u lens   = dr::gather<Vector4u>(m_grid_len,   bw.idx, active);
        Vector4u lo4    = ragged_searchsorted4(m_nodes, starts, lens, cos_t, active);
        return eval_pdf_flat_at(ctx, mei, wo, bw, cos_t, starts, lens, lo4, active);
    }

    template <size_t N>
    struct ColumnEntry {
        std::vector<ScalarFloat> x;
        std::array<std::vector<ScalarFloat>, N> y;
        ScalarFloat weight = 0.0;
    };

    template <size_t N>
    std::array<ColumnEntry<N>, 4>
    gather_column_entries(const BilinearWeights &bw, Mask active,
                           ScalarFloat &lo_bound, ScalarFloat &hi_bound) const {
        Vector4u starts = dr::gather<Vector4u>(m_grid_start, bw.idx, active);
        Vector4u lens   = dr::gather<Vector4u>(m_grid_len,   bw.idx, active);
        Vector4f norms  = dr::gather<Vector4f>(m_norm,       bw.idx, active);

        std::array<ColumnEntry<N>, 4> entries;
        for (size_t i = 0; i < 4; ++i) {
            size_t s = dr::slice(starts[i], 0);
            size_t l = dr::slice(lens[i],   0);

            ColumnEntry<N> &e = entries[i];
            e.weight = dr::slice(bw.w[i], 0) * dr::slice(norms[i], 0);
            e.x.resize(l);
            for (auto &col : e.y)
                col.resize(l);
            for (size_t k = 0; k < l; ++k) {
                e.x[k] = dr::slice(m_nodes, s + k);
                for (size_t c = 0; c < N; ++c)
                    e.y[c][k] = dr::slice(m_mueller, (s + k) * 6u + c);
            }

            ScalarFloat first = e.x.front(), last = e.x.back();
            if (i == 0) { lo_bound = first; hi_bound = last; }
            else        { lo_bound = std::max(lo_bound, first); hi_bound = std::min(hi_bound, last); }
        }
        return entries;
    }

    template <size_t N>
    static void merge_and_interpolate(const std::array<ColumnEntry<N>, 4> &entries,
                                       bool restrict_to_intersection,
                                       ScalarFloat lo_bound, ScalarFloat hi_bound,
                                       std::vector<ScalarFloat> &merged,
                                       std::array<std::vector<ScalarFloat>, N> &Y) {
        merged.clear();
        if (restrict_to_intersection)
            merged = { lo_bound, hi_bound };
        for (const ColumnEntry<N> &e : entries)
            for (ScalarFloat x : e.x)
                if (!restrict_to_intersection || (x >= lo_bound && x <= hi_bound))
                    merged.push_back(x);

        std::sort(merged.begin(), merged.end());
        merged.erase(std::unique(merged.begin(), merged.end()), merged.end());

        size_t n = merged.size();
        for (auto &col : Y)
            col.assign(n, ScalarFloat(0.0));

        for (const ColumnEntry<N> &e : entries) {
            if (e.weight == 0.0 || e.x.empty())
                continue;

            size_t l  = e.x.size();
            size_t hi = 1;
            for (size_t k = 0; k < n; ++k) {
                ScalarFloat x = merged[k];
                while (hi < l - 1 && e.x[hi] < x)
                    ++hi;
                size_t lo = hi - 1;

                ScalarFloat h = e.x[hi] - e.x[lo];
                ScalarFloat t = h > 0.0 ? (x - e.x[lo]) / h : ScalarFloat(0.0);
                t = std::clamp(t, ScalarFloat(0.0), ScalarFloat(1.0));

                for (size_t c = 0; c < N; ++c)
                    Y[c][k] += e.weight * (e.y[c][lo] + t * (e.y[c][hi] - e.y[c][lo]));
            }
        }
    }

    ref<PhaseFunction<Float, Spectrum>>
    reconstruct_tabulated(const BilinearWeights &bw, Mask active) const {
        ScalarFloat lo_bound = 0.0, hi_bound = 0.0;
        auto entries = gather_column_entries<6>(bw, active, lo_bound, hi_bound);

        std::vector<ScalarFloat> merged;
        std::array<std::vector<ScalarFloat>, 6> Y;
        merge_and_interpolate<6>(entries, false, lo_bound, hi_bound, merged, Y);

        auto serialize = [](const std::vector<ScalarFloat> &v, int precision) {
            std::ostringstream ss;
            ss << std::scientific << std::setprecision(precision);
            for (size_t k = 0; k < v.size(); ++k) {
                if (k) ss << ",";
                ss << v[k];
            }
            return ss.str();
        };

        Properties props("tabphase_polarized");
        props.set<std::string>("nodes", serialize(merged, 30));
        props.set<std::string>("m11", serialize(Y[0], 17));
        props.set<std::string>("m12", serialize(Y[1], 17));
        props.set<std::string>("m22", serialize(Y[2], 17));
        props.set<std::string>("m33", serialize(Y[3], 17));
        props.set<std::string>("m34", serialize(Y[4], 17));
        props.set<std::string>("m44", serialize(Y[5], 17));

        return PluginManager::instance()
            ->create_object<PhaseFunction<Float, Spectrum>>(props);
    }

    std::tuple<Vector3f, Spectrum, Float>
    sample(const PhaseFunctionContext &ctx,
           const MediumInteraction3f &mei,
           Float sample1,
           const Point2f &sample2,
           Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::PhaseFunctionSample, active);

        auto [bw, valid] = get_weights(mei, active);
        active = active && valid;

        Float cos_t;
        switch (m_blending_method) {
            case BlendingMethod::Stochastic:
                return sample_stochastic(ctx, mei, sample1, sample2, bw, active);

            case BlendingMethod::Reconstruct: {
                ref<PhaseFunction<Float, Spectrum>> tabulated = reconstruct_tabulated(bw, active);
                return tabulated->sample(ctx, mei, sample1, sample2, active);
            }

            case BlendingMethod::Bisect:
                cos_t = sample_cos_theta_bisect(sample2.x(), bw, active);
                break;
            case BlendingMethod::Tabulate:
                cos_t = sample_cos_theta_tabulate(sample2.x(), bw, active);
                break;

            case BlendingMethod::Search: {
                auto [cos_t_s, starts, lens, lo4, stable] = sample_cos_theta_search(sample2.x(), bw, active);
                Float sin_t = dr::safe_sqrt(1.f - cos_t_s * cos_t_s);
                auto [sin_phi, cos_phi] = dr::sincos(2.f * dr::Pi<ScalarFloat> * sample2.y());
                Vector3f wo{ sin_t * cos_phi, sin_t * sin_phi, cos_t_s };
                wo = -mei.to_world(wo);

                Vector4u lo4_fresh = ragged_searchsorted4(m_nodes, starts, lens, cos_t_s, active && !stable);
                Vector4u lo4_final = dr::select(stable, lo4, lo4_fresh);

                auto [phase_val, pdf] = eval_pdf_flat_at(ctx, mei, wo, bw, cos_t_s, starts, lens, lo4_final, active);
                Spectrum weight = phase_val * dr::rcp(pdf);

                wo     = dr::select(active, wo,     Vector3f(0.f));
                weight = dr::select(active, weight, Spectrum(0.f));
                pdf    = dr::select(active, pdf,    Float(0.f));
                return { wo, weight, pdf };
            }
        }
        Float sin_t = dr::safe_sqrt(1.f - cos_t * cos_t);
        auto [sin_phi, cos_phi] =
            dr::sincos(2.f * dr::Pi<ScalarFloat> * sample2.y());

        Vector3f wo{ sin_t * cos_phi, sin_t * sin_phi, cos_t };
        wo = -mei.to_world(wo);

        auto [phase_val, pdf] = eval_pdf_flat(ctx, mei, wo, bw, active);
        Spectrum weight = phase_val * dr::rcp(pdf);

        wo     = dr::select(active, wo,     Vector3f(0.f));
        weight = dr::select(active, weight, Spectrum(0.f));
        pdf    = dr::select(active, pdf,    Float(0.f));

        return { wo, weight, pdf };
    }

    std::pair<Spectrum, Float>
    eval_pdf(const PhaseFunctionContext &ctx,
             const MediumInteraction3f &mei,
             const Vector3f &wo,
             Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::PhaseFunctionEvaluate, active);

        auto [bw, valid] = get_weights(mei, active);
        active = active && valid;
        if (m_blending_method == BlendingMethod::Reconstruct) {
            ref<PhaseFunction<Float, Spectrum>> tabulated = reconstruct_tabulated(bw, active);
            return tabulated->eval_pdf(ctx, mei, wo, active);
        }

        auto [val, pdf] = eval_pdf_flat(ctx, mei, wo, bw, active);

        val = dr::select(active, val, Spectrum(0.f));
        pdf = dr::select(active, pdf, Float(0.f));

        return { val, pdf };
    }

    void traverse(TraversalCallback *cb) override {
        cb->put("r_eff_volume",  m_r_eff_volume.get(), ParamFlags::NonDifferentiable);
        cb->put("v_eff_volume",  m_v_eff_volume.get(), ParamFlags::NonDifferentiable);
        cb->put("r_eff_grid",    m_r_eff_grid,         ParamFlags::NonDifferentiable);
        cb->put("v_eff_grid",    m_v_eff_grid,         ParamFlags::NonDifferentiable);
        cb->put("nodes",           m_nodes,              ParamFlags::NonDifferentiable);
        cb->put("phase_mueller",   m_mueller,            ParamFlags::NonDifferentiable);
        cb->put("grid_start",      m_grid_start,         ParamFlags::NonDifferentiable);
        cb->put("grid_len",        m_grid_len,           ParamFlags::NonDifferentiable);
        cb->put("sigma_s_weight",  m_sigma_s_weight,     ParamFlags::NonDifferentiable);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "ParticlePhaseFunction[" << std::endl
            << "  n_r = " << m_n_r << "," << std::endl
            << "  n_v = " << m_n_v << "," << std::endl
            << "  total_pts = " << dr::width(m_nodes) << "," << std::endl
            << "  blending_method = " << method_to_string(m_blending_method) << std::endl
            << "]";
        return oss.str();
    }

    FloatStorage get_nodes() const override {
        size_t n = dr::width(m_nodes);
        std::vector<ScalarFloat> tmp(n);
        for (size_t i = 0; i < n; ++i)
            tmp[i] = dr::slice(m_nodes, i);

        std::sort(tmp.begin(), tmp.end());
        auto last = std::unique(tmp.begin(), tmp.end());
        tmp.erase(last, tmp.end());

        return dr::load<FloatStorage>(tmp.data(), tmp.size());
    }


    void eval_max(const FloatStorage &nodes, FloatStorage &values) const override {
        size_t n = dr::width(nodes);
        std::vector<ScalarFloat> vals(n, ScalarFloat(0));

        for (int ir = 0; ir < m_n_r; ++ir) {
            for (int iv = 0; iv < m_n_v; ++iv) {
                size_t entry = ir * m_n_v + iv;
                size_t start = dr::slice(m_grid_start, entry);
                size_t len   = dr::slice(m_grid_len,   entry);
                ScalarFloat norm = dr::slice(m_norm,      entry);

                size_t lo = start;

                for (size_t i = 0; i < n; ++i) {
                    ScalarFloat cos_t = dr::slice(nodes, i);

                    while (lo < start + len - 2 && dr::slice(m_nodes, lo + 1) < cos_t)
                        ++lo;

                    size_t hi     = lo + 1;
                    ScalarFloat x0  = dr::slice(m_nodes, lo);
                    ScalarFloat x1  = dr::slice(m_nodes, hi);
                    ScalarFloat t   = dr::clip((cos_t - x0) / (x1 - x0), ScalarFloat(0), ScalarFloat(1));
                    ScalarFloat y0  = dr::slice(m_mueller, lo * 6u);
                    ScalarFloat y1  = dr::slice(m_mueller, hi * 6u);
                    ScalarFloat m11 = dr::fmadd(t, y1 - y0, y0) * norm * dr::InvTwoPi<ScalarFloat>;
                    vals[i] = dr::maximum(vals[i], m11);
                }
            }
        }

        values = dr::load<FloatStorage>(vals.data(), vals.size());
    }

    MI_DECLARE_CLASS(ParticlePhaseFunction)

private:
    ref<Volume> m_r_eff_volume;
    ref<Volume> m_v_eff_volume;

    FloatStorage  m_r_eff_grid;
    FloatStorage  m_v_eff_grid;
    FloatStorage  m_nodes;
    FloatStorage  m_mueller;
    UInt32Storage m_grid_start;
    UInt32Storage m_grid_len;
    FloatStorage  m_cdf;
    FloatStorage  m_norm;
    FloatStorage  m_sigma_s_weight;

    int m_n_r = 0;
    int m_n_v = 0;

    BlendingMethod m_blending_method;
};

MI_EXPORT_PLUGIN(ParticlePhaseFunction)

NAMESPACE_END(mitsuba)
