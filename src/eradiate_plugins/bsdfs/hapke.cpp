#include <mitsuba/core/frame.h>
#include <mitsuba/core/fwd.h>
#include <mitsuba/core/math.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _plugin-bsdf-hapke:

Hapke surface model (:monosp:`hapke`)
-------------------------------------

The model uses 6 parameters:

 * - w
   - |spectrum| or |texture|
   - :math:`0 \le w \le 1`.
   - |exposed| |differentiable|

 * - b
   - |spectrum| or |texture|
   - :math:`0 \le w \le 1`.
   - |exposed| |differentiable|

 * - c
   - |spectrum| or |texture|
   - :math:`0 \le w \le 1`.
   - |exposed| |differentiable|

 * - theta (degree)
   - |spectrum| or |texture|
   - :math:`0 \le w \le 90`.
   - |exposed| |differentiable|

 * - B_0
   - |spectrum| or |texture|
   - :math:`0 \le w \le 1`.
   - |exposed| |differentiable|

 * - h
   - |spectrum| or |texture|
   - :math:`0 \le w \le 1`.
   - |exposed| |differentiable|

Implement the Hapke BSDF model as proposed by Bruce Hapke in 1984
(https://doi.org/10.1016/0019-1035(84)90054-X).

All parameters are required.

*/

// In particular, and for quick reference, equations 1 to 3, 14 to 18, 31 to 36,
// 45 to 58 from the reference paper are directly implemented in this module.

MI_VARIANT
class HapkeBSDF final : public BSDF<Float, Spectrum> {
public:
    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture)

    HapkeBSDF(const Properties &props) : Base(props) {

        m_w     = props.texture<Texture>("w");
        m_b     = props.texture<Texture>("b");
        m_c     = props.texture<Texture>("c");
        m_theta = props.texture<Texture>("theta");
        m_B_0   = props.texture<Texture>("B_0");
        m_h     = props.texture<Texture>("h");

        m_flags = BSDFFlags::GlossyReflection | BSDFFlags::FrontSide;
        m_components.push_back(m_flags);
    }

    void traverse(TraversalCallback *callback) override {
        callback->put_object("w", m_w.get(), +ParamFlags::Differentiable);
        callback->put_object("b", m_b.get(), +ParamFlags::Differentiable);
        callback->put_object("c", m_c.get(), +ParamFlags::Differentiable);
        callback->put_object("theta", m_theta.get(),
                             +ParamFlags::Differentiable);
        callback->put_object("B_0", m_B_0.get(), +ParamFlags::Differentiable);
        callback->put_object("h", m_h.get(), +ParamFlags::Differentiable);
    }

    std::pair<BSDFSample3f, Spectrum> sample(const BSDFContext &ctx,
                                             const SurfaceInteraction3f &si,
                                             Float /* position_sample */,
                                             const Point2f &direction_sample,
                                             Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);

        Float cos_theta_i = Frame3f::cos_theta(si.wi);
        BSDFSample3f bs   = dr::zeros<BSDFSample3f>();

        active &= cos_theta_i > 0.f;
        if (unlikely(dr::none_or<false>(active) ||
                     !ctx.is_enabled(BSDFFlags::GlossyReflection)))
            return { bs, 0.f };

        bs.wo           = warp::square_to_cosine_hemisphere(direction_sample);
        bs.pdf          = warp::square_to_cosine_hemisphere_pdf(bs.wo);
        bs.eta          = 1.f;
        bs.sampled_type = +BSDFFlags::GlossyReflection;
        bs.sampled_component = 0;

        UnpolarizedSpectrum value =
            eval_hapke(si, bs.wo, active) * Frame3f::cos_theta(bs.wo) / bs.pdf;

        return { bs, depolarizer<Spectrum>(value) & (active && bs.pdf > 0.f) };
    }

    MI_INLINE UnpolarizedSpectrum eval_H(const UnpolarizedSpectrum &w,
                                         const UnpolarizedSpectrum &x) const {
        const UnpolarizedSpectrum gamma = dr::sqrt(1.f - w);
        const UnpolarizedSpectrum ro    = (1.f - gamma) / (1.f + gamma);
        return 1.f / (1.f - w * x *
                                (ro + (1.f - 2.f * ro * x) * 0.5f *
                                          dr::log((1.f + x) / x)));
    }

    MI_INLINE UnpolarizedSpectrum
    eval_chi(const UnpolarizedSpectrum &tan_theta) const {
        return 1 / dr::sqrt(1.f + dr::Pi<ScalarFloat> * dr::square(tan_theta));
    }

    MI_INLINE Float eval_f(const Float &phi) const {
        const Float clamped_phi_div2 =
            dr::clip(phi / 2.f, 0.f, dr::Pi<Float> / 2.f - dr::Epsilon<Float>);
        return dr::exp(-2.f * dr::tan(clamped_phi_div2));
    }

    MI_INLINE UnpolarizedSpectrum eval_E1(const UnpolarizedSpectrum &tan_theta,
                                          const Float &x) const {
        return dr::exp(-2.f * dr::InvPi<Float> / tan_theta / dr::tan(x));
    }

    MI_INLINE UnpolarizedSpectrum eval_E2(const UnpolarizedSpectrum &tan_theta,
                                          const Float &x) const {
        return dr::exp(-dr::InvPi<Float> / dr::square(tan_theta) /
                       dr::square(dr::tan(x)));
    }

    MI_INLINE UnpolarizedSpectrum eval_eta_0e(
        const UnpolarizedSpectrum &chi, const Float &cos_i, const Float &sin_i,
        const UnpolarizedSpectrum &tan_theta, const UnpolarizedSpectrum &E2_i,
        const UnpolarizedSpectrum &E1_i) const {
        return chi * (cos_i + sin_i * tan_theta * E2_i / (2.f - E1_i));
    }

    MI_INLINE UnpolarizedSpectrum eval_eta_e(
        const UnpolarizedSpectrum &chi, const Float &cos_e, const Float &sin_e,
        const UnpolarizedSpectrum &tan_theta, const UnpolarizedSpectrum &E2_e,
        const UnpolarizedSpectrum &E1_e) const {
        return chi * (cos_e + sin_e * tan_theta * E2_e / (2.f - E1_e));
    }

    MI_INLINE UnpolarizedSpectrum eval_mu(const UnpolarizedSpectrum &tan_theta,
                                          const Float &e, const Float &i,
                                          const Float &cos_x,
                                          const Float &sin_x, const Float &phi,
                                          const Float &opt_cos_phi,
                                          const Float &sign) const {

        UnpolarizedSpectrum chi = eval_chi(tan_theta);

        UnpolarizedSpectrum E1_e = eval_E1(tan_theta, e);
        UnpolarizedSpectrum E1_i = eval_E1(tan_theta, i);
        UnpolarizedSpectrum E2_e = eval_E2(tan_theta, e);
        UnpolarizedSpectrum E2_i = eval_E2(tan_theta, i);

        Float sin_phi_div2 = dr::sin(phi * 0.5f);
        Float phi_div_pi   = phi * dr::InvPi<ScalarFloat>;

        return chi * (cos_x + sin_x * tan_theta *
                                  (opt_cos_phi * E2_e +
                                   sign * dr::square(sin_phi_div2) * E2_i) /
                                  (2.f - E1_e - phi_div_pi * E1_i));
    }

    MI_INLINE UnpolarizedSpectrum
    eval_mu_eG(const UnpolarizedSpectrum &tan_theta, const Float &e,
               const Float &i, const Float &phi, const Float &cos_phi) const {

        Float opt_cos_phi = dr::select(e <= i, cos_phi, 1.f);
        Float sign        = dr::select(e <= i, 1.f, -1.f);
        Float cos_e       = dr::cos(e);
        Float sin_e       = dr::sin(e);

        Float a = dr::select(e <= i, i, e);
        Float b = dr::select(e <= i, e, i);

        return eval_mu(tan_theta, a, b, cos_e, sin_e, phi, opt_cos_phi, sign);
    }

    MI_INLINE UnpolarizedSpectrum
    eval_mu_0eG(const UnpolarizedSpectrum &tan_theta, const Float &e,
                const Float &i, const Float &phi, const Float &cos_phi) const {

        Float opt_cos_phi = dr::select(e <= i, 1.f, cos_phi);
        Float sign        = dr::select(e <= i, -1.f, 1.f);
        Float cos_i       = dr::cos(i);
        Float sin_i       = dr::sin(i);

        Float a = dr::select(e <= i, i, e);
        Float b = dr::select(e <= i, e, i);

        return eval_mu(tan_theta, a, b, cos_i, sin_i, phi, opt_cos_phi, sign);
    }

    MI_INLINE UnpolarizedSpectrum
    eval_M(const UnpolarizedSpectrum &w, const UnpolarizedSpectrum &mu_0eG,
           const UnpolarizedSpectrum &mu_eG) const {
        // Multiple scattering function
        return eval_H(w, mu_0eG) * eval_H(w, mu_eG) - 1.f;
    }

    MI_INLINE Float eval_cos_g(const Float &mu_0, const Float &mu,
                               const Float &sin_i, const Float &sin_e,
                               const Float &cos_phi) const {
        return mu_0 * mu + sin_i * sin_e * cos_phi;
    }

    MI_INLINE UnpolarizedSpectrum
    eval_P(const UnpolarizedSpectrum &b, const UnpolarizedSpectrum &c,
           const UnpolarizedSpectrum &cos_g) const {
        // Fonction de phase P
        UnpolarizedSpectrum numerator = 1.f - dr::square(b);
        UnpolarizedSpectrum term1 =
            (1.f - c) * numerator /
            dr::pow(1 + 2 * b * cos_g + dr::square(b), 3.f / 2.f);
        UnpolarizedSpectrum term2 =
            c * numerator / dr::pow(1 - 2 * b * cos_g + dr::square(b), 3.f / 2.f);
        return term1 + term2;
    }

    MI_INLINE UnpolarizedSpectrum eval_S(
        const UnpolarizedSpectrum &eta_e, const UnpolarizedSpectrum &eta_0e,
        const UnpolarizedSpectrum &chi, const Float &e, const Float &i,
        const Float &mu, const UnpolarizedSpectrum &mu_0,
        const UnpolarizedSpectrum &mu_e, const UnpolarizedSpectrum &f) const {

        UnpolarizedSpectrum opt_mu  = dr::select(e < i, mu, mu_0);
        UnpolarizedSpectrum opt_eta = dr::select(e < i, eta_e, eta_0e);

        return (mu_e * mu_0 * chi) /
               (eta_e * eta_0e * (1.f - f + f * chi * opt_mu / opt_eta));
    }

    MI_INLINE UnpolarizedSpectrum eval_B(const UnpolarizedSpectrum &B_0,
                                         const UnpolarizedSpectrum &h,
                                         const UnpolarizedSpectrum &g) const {
        // Opposition effect
        return B_0 / (1.f + 1.f / h * dr::tan(g / 2));
    }

    const UnpolarizedSpectrum eval_hapke(const SurfaceInteraction3f &si,
                                         const Vector3f &wo,
                                         Mask active) const {

        UnpolarizedSpectrum theta = dr::deg_to_rad(m_theta->eval(si, active));
        UnpolarizedSpectrum tan_theta = dr::tan(theta);
        UnpolarizedSpectrum w         = m_w->eval(si, active);

        auto [sin_phi_e, cos_phi_e] = Frame3f::sincos_phi(wo);
        auto [sin_phi_i, cos_phi_i] = Frame3f::sincos_phi(si.wi);
        Float cos_phi = cos_phi_e * cos_phi_i + sin_phi_e * sin_phi_i;
        Float sin_e   = Frame3f::sin_theta(wo);
        Float mu      = Frame3f::cos_theta(wo);
        Float tan_e   = Frame3f::tan_theta(wo);
        Float sin_i   = Frame3f::sin_theta(si.wi);
        Float mu_0    = Frame3f::cos_theta(si.wi);
        Float tan_i   = Frame3f::tan_theta(si.wi);

        Float i      = dr::atan(tan_i);
        Float e      = dr::atan(tan_e);
        Float fr_phi = dr::safe_acos(cos_phi);
        Float phi    = dr::abs(dr::select(fr_phi > dr::Pi<Float>,
                                          2.f * dr::Pi<Float> - fr_phi, fr_phi));

        UnpolarizedSpectrum w_div_4 = w * 0.25f;

        UnpolarizedSpectrum mu_0eG = eval_mu_0eG(tan_theta, e, i, phi, cos_phi);
        UnpolarizedSpectrum mu_eG  = eval_mu_eG(tan_theta, e, i, phi, cos_phi);

        UnpolarizedSpectrum mu_ratio = mu_0eG / (mu_0eG + mu_eG) / mu_0;

        UnpolarizedSpectrum b = m_b->eval(si, active);
        UnpolarizedSpectrum c = m_c->eval(si, active);
        Float cos_g           = eval_cos_g(mu_0, mu, sin_i, sin_e, cos_phi);
        Float g               = dr::safe_acos(cos_g);

        UnpolarizedSpectrum P = eval_P(b, c, cos_g);

        UnpolarizedSpectrum B_0 = m_B_0->eval(si, active);
        UnpolarizedSpectrum h   = m_h->eval(si, active);
        UnpolarizedSpectrum B   = eval_B(B_0, h, g);

        UnpolarizedSpectrum M = eval_M(w, mu_0eG, mu_eG);

        UnpolarizedSpectrum f = eval_f(phi);

        UnpolarizedSpectrum chi = eval_chi(tan_theta);

        UnpolarizedSpectrum E1_e = eval_E1(tan_theta, e);
        UnpolarizedSpectrum E1_i = eval_E1(tan_theta, i);
        UnpolarizedSpectrum E2_e = eval_E2(tan_theta, e);
        UnpolarizedSpectrum E2_i = eval_E2(tan_theta, i);

        UnpolarizedSpectrum eta_0e =
            eval_eta_0e(chi, mu_0, sin_i, tan_theta, E2_i, E1_i);
        UnpolarizedSpectrum eta_e =
            eval_eta_e(chi, mu, sin_e, tan_theta, E2_e, E1_e);

        UnpolarizedSpectrum S =
            eval_S(eta_e, eta_0e, chi, e, i, mu, mu_0, mu_eG, f);

        Log(Trace, "mu ratio %s", mu_ratio);
        Log(Trace, "P %s", P);
        Log(Trace, "B %s", B);
        Log(Trace, "M %s", M);
        Log(Trace, "S %s", S);

        return w_div_4 * mu_ratio * (P * (1 + B) + M) * S;
    }

    Spectrum eval(const BSDFContext & /*ctx*/, const SurfaceInteraction3f &si,
                  const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;
        const Spectrum value = eval_hapke(si, wo, active);

        return dr::select(
            active, depolarizer<Spectrum>(value) * dr::abs(cos_theta_o), 0.f);
    }

    Float pdf(const BSDFContext &, const SurfaceInteraction3f &si,
              const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        Float pdf = warp::square_to_cosine_hemisphere_pdf(wo);

        return dr::select(cos_theta_i > 0.f && cos_theta_o > 0.f, pdf, 0.f);
    }

    std::pair<Spectrum, Float> eval_pdf(const BSDFContext &,
                                        const SurfaceInteraction3f &si,
                                        const Vector3f &wo,
                                        Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;
        Spectrum value = eval_hapke(si, wo, active);
        Float pdf      = warp::square_to_cosine_hemisphere_pdf(wo);

        return { depolarizer<Spectrum>(value) * dr::abs(cos_theta_o) & active,
                 dr::select(active, pdf, 0.f) };
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "HapkeBSDF["
            << "  w = " << string::indent(m_w) << "," << std::endl
            << "  b = " << string::indent(m_b) << "," << std::endl
            << "  c = " << string::indent(m_c) << "," << std::endl
            << "  theta = " << string::indent(m_theta) << "," << std::endl
            << "  B_0 = " << string::indent(m_B_0) << "," << std::endl
            << "  h = " << string::indent(m_h) << "," << std::endl;
        oss << std::endl << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS();

private:
    ref<Texture> m_w;
    ref<Texture> m_b;
    ref<Texture> m_c;
    ref<Texture> m_theta;
    ref<Texture> m_B_0;
    ref<Texture> m_h;
};

MI_IMPLEMENT_CLASS_VARIANT(HapkeBSDF, BSDF)
MI_EXPORT_PLUGIN(HapkeBSDF, "Hapke BSDF")
NAMESPACE_END(mitsuba)
