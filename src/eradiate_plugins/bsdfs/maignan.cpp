#include <mitsuba/core/frame.h>
#include <mitsuba/core/fwd.h>
#include <mitsuba/core/math.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/eradiate/oceanprops.h>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _plugin-bsdf-maignan:

Polarized reflection model by Maignan et al. (:monosp:`maignan`)
--------------------------------------------------------

.. pluginparameters::

 * - rho_0
   - |spectrum| or |texture|
   - :math:`\rho_0 \ge 0`. Default: 0.1
   - |exposed| |differentiable|

 * - k
   - |spectrum| or |texture|
   - :math:`k \in \mathbb{R}`. Default: 0.1
   - |exposed| |differentiable|

 * - g
   - |spectrum| or |texture|
   - :math:`-1 \le g \le 1`. Default: 0.0
   - |exposed| |differentiable|

 * - rho_c
   - |spectrum| or |texture|
   - Default: Equal to `rho_0`
   - |exposed| |differentiable|

 * - ndvi
   - |spectrum| or |texture|
   - :math:`0 \le ndvi \le 1`. Default: 0.0
   - |exposed| |differentiable|

 * - refr_re
   - |spectrum| or |texture|
   - :math:`1 \le refr{_re} \le \infty`. Default: 1.5
   - |exposed| |differentiable|

 * - refr_im
   - |spectrum| or |texture|
   - :math:`0 \le refr{_im} \le \infty`. Default: 0.0
   - |exposed| |differentiable|

This plugin implements the reflection model proposed by
:cite:`Maignan2009PolarizedReflectance`.

Apart from floating point values, model parameters can be defined by nested or
referenced textures which are then mapped onto the shape based on its UV
parameterization (FIXMECE what does this mean?).

The model combines the RPV model for land surface reflection with the Fresnel equations to compute the polarization by reflection at land surfaces.

This plugin also supports the most common extension of the RPV model to four
parameters, namely the :math:`\rho_c` extension, as used in
:cite:`Widlowski2006Rayspread`.

For the fundamental formulae defining the RPV model, please refer to the
Eradiate Scientific Handbook.

Note that this material is one-sided---that is, observed from the
back side, it will be completely black. If this is undesirable,
consider using the ``twosided`` BSDF adapter plugin.
The following snippet describes an RPV material with monochromatic parameters:

.. tab-set-code::

    .. code-block:: python

        "type": "maignan",
        "rho_0": 0.02,
        "k": 0.3,
        "g": -0.12
        "ndvi": 0.0,
        "refr_re": 1.5,
        "refr_im": 0.0,
        "rho_c": 0.02

    .. code-block:: xml

        <bsdf type="maignan">
            <float name="rho_0" value="0.02"/>
            <float name="k" value="0.3"/>
            <float name="g" value="-0.12"/>
            <float name="ndvi" value="0.0"/>
            <float name="refr_re" value="1.5"/>
            <float name="refr_im" value="0.0"/>
            <float name="rho_c" value="0.02"/>
        </bsdf>
*/

MI_VARIANT
class MaignanBSDF final : public BSDF<Float, Spectrum> {
public:
    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture)

    using Complex2u = dr::Complex<UnpolarizedSpectrum>;
    
    MaignanBSDF(const Properties &props) : Base(props) {
        m_rho_0 = props.texture<Texture>("rho_0", 0.1f);
        m_g = props.texture<Texture>("g", 0.f);
        m_k = props.texture<Texture>("k", 0.1f);
        m_ndvi = props.texture<Texture>("ndvi", 0.f);
        m_refr_re = props.texture<Texture>("refr_re", 1.5f);
        m_refr_im = props.texture<Texture>("refr_im", 0.f);
        m_rho_c = props.texture<Texture>("rho_c", m_rho_0);
        m_flags = BSDFFlags::GlossyReflection | BSDFFlags::FrontSide;
        dr::set_attr(this, "flags", m_flags);
        m_components.push_back(m_flags);
    }

    void traverse(TraversalCallback *callback) override {
        callback->put_object("rho_0", m_rho_0.get(),
                             +ParamFlags::Differentiable);
        callback->put_object("g", m_g.get(), +ParamFlags::Differentiable);
        callback->put_object("k", m_k.get(), +ParamFlags::Differentiable);
        callback->put_object("rho_c", m_rho_c.get(),
                             +ParamFlags::Differentiable);
        callback->put_object("ndvi", m_ndvi.get(),
                             +ParamFlags::Differentiable);
        callback->put_object("refr_re", m_refr_re.get(),
                             +ParamFlags::Differentiable);
        callback->put_object("refr_im", m_refr_im.get(),
                             +ParamFlags::Differentiable);
    }

    std::pair<BSDFSample3f, Spectrum> sample(const BSDFContext &ctx,
                                             const SurfaceInteraction3f &si,
                                             Float /* position_sample */,
                                             const Point2f &direction_sample,
                                             Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);

        Float cos_theta_i = Frame3f::cos_theta(si.wi);
        BSDFSample3f bs = dr::zeros<BSDFSample3f>();

        active &= cos_theta_i > 0.f;
        if (unlikely(dr::none_or<false>(active) ||
                     !ctx.is_enabled(BSDFFlags::GlossyReflection)))
            return { bs, 0.f };

        bs.wo = warp::square_to_cosine_hemisphere(direction_sample);
        bs.pdf = warp::square_to_cosine_hemisphere_pdf(bs.wo);
        bs.eta = 1.f;
        bs.sampled_type = +BSDFFlags::GlossyReflection;
        bs.sampled_component = 0;

        UnpolarizedSpectrum ndvi = m_ndvi->eval(si, active);

        //CE copied from ocean_mishchenko probably needed with polarization
        // `TransportMode` has two states:
        //     - `Radiance`, trace from the sensor to the light sources
        //     - `Importance`, trace from the light sources to the sensor
        Vector3f wo_hat = ctx.mode == TransportMode::Radiance ? bs.wo : si.wi,
                 wi_hat = ctx.mode == TransportMode::Radiance ? si.wi : bs.wo;

        //CE: replace ref_re amd ref_im by eta, k as in bsdf-ocean-reflection?
        Complex2u n_air(bs.eta, 0.);
        // refractive index of "vegetation" or "surface", input to Maignan BPDF
        Complex2u n_water(m_refr_re->eval(si, active), m_refr_im->eval(si, active)); 
        
        Spectrum F;
        if constexpr (is_polarized_v<Spectrum>) {
          F = fresnel_sunglint_polarized(n_air, n_water, -wo_hat, wi_hat);

          //CE: Replace F11 by RPV BRDF, how to assess elements?
          
          /* The Stokes reference frame vector of this matrix lies in the
             meridian plane spanned by wi and n. */
          Vector3f n(0, 0, 1);
          Vector3f p_axis_in  =
            dr::normalize(dr::cross(dr::normalize(dr::cross(n,-wo_hat)),-wo_hat));
          Vector3f p_axis_out =
            dr::normalize(dr::cross(dr::normalize(dr::cross(n,wi_hat)),wi_hat));
          
          dr::masked(p_axis_in, dr::any(dr::isnan(p_axis_in))) =
            Vector3f(0.f, 1.f, 0.f);
          dr::masked(p_axis_out, dr::any(dr::isnan(p_axis_out))) =
            Vector3f(0.f, 1.f, 0.f);
          
          // Rotate in/out reference vector of `value` s.t. it aligns with the
          // implicit Stokes bases of -wo_hat & wi_hat. */
          F = mueller::rotate_mueller_basis(
              F, -wo_hat, p_axis_in, mueller::stokes_basis(-wo_hat),
              wi_hat, p_axis_out, mueller::stokes_basis(wi_hat));
          }
        else {
          F =  eval_rpv(si, bs.wo, active) * exp(-ndvi) *  Frame3f::cos_theta(bs.wo) / bs.pdf;
          
        }
        
        return   { bs, F & (active) };
        //return { bs, depolarizer<Spectrum>(value) & (active && bs.pdf > 0.f) };
    }

    /* Evaluation of the RPV BRDF (without foreshortening factor) as per the
       Eradiate Scientific Handbook. */
    UnpolarizedSpectrum eval_rpv(const SurfaceInteraction3f &si,
                                 const Vector3f &wo, Mask active) const {
        UnpolarizedSpectrum rho_0 = m_rho_0->eval(si, active);
        UnpolarizedSpectrum rho_c = m_rho_c->eval(si, active);
        UnpolarizedSpectrum g = m_g->eval(si, active);
        UnpolarizedSpectrum k = m_k->eval(si, active);
        //UnpolarizedSpectrum ndvi = m_ndvi->eval(si, active);
        //UnpolarizedSpectrum refr_re = m_refr_re->eval(si, active);
        //UnpolarizedSpectrum refr_im = m_refr_im->eval(si, active);

        auto [sin_phi_i, cos_phi_i] = Frame3f::sincos_phi(si.wi);
        auto [sin_phi_o, cos_phi_o] = Frame3f::sincos_phi(wo);
        Float cos_phi_i_minus_phi_o =
            cos_phi_i * cos_phi_o + sin_phi_i * sin_phi_o;
        Float sin_theta_i = Frame3f::sin_theta(si.wi);
        Float cos_theta_i = Frame3f::cos_theta(si.wi);
        Float tan_theta_i = Frame3f::tan_theta(si.wi);
        Float sin_theta_o = Frame3f::sin_theta(wo);
        Float cos_theta_o = Frame3f::cos_theta(wo);
        Float tan_theta_o = Frame3f::tan_theta(wo);

        // Henyey-Greenstein component
        Float cos_Theta = cos_theta_i * cos_theta_o +
                          sin_theta_i * sin_theta_o * cos_phi_i_minus_phi_o;
        // The following uses cos(pi-x) = -cos(x)
        UnpolarizedSpectrum F = (1.f - dr::sqr(g)) /
                     dr::pow((1.f + dr::sqr(g) + 2.f * g * cos_Theta), 1.5f);

        // Hot spot component
        Float G = dr::safe_sqrt(dr::sqr(tan_theta_i) + dr::sqr(tan_theta_o) -
                                2.f * tan_theta_i * tan_theta_o *
                                    cos_phi_i_minus_phi_o);
        UnpolarizedSpectrum H = 1.f + (1.f - rho_c) / (1.f + G);

        // Minnaert component
        UnpolarizedSpectrum M = dr::pow(
            cos_theta_i * cos_theta_o * (cos_theta_i + cos_theta_o), k - 1.f);

        
        // Total value 
        UnpolarizedSpectrum value = rho_0 * M * F * H * dr::InvPi<Float>;

        

        
        
        return value;
    }
    
    void debug_message() const {
        std::cout << "Debugging MaignanBSDF: " << to_string() << std::endl;
    }

    Spectrum eval(const BSDFContext & /*ctx*/, const SurfaceInteraction3f &si,
                  const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;
        Spectrum value = eval_rpv(si, wo, active);

        return dr::select(
            active, depolarizer<Spectrum>(value) * dr::abs(cos_theta_o), 0.f);
    }

    Float pdf(const BSDFContext & /* ctx */, const SurfaceInteraction3f &si,
              const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        Float pdf = warp::square_to_cosine_hemisphere_pdf(wo);

        return dr::select(cos_theta_i > 0.f && cos_theta_o > 0.f, pdf, 0.f);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "MaignanBSDF[" << std::endl
            << "  rho_0 = " << string::indent(m_rho_0) << "," << std::endl
            << "  g = " << string::indent(m_g) << "," << std::endl
            << "  k = " << string::indent(m_k) << "," << std::endl
            << "  ndvi = " << string::indent(m_ndvi) << "," << std::endl
            << "  refr_re = " << string::indent(m_refr_re) << "," << std::endl
            << "  refr_im = " << string::indent(m_refr_im);

        if (m_rho_0 != m_rho_c) {
            oss << "," << std::endl << "  rho_c = " << string::indent(m_rho_c);
        }
        oss << std::endl << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS()

private:
    ref<Texture> m_rho_0;
    ref<Texture> m_g;
    ref<Texture> m_k;
    ref<Texture> m_rho_c;
    ref<Texture> m_ndvi;
    ref<Texture> m_refr_re;
    ref<Texture> m_refr_im;
};

MI_IMPLEMENT_CLASS_VARIANT(MaignanBSDF, BSDF)
MI_EXPORT_PLUGIN(MaignanBSDF, "Maignan BSDF")
NAMESPACE_END(mitsuba)
