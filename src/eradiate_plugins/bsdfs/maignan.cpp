#include <mitsuba/core/frame.h>
#include <mitsuba/core/fwd.h>
#include <mitsuba/core/math.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/eradiate/oceanprops.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _plugin-bsdf-maignan:

Polarized reflection model by Maignan et al. (:monosp:`maignan`)
--------------------------------------------------------

.. pluginparameters::

 * - C
   - |spectrum| or |texture|
   - :math:`\rho_0 \ge 0`. Default: 5.0 
   - |exposed| |differentiable|

 * - ndvi
   - |spectrum| or |texture|
   - :math:`0 \le ndvi \le 1`. Default: 0.8
   - |exposed| |differentiable|

 * - refr_re
   - |spectrum| or |texture|
   - :math:`1 \le refr{_re} \le \infty`. Default: 1.5
   - |exposed| |differentiable|

 * - refr_im
   - |spectrum| or |texture|
   - :math:`0 \le refr{_im} \le \infty`. Default: 0.0
   - |exposed| |differentiable|

 * - ext_ior
   - |spectrum| or |texture|
   - Exterior index of refraction specified numerically or using a known
     material name. Note that the complex component is assumed to be 0
     (Default: 1.000277).

This plugin implements the reflection model proposed by
:cite:`Maignan2009PolarizedReflectance`.

Apart from floating point values, model parameters can be defined by nested or
referenced textures which are then mapped onto the shape based on its UV
parameterization (FIXMECE what does this mean?).

The model is based on fits to POLDER observations and combines a BRDF with the Fresnel 
reflectance matrix (see Eq. 21 in :cite:`Maignan2009PolarizedReflectance`).

Note that this material is one-sided---that is, observed from the
back side, it will be completely black. If this is undesirable,
consider using the ``twosided`` BSDF adapter plugin.
The following snippet describes an RPV material with monochromatic parameters:

.. tab-set-code::

    .. code-block:: python

        "type": "maignan",
        "C": 5.,
        "ndvi": 0.8,
        "refr_re": 1.5,
        "refr_im": 0.0,
        
    .. code-block:: xml

        <bsdf type="maignan">
            <float name="C" value="5"/>
            <float name="ndvi" value="0.8"/>
            <float name="refr_re" value="1.5"/>
            <float name="refr_im" value="0.0"/>
            <float name="ext_ior" value="1.0"/>
        </bsdf>
*/

MI_VARIANT
class MaignanBSDF final : public BSDF<Float, Spectrum> {
public:
    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture)

    using Complex2u = dr::Complex<UnpolarizedSpectrum>;

    MaignanBSDF(const Properties &props) : Base(props) {
        m_C   = props.texture<Texture>("C", 0.1f);
        m_ndvi    = props.texture<Texture>("ndvi", 0.f);
        m_refr_re = props.texture<Texture>("refr_re", 1.5f);
        m_refr_im = props.texture<Texture>("refr_im", 0.f);
        m_ext_eta = props.texture<Texture>("ext_ior", 1.000277f);
        m_flags   = BSDFFlags::GlossyReflection | BSDFFlags::FrontSide;
        dr::set_attr(this, "flags", m_flags);
        m_components.push_back(m_flags);
    }

    void traverse(TraversalCallback *callback) override {
        callback->put_object("C", m_C.get(),
                             +ParamFlags::Differentiable);
        callback->put_object("ndvi", m_ndvi.get(), +ParamFlags::Differentiable);
        callback->put_object("refr_re", m_refr_re.get(),
                             +ParamFlags::Differentiable);
        callback->put_object("refr_im", m_refr_im.get(),
                             +ParamFlags::Differentiable);
        callback->put_object("ext_ior", m_ext_eta, +ParamFlags::Differentiable);                   
    }

    

    std::pair<BSDFSample3f, Spectrum> sample(const BSDFContext &ctx,
                                             const SurfaceInteraction3f &si,
                                             Float /* position_sample */,
                                             const Point2f &direction_sample,
                                             Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);

        BSDFSample3f bs   = dr::zeros<BSDFSample3f>();
        Float cos_theta_i = Frame3f::cos_theta(si.wi);  
        Float sin_theta_i = Frame3f::sin_theta(si.wi);

        active &= cos_theta_i > 0.f;
        if (unlikely(dr::none_or<false>(active) ||
                     !ctx.is_enabled(BSDFFlags::GlossyReflection)))
            return { bs, 0.f };

        bs.wo           = warp::square_to_cosine_hemisphere(direction_sample);
        bs.pdf          = warp::square_to_cosine_hemisphere_pdf(bs.wo);
        bs.eta          = 1.f;
        bs.sampled_type = +BSDFFlags::GlossyReflection;
        bs.sampled_component = 0;

        
        // CE copied from ocean_mishchenko probably needed with polarization
        //  `TransportMode` has two states:
        //      - `Radiance`, trace from the sensor to the light sources
        //      - `Importance`, trace from the light sources to the sensor
        Vector3f wo_hat = ctx.mode == TransportMode::Radiance ? bs.wo : si.wi,
                 wi_hat = ctx.mode == TransportMode::Radiance ? si.wi : bs.wo;

        Complex2u n_air(bs.eta, 0.);
        // refractive index of "vegetation" or "surface", input to Maignan BPDF
        Complex2u n_surface(m_refr_re->eval(si, active),
                        m_refr_im->eval(si, active));
                        
        auto [sin_phi_i, cos_phi_i] = Frame3f::sincos_phi(si.wi);
        auto [sin_phi_o, cos_phi_o] = Frame3f::sincos_phi(bs.wo);
        Float cos_phi_i_minus_phi_o =
            cos_phi_i * cos_phi_o + sin_phi_i * sin_phi_o;
        Float sin_theta_o = Frame3f::sin_theta(bs.wo);
        Float cos_theta_o = Frame3f::cos_theta(bs.wo);
                    
        Float cos_Theta = cos_theta_i * cos_theta_o +
                        sin_theta_i * sin_theta_o * cos_phi_i_minus_phi_o;

        Float tan_alpha_i = dr::sqrt((1.0-cos_Theta)/(1.0+cos_Theta));
        
        // Fresnel matrix as implemented in MYSTIC
        // Float cos_alpha_i = dr::sqrt((1.0+cos_Theta)/2.0);                                        
        // Complex2u r1 = n_surface * n_surface * cos_alpha_i; 
        // Complex2u r2 = dr::sqrt(n_surface * n_surface - 1 + cos_alpha_i * cos_alpha_i);
        // Complex2u r_par= (r1 - r2) / (r1 + r2);
        // Complex2u r_perp = (cos_alpha_i - r2) / (cos_alpha_i +r2); 
        // Complex2u r_par_sq= dr::real(r_par*dr::conj(r_par));
        // Complex2u r_perp_sq = dr::real(r_perp*dr::conj(r_perp));    
        // Complex2u r_par_per1 = r_par*dr::conj(r_perp);
        // Complex2u r_par_per2 = dr::conj(r_par_per1);

        // const UnpolarizedSpectrum F00 = r_par_sq + r_perp_sq;
        // const UnpolarizedSpectrum F11 = F00; 
        // const UnpolarizedSpectrum F10 = r_par_sq - r_perp_sq;
        // const UnpolarizedSpectrum F01 = F10;
        // const UnpolarizedSpectrum F22 = dr::real(r_par_per1 + r_par_per2);
        // const UnpolarizedSpectrum F33 = F22;
        // const UnpolarizedSpectrum F23 = dr::real(r_par_per1 - r_par_per2);
        // const UnpolarizedSpectrum F32 = F23;
 
        //Factor of Maignan parameterization, only difference to ocean_bpdf I suppose 
        // Maignan & Breon 2009, Eq. 21
        UnpolarizedSpectrum C =  m_C->eval(si, active) * dr::exp(-tan_alpha_i) * dr::exp(-m_ndvi->eval(si, active)) /
                (4* (cos_theta_i + cos_theta_o));

        Spectrum F;

        if constexpr (is_polarized_v<Spectrum>) {   
                  
            F = fresnel_sunglint_polarized(n_air, n_surface, -wo_hat, wi_hat);      
            /* The Stokes reference frame vector of this matrix lies in the
               meridian plane spanned by wi and n. */
            Vector3f n(0, 0, 1);
            Vector3f p_axis_in = dr::normalize(
                dr::cross(dr::normalize(dr::cross(n, -wo_hat)), -wo_hat));
            Vector3f p_axis_out = dr::normalize(
                dr::cross(dr::normalize(dr::cross(n, wi_hat)), wi_hat));

            dr::masked(p_axis_in, dr::any(dr::isnan(p_axis_in))) =
                Vector3f(0.f, 1.f, 0.f);
            dr::masked(p_axis_out, dr::any(dr::isnan(p_axis_out))) =
                Vector3f(0.f, 1.f, 0.f);

            // Rotate in/out reference vector of `value` s.t. it aligns with the
            // implicit Stokes bases of -wo_hat & wi_hat. */
        
            F = mueller::rotate_mueller_basis(  
                F, -wo_hat, p_axis_in, mueller::stokes_basis(-wo_hat), wi_hat,
                p_axis_out, mueller::stokes_basis(wi_hat));
        } else {
           UnpolarizedSpectrum(fresnel_sunglint_polarized(n_air, n_surface, -wo_hat, wi_hat)[0][0]);
        }

        return { bs, C*F & (active) };
        
    }

    Spectrum eval(const BSDFContext &ctx, const SurfaceInteraction3f &si,
                  const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        //bool has_glint = ctx.is_enabled(BSDFFlags::GlossyReflection);

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);
        Float sin_theta_i = Frame3f::sin_theta(si.wi),
              sin_theta_o = Frame3f::sin_theta(wo);  

        // Ensure incoming and outgoing directions are in the upper hemisphere
        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

        // `TransportMode` has two states:
        //     - `Radiance`, trace from the sensor to the light sources
        //     - `Importance`, trace from the light sources to the sensor
        Vector3f wo_hat = ctx.mode == TransportMode::Radiance ? wo : si.wi,
                 wi_hat = ctx.mode == TransportMode::Radiance ? si.wi : wo;

        
        // Refractive index of air        
        Complex2u n_air(m_ext_eta->eval(si, active), 0.);
        // refractive index of "vegetation" or "surface", input to Maignan BPDF
        Complex2u n_surface(m_refr_re->eval(si, active),
                         m_refr_im->eval(si, active));

        Spectrum F;   
                     
        auto [sin_phi_i, cos_phi_i] = Frame3f::sincos_phi(si.wi);
        auto [sin_phi_o, cos_phi_o] = Frame3f::sincos_phi(wo);
        Float cos_phi_i_minus_phi_o =
            cos_phi_i * cos_phi_o + sin_phi_i * sin_phi_o;
                            
        Float cos_Theta =  cos_theta_i * cos_theta_o +
                        sin_theta_i * sin_theta_o * cos_phi_i_minus_phi_o;

        Float tan_alpha_i = dr::sqrt((1.0-cos_Theta)/(1.0+cos_Theta));  
        
        // Maignan & Breon 2009, Eq. 21
        UnpolarizedSpectrum C =  m_C->eval(si, active)* dr::exp(-tan_alpha_i) *dr::exp(-m_ndvi->eval(si, active)) /(4* (cos_theta_i + cos_theta_o));

                    
        if constexpr (is_polarized_v<Spectrum>) { 

            F = fresnel_sunglint_polarized(n_air, n_surface, -wo_hat, wi_hat);     
            /* The Stokes reference frame vector of this matrix lies in the
               meridian plane spanned by wi and n. */
            Vector3f n(0, 0, 1);
            Vector3f p_axis_in = dr::normalize(
                dr::cross(dr::normalize(dr::cross(n, -wo_hat)), -wo_hat));
            Vector3f p_axis_out = dr::normalize(
                dr::cross(dr::normalize(dr::cross(n, wi_hat)), wi_hat));

            dr::masked(p_axis_in, dr::any(dr::isnan(p_axis_in))) =
                Vector3f(0.f, 1.f, 0.f);
            dr::masked(p_axis_out, dr::any(dr::isnan(p_axis_out))) =
                Vector3f(0.f, 1.f, 0.f);

            // Rotate in/out reference vector of `value` s.t. it aligns with the
            // implicit Stokes bases of -wo_hat & wi_hat. */
        
            F = mueller::rotate_mueller_basis(  
                F, -wo_hat, p_axis_in, mueller::stokes_basis(-wo_hat), wi_hat,
                p_axis_out, mueller::stokes_basis(wi_hat));
        } else {
            UnpolarizedSpectrum(fresnel_sunglint_polarized(n_air, n_surface, -wo_hat, wi_hat)[0][0]);
        }
        
        // if constexpr (is_polarized_v<Spectrum>) {
        //     F[0][0] = C* F[0][0];
        // } else {
        //     F = C * F;
        // }
        return C* F & active;
         
    }                 
        

  
    void debug_message() const {
        std::cout << "Debugging MaignanBSDF: " << to_string() << std::endl;
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
            << "  C = " << string::indent(m_C) << "," << std::endl
            << "  ext_ior = " << string::indent(m_ext_eta) << "," << std::endl
            << "  ndvi = " << string::indent(m_ndvi) << "," << std::endl
            << "  refr_re = " << string::indent(m_refr_re) << "," << std::endl
            << "  refr_im = " << string::indent(m_refr_im);
        
        oss << std::endl << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS()

private:
    ref<Texture> m_C;
    ref<Texture> m_ndvi;
    ref<Texture> m_refr_re;
    ref<Texture> m_refr_im;
    ref<Texture> m_ext_eta;
};

MI_IMPLEMENT_CLASS_VARIANT(MaignanBSDF, BSDF)
MI_EXPORT_PLUGIN(MaignanBSDF, "Maignan BSDF")
NAMESPACE_END(mitsuba)
