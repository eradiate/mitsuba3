#include <array>
#include <drjit/dynamic.h>
#include <drjit/texture.h>
#include <mitsuba/core/distr_1d.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/quad.h>
#include <mitsuba/core/random.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/eradiate/oceanprops.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/ior.h>
#include <mitsuba/render/texture.h>
#include <tuple>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _plugin-bsdf-ocean_mishchenko:

Oceanic reflection model (:monosp:`ocean-mishchenko`)
-----------------------------------------------------

.. pluginparameters::

 * - wind_speed
   - |float|
   - :math:`k \in [0, 37.54]` m/s.
   - Specifies the wind speed at which to evaluate the oceanic reflectance
     (Default: :monosp:`0.1 m/s`).

 * - eta, k
   - |spectrum| or |texture|
   - Real and imaginary components of the water's index of refraction.
     (Default: :monosp:`1.33, 0.`)
   - |exposed|, |differentiable|, |discontinuous|

 * - ext_ior
   - |spectrum| or |texture|
   - Exterior index of refraction specified numerically or using a known
     material name. Note that the complex component is assumed to be 0 
     (Default: 1.000277).

 * - shininess
   - |float|
   - :math:`k \in [0, \infty]`.
   - Specifies the shininess which is used as the exponent for Blinn-Phong MIS
     (Default: :monosp:`50.`).

 * - shadowing
   - |bool|
   - Flag that indicates whether shadowing is calculated or not
     (Default: :monosp:`true`).

This plugin implements the polarized oceanic reflection model originally
implemented by :cite:`Mishchenko1997AerosolRetrievalPolarization`. This model
focuses on the sunglint reflectance, which follows the Cox and Munk surface
slope probability distribution. 

Note that this material is one-sided---that is, observed from the
back side, it will be completely black. If this is undesirable,
consider using the ``twosided`` BSDF adapter plugin.
The following snippet describes an oceanic surface material with monochromatic
parameters:

.. tab-set-code::

    .. code-block:: python

        "type": "ocean_mishchenko",
        "wind_speed": 10,
        "eta": 1.33,
        "k": 0.,
        "ext_ior": 1.0,
        "shininess": 50

    .. code-block:: xml

        <bsdf type="ocean_mishchenko">
            <float name="wind_speed" value="10"/>
            <float name="eta" value="1.33"/>
            <float name="k" value="0."/>
            <float name="ext_ior" value=1.0/>
            <float name="shininess" value="50"/>
        </bsdf>

.. note:: This model only implements the sunglint reflection. See :ref:`ocean
    legacy <plugin-bsdf-ocean_legacy>` for a bsdf that includes whitecap, sunglint,
    and underlight reflectance.

*/

template <typename Float, typename Spectrum>
class MishchenkoOceanBSDF final : public BSDF<Float, Spectrum> {
public:
    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture)

    using Complex2u = dr::Complex<UnpolarizedSpectrum>;

    /**
     * @brief Construct a new OceanBSDF object.
     *
     * @param props A set of properties to initialize the oceanic BSDF.
     */
    MishchenkoOceanBSDF(const Properties &props) : Base(props) {
        // Retrieve parameters
        m_wind_speed = props.get<ScalarFloat>("wind_speed", 0.1f);
        m_eta        = props.texture<Texture>("eta", 1.33f);
        m_k          = props.texture<Texture>("k", 0.f);
        m_ext_eta    = props.texture<Texture>("ext_ior", 1.000277f);
        m_shininess  = props.get<ScalarFloat>("shininess", 50.f);
        m_shadowing  = props.get<bool>("shadowing", true);

        update();

        // => Sun glint reflectance at the water surface is "specular"
        m_components.push_back(BSDFFlags::GlossyReflection |
                               BSDFFlags::FrontSide);

        // Set all the flags
        for (auto c : m_components)
            m_flags |= c;
        dr::set_attr(this, "flags", m_flags);
    }

    void traverse(TraversalCallback *callback) override {
        callback->put_parameter("wind_speed", m_wind_speed, +ParamFlags::Differentiable);
        callback->put_object("eta", m_eta, +ParamFlags::Differentiable);
        callback->put_object("k", m_k, +ParamFlags::Differentiable);
        callback->put_object("ext_ior", m_ext_eta, +ParamFlags::Differentiable);
    }

    void
    parameters_changed(const std::vector<std::string> & /*keys*/) override {
        update();
    }

    /**
     * @brief Update the variables that can be preprocessed.
     * This includes the complex index of refraction and mean square slope
     */
    void update() {

        // compute the index of refraction
        m_sigma2 = 0.5f*cox_munk_MMS(m_wind_speed);
    }

    std::pair<BSDFSample3f, Spectrum> sample(const BSDFContext &ctx,
                                             const SurfaceInteraction3f &si,
                                             Float /*sample1*/,
                                             const Point2f &sample2,
                                             Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);

        Float cos_theta_i = Frame3f::cos_theta(si.wi);
        active &= cos_theta_i > 0.f;

        BSDFSample3f bs = dr::zeros<BSDFSample3f>();
        Spectrum value(0.f);

        if (unlikely(dr::none_or<false>(active)) ||
            !ctx.is_enabled(BSDFFlags::GlossyReflection))
            return { bs, value };

        // For Blinn-Phong, we need to sample the half-vector
        Float ksi_1 = sample2.x(), ksi_2 = sample2.y();
        Float phi_h   = dr::TwoPi<Float> * ksi_1;
        Float theta_h = dr::acos(dr::pow(ksi_2, 1.f / (m_shininess + 2.f)));
        Vector3f half = dr::normalize(
            Vector3f(dr::sin(theta_h) * dr::cos(phi_h),
                     dr::sin(theta_h) * dr::sin(phi_h), dr::cos(theta_h)));

        Vector3f wo = 2.f * dr::dot(si.wi, half) * half - si.wi;

        // In the case of sampling the glint component, the outgoing
        // direction is sampled using the Blinn-Phong distribution.
        bs.wo                = wo;
        bs.sampled_component = 0;
        bs.sampled_type      = +BSDFFlags::GlossyReflection;

        bs.pdf = pdf(ctx, si, bs.wo, active);
        bs.eta = 1.f;

        Mask shadowing  = m_shadowing & active;
        Complex2u n_air(m_ext_eta->eval(si, active), 0.);
        Complex2u n_water(m_eta->eval(si, active), 
                          m_k->eval(si, active));

        // `TransportMode` has two states:
        //     - `Radiance`, trace from the sensor to the light sources
        //     - `Importance`, trace from the light sources to the sensor
        Vector3f wo_hat = ctx.mode == TransportMode::Radiance ? bs.wo : si.wi,
                 wi_hat = ctx.mode == TransportMode::Radiance ? si.wi : bs.wo;

        if constexpr (is_polarized_v<Spectrum>) {

            Spectrum val = eval_mishchenko(n_air, n_water, -wo_hat, wi_hat,
                                           m_sigma2, shadowing);
            value = val * Frame3f::cos_theta(wo_hat) / (bs.pdf * dr::Pi<Float>);

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
            value = mueller::rotate_mueller_basis(
                value, -wo_hat, p_axis_in, mueller::stokes_basis(-wo_hat),
                wi_hat, p_axis_out, mueller::stokes_basis(wi_hat));

        } else {
            value = eval_mishchenko(n_air, n_water, -wo_hat, wi_hat, m_sigma2,
                                    shadowing)[0][0];
            value *= Frame3f::cos_theta(wo_hat) / (bs.pdf * dr::Pi<Float>);
        }

        return { bs, value & (active && bs.pdf > 0.f) };
    }

    Spectrum eval(const BSDFContext &ctx, const SurfaceInteraction3f &si,
                  const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        bool has_glint = ctx.is_enabled(BSDFFlags::GlossyReflection);

        if (unlikely(dr::none_or<false>(active) || !has_glint))
            return 0.f;

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        // Ensure incoming and outgoing directions are in the upper hemisphere
        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

        // Compute the whitecap reflectance
        Spectrum value(0.f);

        Mask shadowing  = m_shadowing & active;
        Complex2u n_air(m_ext_eta->eval(si, active), 0.);
        Complex2u n_water(m_eta->eval(si, active), 
                          m_k->eval(si, active));

        // `TransportMode` has two states:
        //     - `Radiance`, trace from the sensor to the light sources
        //     - `Importance`, trace from the light sources to the sensor
        Vector3f wo_hat = ctx.mode == TransportMode::Radiance ? wo : si.wi,
                 wi_hat = ctx.mode == TransportMode::Radiance ? si.wi : wo;

        if constexpr (is_polarized_v<Spectrum>) {
            value = eval_mishchenko(n_air, n_water, -wo_hat, wi_hat, m_sigma2,
                                    shadowing);
            value *= Frame3f::cos_theta(wo_hat) * dr::InvPi<Float>;

            // Compute the Stokes reference frame vectors of this matrix that
            // are perpendicular to the plane of reflection.
            Vector3f n(0.f, 0.f, 1.f);

            /* The Stokes reference frame vector of this matrix lies in the
               scattering plane spanned by wi and wo. */
            Vector3f p_axis_in  = 
                dr::normalize(dr::cross(dr::normalize(dr::cross(n,-wo_hat)),-wo_hat));
            Vector3f p_axis_out = 
                dr::normalize(dr::cross(dr::normalize(dr::cross(n,wi_hat)),wi_hat));

            dr::masked(p_axis_in, dr::any(dr::isnan(p_axis_in)))   = 
                Vector3f(0.f, 1.f, 0.f);
            dr::masked(p_axis_out, dr::any(dr::isnan(p_axis_out))) = 
                Vector3f(0.f, 1.f, 0.f);

            // Rotate in/out reference vector of `value` s.t. it aligns with the
            // implicit Stokes bases of -wo_hat & wi_hat. */
            value = mueller::rotate_mueller_basis(
                value, -wo_hat, p_axis_in, mueller::stokes_basis(-wo_hat),
                wi_hat, p_axis_out, mueller::stokes_basis(wi_hat));
        } else {
            value = eval_mishchenko(n_air, n_water, -wo_hat, wi_hat, m_sigma2,
                                    shadowing)[0][0];
            value *= Frame3f::cos_theta(wo_hat) * dr::InvPi<Float>;
        }

        return value & active;
    }

    Float pdf(const BSDFContext &ctx, const SurfaceInteraction3f &si,
              const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        if (unlikely(!ctx.is_enabled(BSDFFlags::GlossyReflection) ||
                     dr::none_or<false>(active)))
            return 0.f;

        // Check if the normal has only zeros. If this is the case, use a
        // default normal
        Vector3f normal   = si.n;
        Mask degen_normal = dr::all(dr::eq(normal, Vector3f(0.f)));
        dr::masked(normal, degen_normal) = Vector3f(0.f, 0.f, 1.f);

        Vector3f half    = dr::normalize(si.wi + wo);
        Float projection = dr::dot(half, normal);
        Float D          = ((m_shininess + 2.f) /
                   dr::TwoPi<Float>) *dr::pow(projection, m_shininess);

        // We multiply the probability of the specular lobe with the pdf of
        // the Blinn-Phong distribution and the probability of the diffuse lobe
        // with the pdf of the cosine-weighted hemisphere.
        Float pdf_specular = (D * projection) / (4.0f * dr::dot(si.wi, half));

        Float pdf = pdf_specular;

        // If the outgoing direction is in the lower hemisphere, we return zero
        Float cos_theta_o = Frame3f::cos_theta(wo);

        return dr::select(cos_theta_o > 0.f, pdf, 0.f);
    }

    /**
     * @brief Evaluate ocean BRDF as implemented by Mishchenko.
     *
     * Evaluate the polarized BRDF of the ocean implemented by Mishchenko.
     * This model accounts for mean square slope probability, polarized
     * fresnel reflectance, and shadowing.
     *
     * @param n_ext Complex index of refraction outside of the medium.
     * @param n_water Complex index of refraction of the water.
     * @param wi Incident direction of light (physics convention).
     * @param wo Outgoing direction of light (physics convention).
     * @param sigma2 Mean square slope squared.
     * @param shadowing Flag that indicates the inclusion of shadowing.
     *
     * @return Ocean BPDF as a mueller matrix.
     */
    MuellerMatrix<UnpolarizedSpectrum> eval_mishchenko(
        const Complex2u &n_ext, const Complex2u &n_water, 
        Vector3f wi, Vector3f wo,
        const Float sigma2, Mask shadowing
    ) const {

        MuellerMatrix<UnpolarizedSpectrum> F =
            fresnel_sunglint_polarized(n_ext, n_water, wi, wo);

        Float mu_i = dr::abs(wi.z());
        Float mu_o = dr::abs(wo.z());

        dr::masked(wi, mu_i > 0.9999999f) =
            dr::normalize(wi + Vector3f(0.00001f, 0.f, 0.f));
        dr::masked(wo, mu_o > 0.9999999f) =
            dr::normalize(wo + Vector3f(0.00001f, 0.f, 0.f));
        dr::masked(mu_i, mu_i > 0.9999999f) = 0.9999999f;
        dr::masked(mu_o, mu_o > 0.9999999f) = 0.9999999f;

        // local surface normal k_d
        const Vector3f k_d    = wi - wo;
        const Float k_d_norm2 = dr::dot(k_d, k_d);

        // Slope probability distribution (Cox and Munk, 1955)
        Float coeff =
            1.f / (dr::pow(k_d.z(), 4.f) * 4.f * mu_i * mu_o * 2.f * sigma2);
        const Float exp = dr::exp(-(k_d.x() * k_d.x() + k_d.y() * k_d.y()) 
                                  / (k_d.z() * k_d.z() * 2.f * sigma2));
        coeff *= exp * (k_d_norm2 * k_d_norm2);

        F = F * coeff;

        // shadowing
        if (dr::any_or<true>(shadowing)) {

            Float shadow = 1.f / (1.f + gamma_shadowing(mu_i, sigma2) + gamma_shadowing(mu_o, sigma2));
            dr::masked(F, shadowing) = F * shadow;
        }

        return F;
    }

    /**
     * @brief Shadowing function
     *
     * @param mu absolute value of the cosine of theta.
     * @param sigma2 Mean square slope squared.
     */
    Float gamma_shadowing(const Float mu, const Float sigma2) const {
        Float cot = mu / dr::sqrt(1.f - mu * mu);
        Float result = 
            0.5f *
            (dr::sqrt(2.f * (1.f - mu * mu) / dr::Pi<Float>) 
             * dr::sqrt(sigma2) / mu 
             * dr::exp(-cot * cot / (2.f * sigma2)) 
             - (1.f - dr::erf(cot / dr::sqrt(2.f * sigma2))));
        return result;
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "OceanMishchenko[" << std::endl
            << "  wind_speed = " << string::indent(m_wind_speed) << std::endl
            << "  eta = " << string::indent(m_eta) << std::endl
            << "  k = " << string::indent(m_k) << std::endl
            << "  ext_ior = " << string::indent(m_ext_eta) << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS()
private:
    //  User-provided fields
    ScalarFloat m_wind_speed;
    ref<Texture> m_eta;
    ref<Texture> m_k;
    ref<Texture> m_ext_eta;
    ScalarFloat m_shininess;
    bool m_shadowing;

    ScalarFloat m_sigma2;
    OceanProperties<Float, Spectrum> m_ocean_props;
};

MI_IMPLEMENT_CLASS_VARIANT(MishchenkoOceanBSDF, BSDF)
MI_EXPORT_PLUGIN(MishchenkoOceanBSDF, "Mishchenko Ocean material")
NAMESPACE_END(mitsuba)
