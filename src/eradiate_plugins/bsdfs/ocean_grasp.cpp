#include <array>
#include <drjit/dynamic.h>
#include <drjit/texture.h>
#include <mitsuba/core/distr_1d.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/random.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/microfacet.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/eradiate/oceanprops.h>
#include <tuple>

#define MI_ROUGH_TRANSMITTANCE_RES 64

NAMESPACE_BEGIN(mitsuba)

/**!

.. _plugin-bsdf-ocean_grasp:

(GRASP) Oceanic reflection model (:monosp:`ocean-grasp`)
-------------------------------------------------------------

.. pluginparameters::

 * - wavelength
   - |float|
   - Specifies the wavelength at which to evaluate the oceanic
     reflectance. Range: [200, 4000] nm. Required
   - |exposed|

 * - wind_speed
   - |spectrum| or |texture|
   - Specifies the wind speed at which to evaluate the oceanic
     reflectance. Range: [0, 37.54] m/s. Default: 0.1 m/s
   - |exposed| |differentiable|

 * - eta, k
   - |spectrum| or |texture|
   - Real and imaginary components of the water's index of refraction.
     Default: 1.33, 0.0
   - |exposed| |differentiable|

 * - ext_ior
   - |spectrum| or |texture|
   - Exterior index of refraction specified numerically or using a known
     material name. Note that the complex component is assumed to be 0.
     Default: 1.000277
   - |exposed| |differentiable|

 * - water_body_reflectance
   - |spectrum| or |texture|
   - Diffuse reflectance of radiations that entered and exited the water body.
     (Default: 0.).
   - |exposed| |differentiable|

 * - component
   - |int|
   - Debug: specifies which component of the oceanic reflection model to evaluate.
     Default: 0.

     * Component 0 is used to evaluate the total oceanic reflectance.
     * Component 1 evaluates the whitecap reflectance.
     * Component 2 evaluates the sun glint reflectance.
     * Component 3 evaluates the underlight reflectance.
     * Component 4 evaluates the whitecap and underlight reflectance together.

   - —

This plugin implements the oceanic reflection model originally detailed in
:cite:`Litvinov2024AerosolSurfaceCharacterization`. Note that this model
is monochromatic.

For the fundamental formulae defining the oceanic reflectance model, please
refer to the Eradiate Scientific Handbook.

Note that this material is one-sided---that is, observed from the
back side, it will be completely black. If this is undesirable,
consider using the ``twosided`` BSDF adapter plugin.
The following snippet describes an oceanic surface material with monochromatic
parameters:

.. tab-set-code::

    .. code-block:: python

        "type": "ocean_grasp",
        "wavelength": 550,
        "wind_speed": 10,
        "eta": 1.33,
        "k": 0.,
        "ext_ior": 1.0,
        "water_body_reflectance": 0.02,
        "component": 0,

    .. code-block:: xml

        <bsdf type="ocean_grasp">
            <float name="wavelength" value="550"/>
            <float name="wind_speed" value="10"/>
            <float name="eta" value="1.33"/>
            <float name="k" value="0."/>
            <float name="ext_ior" value=1.0/>
            <float name="water_body_reflectance" value=0.02/>
            <int name="component" value="0"/>
        </bsdf>
*/

template <typename Float, typename Spectrum>
class GRASPOceanBSDF final : public BSDF<Float, Spectrum> {
public:
    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture, MicrofacetDistribution)

    using Complex2u = dr::Complex<UnpolarizedSpectrum>;

    /**
     * @brief Construct a new OceanBSDF object.
     *
     * @param props A set of properties to initialize the oceanic BSDF.
     */
    GRASPOceanBSDF(const Properties &props) : Base(props) {
        // Retrieve parameters
        m_wavelength = props.get<ScalarFloat>("wavelength");
        m_eta        = props.texture<Texture>("eta", 1.33f);
        m_k          = props.texture<Texture>("k", 0.f);
        m_ext_ior    = props.texture<Texture>("ext_ior", 1.000277f);
        m_wind_speed = props.texture<Texture>("wind_speed", 0.1f);
        m_component  = props.get<ScalarInt32>("component", 0);
        m_water_body_reflectance =
            props.texture<Texture>("water_body_reflectance", 0.f);

        // Set the BSDF flags
        // => Whitecap and underlight reflectance is "diffuse"
        m_components.push_back(BSDFFlags::DiffuseReflection |
                               BSDFFlags::FrontSide);

        // => Sun glint reflectance at the water surface is "specular"
        m_components.push_back(BSDFFlags::GlossyReflection |
                               BSDFFlags::FrontSide);

        // Set all the flags
        for (auto c : m_components)
            m_flags |= c;

        parameters_changed();
    }

    void traverse(TraversalCallback *callback) override {
        callback->put_parameter("wavelength", m_wavelength,
                                +ParamFlags::NonDifferentiable);
        callback->put_object("wind_speed", m_wind_speed.get(),
                                +ParamFlags::Differentiable);
        callback->put_object("eta", m_eta.get(), +ParamFlags::Differentiable);
        callback->put_object("k", m_k.get(), +ParamFlags::Differentiable);
        callback->put_object("ext_ior", m_ext_ior.get(),
                                +ParamFlags::Differentiable);
        callback->put_object("water_body_reflectance",
                                m_water_body_reflectance.get(),
                                +ParamFlags::Differentiable);
    }

    void parameters_changed(
        const std::vector<std::string> & /*keys*/ = {}) override {

        /* Compute weights that further steer samples towards
        the specular or diffuse components */
        Float d_mean = m_water_body_reflectance->mean(), s_mean = 1.f;

        m_specular_sampling_weight = s_mean / (d_mean + s_mean);

        // Precompute rough reflectance (vectorized)
        // Retrieve slope RMS
        ScalarFloat wind_speed = dr::slice(m_wind_speed->mean());
        ScalarFloat sigma      = dr::slice(eval_sigma(Float(wind_speed)));

        // Evaluate the water and air index of refractions
        ScalarFloat eta = dr::slice(m_eta->mean()) / dr::slice(m_ext_ior->mean());

        using FloatX = DynamicBuffer<ScalarFloat>;
        using Vector3fX = Vector<FloatX, 3>;

        using FloatP = dr::Packet<dr::scalar_t<Float>>;
         mitsuba::MicrofacetDistribution<FloatP, Spectrum> distr(
            MicrofacetType::Beckmann,
            dr::SqrtTwo<ScalarFloat> * sigma,
            true
        );
        FloatX mu = dr::maximum(1e-6f, dr::linspace<FloatX>(0, 1, MI_ROUGH_TRANSMITTANCE_RES));
        FloatX zero = dr::zeros<FloatX>(MI_ROUGH_TRANSMITTANCE_RES);

        Vector3fX wi = Vector3fX(dr::sqrt(1 - mu * mu), zero, mu);

        auto external_transmittance = eval_transmittance(distr, wi, eta);

        m_external_transmittance = dr::load<DynamicBuffer<Float>>(
            external_transmittance.data(),
            dr::width(external_transmittance));

        m_internal_reflectance =
            dr::mean(eval_reflectance(distr, wi, 1.f / eta) * wi.z()) * 2.f;

        dr::make_opaque(m_specular_sampling_weight, m_internal_reflectance);
    }

    /**
     * @brief Evaluate the isotropic Root Mean Square Slope given a wind speed.
     * @param wind_speed Speed of wind at mast height [m/s].
     * @return Float Root Mean Square Slope.
     */
    Float eval_sigma(Float wind_speed) const {
        return dr::sqrt(0.5f * cox_munk_msslope_squared(wind_speed));
    }

    /**
     * @brief Evaluate the whitecap reflectance.
     *
     * Evaluates the whitecap reflectance at the provided wavelength and wind
     * speed.
     *
     * @param wind_speed Speed of wind at mast height [m/s].
     * @return Float The whitecap reflectance.
     */
    Float eval_whitecaps(Float wind_speed) const {
        return whitecap_reflectance_frouin(Float(m_wavelength), wind_speed);
    }

    /**
     * @brief Evaluate the sun glint reflectance.
     *
     * Evaluates the sun glint reflectance at the provided wavelength, incident
     * and outgoing directions, wind direction, wind speed, and chlorinity.
     *
     * @param si Surface interaction.
     * @param wi The incident direction of the light in graphics convention.
     * @param wo The outgoing direction of the light in graphics convention.
     * @return Float The sun glint reflectance.
     */
    Spectrum eval_glint(const SurfaceInteraction3f &si, const Vector3f &wi,
                        const Vector3f &wo, Mask active) const {

        Float cos_theta_i = Frame3f::cos_theta(wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        const Vector3f H = dr::normalize(wi + wo);

        // Retrieve slope RMS
        Float wind_speed = m_wind_speed->eval(si, active)[0];
        Float sigma      = eval_sigma(wind_speed);

        // Evaluate the water and air index of refractions
        Complex2u n_air(m_ext_ior->eval(si, active), 0.f);
        Complex2u n_water(m_eta->eval(si, active), m_k->eval(si, active));

        MicrofacetDistribution distr(MicrofacetType::Beckmann,
                                     dr::SqrtTwo<Float> * sigma, true);

        Float D = distr.eval(H);

        Spectrum F = dr::zeros<Spectrum>();
        if constexpr (is_polarized_v<Spectrum>) {
            F = fresnel_sunglint_polarized(n_air, n_water, -wi, wo);
        } else {
            F = fresnel_sunglint_polarized(n_air, n_water, -wi, wo)[0][0];
        }

        Float G = 1.f / (1.f + lambda(wi, sigma) + lambda(wo, sigma));
        dr::masked(G, dr::dot(wi, H) * cos_theta_i <= 0.f) = 0.f;
        dr::masked(G, dr::dot(wo, H) * cos_theta_o <= 0.f) = 0.f;

        Spectrum value = dr::Pi<Float> * F * D * G / (4.f * cos_theta_i * cos_theta_o);

        return value;
    }

    /**
     * @brief Smith shadowing function.
     *
     * @param v Incident or Outgoing vector.
     * @param sigma Root Mean square slope.
     */
    Float lambda(const Vector3f &v, const Float sigma) const {
        Float sigma_tan = sigma * dr::sqrt(1.f - v.z() * v.z()) / v.z();
        Float result =
            0.5f * (dr::sqrt(2.f / dr::Pi<Float>) * sigma_tan *
                        dr::exp(-dr::rcp(2.f * sigma_tan * sigma_tan)) -
                    (1.f - dr::erf(dr::rcp(dr::sqrt(2.f) * sigma_tan))));

        return result;
    }

    /**
     * @brief Evaluate the underwater light reflectance at the interaction
     * point.
     */
    Float eval_underlight(const SurfaceInteraction3f &si, Mask active) const {
        return m_water_body_reflectance->eval(si, active)[0];

    }

    Float lerp_gather(const DynamicBuffer<Float> &data, Float x, size_t size,
                      Mask active = true) const {
        using UInt32 = dr::uint32_array_t<Float>;
        x *= Float(size - 1);
        UInt32 index = dr::minimum(UInt32(x), uint32_t(size - 2));

        Float v0 = dr::gather<Float>(data, index, active),
              v1 = dr::gather<Float>(data, index + 1, active);

        return dr::lerp(v0, v1, x - Float(index));
    }

    std::pair<BSDFSample3f, Spectrum>
    sample(const BSDFContext &ctx, const SurfaceInteraction3f &si,
           Float sample1, const Point2f &sample2, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);

        bool has_diffuse  = ctx.is_enabled(BSDFFlags::DiffuseReflection, 0),
             has_specular = ctx.is_enabled(BSDFFlags::GlossyReflection, 1);

        Float cos_theta_i = Frame3f::cos_theta(si.wi);
        active &= cos_theta_i > 0.f;

        BSDFSample3f bs = dr::zeros<BSDFSample3f>();
        Spectrum result(0.f);

        if (unlikely(dr::none_or<false>(active)) ||
            (!has_diffuse && !has_specular))
            return { bs, 0.f };

        Float prob_whitecap       = Float(0.),
              prob_water_specular = Float(0.),
              prob_water_diffuse  = Float(0.);

        Float wind_speed = m_wind_speed->eval(si, active)[0];

        // Cater for case where only one lobe is activated
        if (unlikely(has_specular != has_diffuse)) {

            prob_water_specular = has_specular ? 1.f : 0.f;
            prob_water_diffuse  = 1.f - prob_water_specular,
            prob_whitecap       = 1.f - prob_water_specular;
        }
        else {
            // Determine which component should be sampled
            Float coverage   = whitecap_coverage_monahan(wind_speed);
            Float whitecap    = eval_whitecaps(wind_speed);
            Float t_i = lerp_gather(m_external_transmittance, cos_theta_i,
                            MI_ROUGH_TRANSMITTANCE_RES, active);

            prob_whitecap       = whitecap;
            prob_water_specular = (1.f - coverage) * (1.f - t_i) * m_specular_sampling_weight;
            prob_water_diffuse  = (1.f - coverage) * t_i * (1.f - m_specular_sampling_weight);

            Float sum = prob_whitecap + prob_water_diffuse + prob_water_specular;
            prob_whitecap /= sum;
            prob_water_diffuse /= sum;
            prob_water_specular /= sum;
        }

        Mask sample_diffuse  = active && (sample1 < (prob_whitecap + prob_water_diffuse)),
             sample_specular = active && !sample_diffuse;

        if (dr::any_or<true>(sample_diffuse)) {
            // In the case of sampling the diffuse component, the outgoing
            // direction is sampled from a cosine-weighted hemisphere.
            dr::masked(bs.wo, sample_diffuse) =
                warp::square_to_cosine_hemisphere(sample2);
            dr::masked(bs.sampled_component, sample_diffuse) = 0;
            dr::masked(bs.sampled_type, sample_diffuse) =
                +BSDFFlags::DiffuseReflection;
        }

        if (dr::any_or<true>(sample_specular)) {

            Float sigma = eval_sigma(wind_speed);
            MicrofacetDistribution distr(MicrofacetType::Beckmann,
                                         dr::SqrtTwo<Float> * sigma, true);

            auto [H, weight] = distr.sample(si.wi, sample2);
            Vector3f wo      = reflect(si.wi, H);

            dr::masked(bs.wo, sample_specular) = wo;

            dr::masked(bs.sampled_component, sample_specular) = 1;
            dr::masked(bs.sampled_type, sample_specular) =
                +BSDFFlags::GlossyReflection;
        }

        bs.pdf = pdf(ctx, si, bs.wo, active);
        bs.eta = 1.f;

        active &= bs.pdf > 0.f;
        result = eval(ctx, si, bs.wo, active);

        return { bs, result / bs.pdf & active };
    }

    Spectrum eval(const BSDFContext &ctx, const SurfaceInteraction3f &si,
                  const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        bool has_diffuse  = ctx.is_enabled(BSDFFlags::DiffuseReflection, 0),
             has_specular = ctx.is_enabled(BSDFFlags::GlossyReflection, 1);

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        // Ensure incoming and outgoing directions are in the upper hemisphere
        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

        if (unlikely(dr::none_or<false>(active) ||
                     (!has_specular && !has_diffuse)))
            return 0.f;

        // Compute the whitecap reflectance
        Spectrum result(0.f);

        // 6SV considers the incident direction to come from the light source
        //     - `Radiance`, trace from the sensor to the light sources
        //     - `Importance`, trace from the light sources to the sensor
        Vector3f wo_hat = ctx.mode == TransportMode::Radiance ? wo : si.wi,
                 wi_hat = ctx.mode == TransportMode::Radiance ? si.wi : wo;

        // Combine the results

        // Make reflectances available for debug purposes
        UnpolarizedSpectrum whitecap_reflectance(0.f);
        UnpolarizedSpectrum underlight_reflectance(0.f);
        Spectrum glint_reflectance(0.f);

        Float wind_speed = m_wind_speed->eval(si, active)[0];
        Float coverage   = whitecap_coverage_monahan(wind_speed);

        if (has_diffuse) {
            Float t_i = lerp_gather(m_external_transmittance, cos_theta_i,
                                    MI_ROUGH_TRANSMITTANCE_RES, active),
                  t_o = lerp_gather(m_external_transmittance, cos_theta_o,
                                    MI_ROUGH_TRANSMITTANCE_RES, active);
            UnpolarizedSpectrum n_air = m_ext_ior->eval(si, active);
            UnpolarizedSpectrum n_water = m_eta->eval(si, active);
            UnpolarizedSpectrum m_inv_eta_2 = dr::square(n_air / n_water);

            whitecap_reflectance   = eval_whitecaps(wind_speed);
            underlight_reflectance = eval_underlight(si, active);

            underlight_reflectance /= 1.f - UnpolarizedSpectrum(m_internal_reflectance);
            underlight_reflectance *= (m_inv_eta_2 * t_i * t_o);

            // Diffuse scattering implies no polarization.
            result += depolarizer<Spectrum>(whitecap_reflectance +
                                            (1.f - coverage) *
                                                underlight_reflectance);
        }

        if (has_specular) {

            if constexpr (is_polarized_v<Spectrum>) {
                // If sun glint is enabled, compute the glint reflectance
                glint_reflectance = eval_glint(si, wo_hat, wi_hat, active);

                /* The Stokes reference frame vector of this matrix lies in the
                meridian plane spanned by wi and n. */
                Vector3f n(0.f, 0.f, 1.f);
                Vector3f p_axis_in = dr::normalize(
                    dr::cross(dr::normalize(dr::cross(n, -wo_hat)), -wo_hat));
                Vector3f p_axis_out = dr::normalize(
                    dr::cross(dr::normalize(dr::cross(n, wi_hat)), wi_hat));

                dr::masked(p_axis_in, dr::any(dr::isnan(p_axis_in))) =
                    Vector3f(0.f, 1.f, 0.f);
                dr::masked(p_axis_out, dr::any(dr::isnan(p_axis_out))) =
                    Vector3f(0.f, 1.f, 0.f);

                // Rotate in/out reference vector of `value` s.t. it aligns with
                // the implicit Stokes bases of -wo_hat & wi_hat. */
                glint_reflectance = mueller::rotate_mueller_basis(
                    glint_reflectance, -wo_hat, p_axis_in,
                    mueller::stokes_basis(-wo_hat), wi_hat, p_axis_out,
                    mueller::stokes_basis(wi_hat));

            } else {
                // If sun glint is enabled, compute the glint reflectance
                glint_reflectance = eval_glint(si, wo_hat, wi_hat, active);
            }

            result += (1.f - coverage) * glint_reflectance;
        }

        dr::masked(result, active) *= cos_theta_o * dr::InvPi<Float>;

        // For debugging purposes, channel indicates which BRDF term to evaluate
        switch (m_component) {
            case 1:
                result[active] = depolarizer<Spectrum>(whitecap_reflectance);
                break;
            case 2:
                result[active] = (1.f - coverage) * glint_reflectance;
                break;
            case 3:
            {
                result[active] = depolarizer<Spectrum>((1.f - coverage) *
                                                       underlight_reflectance);
                break;

            }
            case 4:
            {
                result[active] = depolarizer<Spectrum>(
                    whitecap_reflectance +
                    (1.f - coverage) * underlight_reflectance);
                break;
            }
            default:
                break;
        }

        return result & active;
    }

    Float pdf(const BSDFContext &ctx, const SurfaceInteraction3f &si,
              const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        bool has_diffuse  = ctx.is_enabled(BSDFFlags::DiffuseReflection, 0),
             has_specular = ctx.is_enabled(BSDFFlags::GlossyReflection, 1);

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

        if (unlikely((!has_diffuse && !has_specular) ||
                     dr::none_or<false>(active)))
            return 0.f;

        Float prob_whitecap       = Float(0.),
              prob_water_specular = Float(0.),
              prob_water_diffuse  = Float(0.);

        Float wind_speed = m_wind_speed->eval(si, active)[0];
        Float sigma = eval_sigma(wind_speed);

        // Cater for case where only one lobe is activated
        if (unlikely(has_specular != has_diffuse)) {

            prob_water_specular = has_specular ? 1.f : 0.f;
            prob_water_diffuse  = 1.f - prob_water_specular,
            prob_whitecap       = 1.f - prob_water_specular;
        }
        else {
            // Determine which component should be sampled
            Float coverage   = whitecap_coverage_monahan(wind_speed);
            Float whitecap   = eval_whitecaps(wind_speed);
            Float t_i = lerp_gather(m_external_transmittance, cos_theta_i,
                            MI_ROUGH_TRANSMITTANCE_RES, active);

            prob_whitecap       = whitecap;
            prob_water_specular = (1.f - coverage) * (1.f - t_i) * m_specular_sampling_weight;
            prob_water_diffuse  = (1.f - coverage) * t_i * (1.f - m_specular_sampling_weight);

            Float sum = prob_whitecap + prob_water_diffuse + prob_water_specular;
            prob_whitecap /= sum;
            prob_water_diffuse /= sum;
            prob_water_specular /= sum;
        }

        // Weight the diffuse component by the cosine hemisphere pdf
        Float prob_cosine = warp::square_to_cosine_hemisphere_pdf(wo);
        prob_whitecap      *= prob_cosine;
        prob_water_diffuse *= prob_cosine;

        // Weight the specular component using the visible beckmann pdf
        Vector3f H = dr::normalize(wo + si.wi);

        MicrofacetDistribution distr(MicrofacetType::Beckmann,
                                     dr::SqrtTwo<Float> * sigma, true);
        // following roughplastic definition, not sure why it differs from the
        // microfacter
        prob_water_specular *=
            distr.eval(H) * distr.smith_g1(si.wi, H) / (4.f * cos_theta_i);

        // bring the three together.
        Float result = prob_whitecap + prob_water_diffuse + prob_water_specular;

        return dr::select(active, result, 0.f);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "GRASPOcean[" << std::endl
            << "  component = " << string::indent(m_component) << ","
            << std::endl
            << "  wavelength = " << string::indent(m_wavelength) << ","
            << std::endl
            << "  wind_speed = " << string::indent(m_wind_speed) << ","
            << std::endl
            << "  eta = " << string::indent(m_eta) << "," << std::endl
            << "  k = " << string::indent(m_k) << "," << std::endl
            << "  ext_eta = " << string::indent(m_ext_ior) << "," << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS()
private:
    //  User-provided fields
    ScalarInt32 m_component;
    ScalarFloat m_wavelength;
    Float m_specular_sampling_weight;
    ref<Texture> m_wind_speed;
    ref<Texture> m_eta;
    ref<Texture> m_k;
    ref<Texture> m_ext_ior;
    ref<Texture> m_water_body_reflectance;

    DynamicBuffer<Float> m_external_transmittance;
    Float m_internal_reflectance;
};

MI_IMPLEMENT_CLASS_VARIANT(GRASPOceanBSDF, BSDF)
MI_EXPORT_PLUGIN(GRASPOceanBSDF, "GRASP Ocean material")
NAMESPACE_END(mitsuba)
