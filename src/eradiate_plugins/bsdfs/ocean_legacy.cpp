#include <array>
#include <drjit/dynamic.h>
#include <drjit/texture.h>
#include <mitsuba/core/distr_1d.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/quad.h>
#include <mitsuba/core/random.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/eradiate/oceanprops.h>
#include <tuple>

// Transmittance angular resolution
// i.e. texture resolution
#define MI_OCEAN_TRANSMITTANCE_RES 64

NAMESPACE_BEGIN(mitsuba)

/**!

.. _plugin-bsdf-ocean_legacy:

(Legacy 6S) Oceanic reflection model (:monosp:`ocean-legacy`)
-------------------------------------------------------------

.. pluginparameters::
 :extra-rows: 8

 * - wavelength
   - |float|
   - :math:`k \in [200, 4000]` nm.
   - Specifies the wavelength at which to evaluate the oceanic reflectance.

 * - wind_speed
   - |float|
   - :math:`k \in [0, 37.54]` m/s.
   - Specifies the wind speed at which to evaluate the oceanic reflectance
     (Default: :monosp:`0.1 m/s`).

 * - wind_direction
   - |float|
   - :math:`k \in [0, 360]` deg.
   - Specifies the wind direction at which to evaluate the oceanic reflectance
     in North Left convention (Default: :monosp:`0. deg`).

 * - chlorinity
   - |float|
   - Specifies the chlorinity of the water at which to evaluate the oceanic
     reflectance (Default: :monosp:`19. g/kg`).

 * - pigmentation
   - |float|
   - :math:`k \in [0.3, \infty]`.
   - Specifies the pigmentation of the water at which to evaluate the oceanic
     reflectance (Default: :monosp:`0.3 mg/m^3`).

 * - shininess
   - |float|
   - :math:`k \in [0, \infty]`.
   - Specifies the shininess which is used as the exponent for Blinn-Phong MIS
     (Default: :monosp:`50.`). If set to -1, the glint component is sampled
     using a cosine-hemisphere strategy.

 * - component
   - |int|
   - Debug: specifies which component of the oceanic reflection model to evaluate.
     Default: 0 Component 0 is used to evaluate the total oceanic reflectance.
     Component 1 evaluates the whitecap reflectance. Component 2 evaluates the sun
     glint reflectance. Component 3 evaluates the underlight reflectance.
     Component 4 evaluates the whitecap and underlight reflectance together.

 * - coverage
   - |float|
   - Fraction of the surface occupied by whitecaps. Modifying this parameter has
     no effect: it is automatically computed from the wind speed.

This plugin implements the oceanic reflection model originally
implemented in the 6S radiative transfer model. Note that this model
is monochromatic.

For the fundamental formulae defining the oceanic reflectance model, please
refer to the Eradiate Scientific Handbook.

Note that this material is one-sided---that is, observed from the
back side, it will be completely black. If this is undesirable,
consider using the ``twosided`` BSDF adapter plugin.
The following snippet describes an oceanic surface material with monochromatic
parameters:

.. warning::
    The wind direction is given in degrees and follows the North Left convention
    as in 6SV.

.. tab-set-code::

    .. code-block:: python

        "type": "ocean_legacy",
        "wavelength": 550,
        "wind_speed": 10,
        "wind_direction": 0,
        "chlorinity": 19,
        "pigmentation": 0.3,
        "shininess": 50
        "component": 0,

    .. code-block:: xml

        <bsdf type="ocean_legacy">
            <float name="wavelength" value="550"/>
            <float name="wind_speed" value="10"/>
            <float name="wind_direction" value="0"/>
            <float name="chlorinity" value="19"/>
            <float name="pigmentation" value="0.3"/>
            <float name="shininess" value="50"/>
            <int name="component" value="0"/>
        </bsdf>
*/

template <typename Float, typename Spectrum> class OceanUtilities {
public:
    MI_IMPORT_TYPES()

    using Complex2u = dr::Complex<UnpolarizedSpectrum>;

    /**
     * @brief Construct a new Ocean Utilities object.
     *
     * Construct a new Ocean Utilities object and initializes the ocean
     * properties.
     */
    OceanUtilities() : m_ocean_props() {}

    /**
     * @brief Update the cox-munk variables and R_omega,
     * both depend on scene's properties only and only need
     * to be called on construct and parameter update.
     **/
    void update(const ScalarFloat &wavelength, const ScalarFloat &wind_speed,
                const ScalarFloat &pigmentation, const bool &shadowing) {
        // Update Cox-Munk variables
        // std::tie(m_sigma_c, m_sigma_u, std::ignore) = mean_square_slope_cox_munk(wind_speed);
        std::tie(m_sigma_c, m_sigma_u) = cox_munk_crosswind_upwind(wind_speed);
        m_sigma_c = dr::sqrt(m_sigma_c);
        m_sigma_u = dr::sqrt(m_sigma_u);

        m_c_21    = 0.01f - 0.0086f * wind_speed;
        m_c_03    = 0.04f - 0.033f * wind_speed;

        m_r_omega = r_omega(wavelength, pigmentation);

        m_shadowing = shadowing;
    }

    /**
     * @brief Evaluate the index of refraction of the ocean.
     *
     * @param wavelength The wavelength at which to evaluate the reflectance.
     * @param chlorinity The chlorinity of the water.
     * @return pair(ScalarFloat,ScalarFloat)
     */
    std::pair<ScalarFloat, ScalarFloat>
    eval_index_refraction(const ScalarFloat &wavelength,
                          const ScalarFloat &chlorinity) const {
        auto [n_real, n_imag] =
            water_ior<Float, Spectrum, ScalarFloat>(m_ocean_props,wavelength,chlorinity);
        return { n_real, n_imag };
    }

    /**
     * @brief Evaluate the fractional coverage of whitecaps.
     *
     * Evaluates the fractional coverage of whitecaps at the given wind speed,
     * using the Monahan et al. (1986) model. The coverage is clamped to the
     * range [0, 1] (i.e. wind speed can be within the range [0, 37.54]).
     *
     * @param wind_speed The wind speed at which to evaluate the coverage.
     * @return ScalarFloat The fractional coverage of whitecaps.
     */
    ScalarFloat eval_whitecap_coverage(const ScalarFloat &wind_speed) const {
        return dr::clamp(m_monahan_alpha *
                         dr::pow(wind_speed, m_monahan_lambda),
                         0.0f, 1.0f);
    }

    /**
     * @brief Evaluate the reflectance of whitecaps.
     *
     * Evaluates the reflectance of whitecaps at the given wavelength and wind
     * speed. The reflectance is computed as the product of the effective
     * reflectance of whitecaps and the fractional coverage of whitecaps.
     *
     * @param wavelength The wavelength at which to evaluate the reflectance.
     * @param wind_speed The wind speed at which to evaluate the reflectance.
     * @return ScalarFloat The reflectance of whitecaps.
     */
    ScalarFloat eval_whitecaps(const ScalarFloat &wavelength,
                               const ScalarFloat &wind_speed) const {
        // Compute the fractional coverage of whitecaps
        ScalarFloat coverage = eval_whitecap_coverage(wind_speed);

        // Proper interpolation of the effective reflectance
        ScalarFloat eff_reflectance =
            m_ocean_props.effective_reflectance(wavelength);

        // Compute the whitecap reflectance
        ScalarFloat whitecap_reflectance = coverage * eff_reflectance;

        return whitecap_reflectance;
    }

    /**
     * @brief Evaluate the sun glint reflectance.
     *
     * Evaluates the sun glint reflectance at the given wavelength, incident and
     * outgoing directions, and wind direction. The effects of wind speed and
     * chlorinity are precomputed in update(). The reflectance is computed using
     * the Cox-Munk distribution, the Fresnel equations, and relative tilt of
     * the oceanic surface.
     *
     * @param n_real The real part of the index of refraction of water.
     * @param n_imag The imaginary part of the index of refraction of water.
     * @param wi The incident direction of the light in graphics convention.
     * @param wo The outgoing direction of the light in graphics convention.
     * @param wind_direction The direction of the wind relative to the incident
     * light.
     * @return Float The sun glint reflectance.
     */
    Spectrum eval_glint(const ScalarFloat &n_real, const ScalarFloat &n_imag,
                     const Vector3f &wi, const Vector3f &wo,
                     const ScalarFloat &wind_direction) const {
        // Transform directions into azimuthal and zenithal angles
        Float s_theta_i = dr::sqrt(1.f - wi.z() * wi.z());
        Float c_theta_i = wi.z();

        Float s_theta_o = dr::sqrt(1.f - wo.z() * wo.z());
        Float c_theta_o = wo.z();

        Float phi_i         = dr::atan2(wi.y(), wi.x());
        Float phi_o         = dr::atan2(wo.y(), wo.x());
        Float phi           = phi_i - phi_o;
        auto [s_phi, c_phi] = dr::sincos(phi);

        Float phi_w             = phi_i - wind_direction;
        auto [s_phi_w, c_phi_w] = dr::sincos(phi_w);

        Spectrum value = dr::zeros<Spectrum>();

        if constexpr (is_polarized_v<Spectrum>){
            value = eval_sun_glint_polarized( n_real, n_imag,
                                              wi, wo,
                                              s_theta_i, c_theta_i,
                                              s_theta_o, c_theta_o,
                                              s_phi, c_phi,
                                              s_phi_w, c_phi_w);

        } else {
            value = eval_sun_glint( n_real, n_imag,
                                    s_theta_i, c_theta_i,
                                    s_theta_o, c_theta_o,
                                    s_phi, c_phi,
                                    s_phi_w, c_phi_w );
        }

        return value;
    }

    /**
     * @brief Evaluate the sun glint reflectance as implemented in 6SV.
     *
     * Evaluates the sun glint reflectance at the given index of refraction,
     * incident and outgoing angles, relative azimuthal angle for
     * ingoing and outgoing directions, and relative azimuthal angle for
     * wind direction. The reflectance is computed using the Cox-Munk
     * distribution, the Fresnel equations, and relative tilt of the
     * oceanic surface.
     *
     * @param n_real The real part of the index of refraction of water.
     * @param n_imag The imaginary part of the index of refraction of water.
     * @param s_i The sine of the incident zenith angle.
     * @param c_i The cosine of the incident zenith angle.
     * @param s_o The sine of the outgoing zenith angle.
     * @param c_o The cosine of the outgoing zenith angle.
     * @param s_phi The sine of the incident azimuthal angle.
     * @param c_phi The cosine of the incident azimuthal angle.
     * @param s_phi_w The sine of the relative azimuthal angle of the wind.
     * @param c_phi_w The cosine of the relative azimuthal angle of the wind.
     * @param invert_real Whether to invert the real part of the IOR.
     * @return Float The sun glint reflectance.
     */
    Float eval_sun_glint( ScalarFloat n_real, ScalarFloat n_imag,
                          const Float &s_i, Float c_i, const Float &s_o,
                          Float c_o, const Float &s_phi, const Float &c_phi,
                          const Float &s_phi_w, const Float &c_phi_w,
                          const bool invert_real = false)  const {

        dr::masked(c_i, c_i < 1e-6f) = 1e-6f;
        dr::masked(c_o, c_o < 1e-6f) = 1e-6f;

        // Implementation analog to 6SV
        Float z_x = (-s_o * s_phi) / (c_i + c_o);
        Float z_y = (s_i + s_o * c_phi) / (c_i + c_o);

        // Tilt angle (rad)
        const Float tan_tilt = dr::sqrt(z_x * z_x + z_y * z_y);
        const Float tilt     = dr::atan(tan_tilt);

        // Cox-Munk specular probability
        Float specular_prob = eval_cox_munk(s_phi_w, c_phi_w, z_x, z_y);
        auto mask           = Mask(specular_prob < 0.f);
        specular_prob       = dr::select(mask, 0.f, specular_prob);

        Float cos_2_chi = c_o * c_i + s_o * s_i * c_phi;

        cos_2_chi = dr::select(cos_2_chi > 1.f, 0.999999999f, cos_2_chi);
        cos_2_chi = dr::select(cos_2_chi < -1.f, -0.999999999f, cos_2_chi);

        Float coschi = dr::sqrt(0.5f * (1.f + cos_2_chi));
        Float sinchi = dr::sqrt(0.5f * (1.f - cos_2_chi));

        coschi = dr::select(coschi > 1.f, 0.999999999f, coschi);
        coschi = dr::select(coschi < -1.f, -0.999999999f, coschi);
        sinchi = dr::select(sinchi > 1.f, 0.999999999f, sinchi);
        sinchi = dr::select(sinchi < -1.f, -0.999999999f, sinchi);

        // Invert the real part of the IOR
        if (invert_real) {
            n_real = 1 / n_real;
            n_imag = 0.f;
        }

        Float fresnel_coeff = fresnel_sunglint_legacy<Float>(n_real, n_imag, coschi, sinchi);
        Float num   = dr::Pi<Float> * fresnel_coeff * specular_prob;
        Float denom = 4.f * c_i * c_o * dr::pow(dr::cos(tilt), 4.f);

        Float result = num/denom;

        return result;
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

    /**
     * @brief Evaluate the polarized sun glint reflectance.
     *
     * Evaluates the polarized sun glint reflectance at the given index of
     * refraction, incident and outgoing angles, relative azimuthal angle for
     * ingoing and outgoing directions, and relative azimuthal angle for
     * wind direction. The reflectance is computed using the Cox-Munk
     * distribution, the polarized Fresnel equations as implemented by Mishchenko,
     * and relative tilt of the oceanic surface.
     *
     * @param n_real The real part of the index of refraction of water.
     * @param n_imag The imaginary part of the index of refraction of water.
     * @param wi The incident light direction in graphics convention.
     * @param wo The outgoing light direction in graphics convention.
     * @param s_i The sine of the incident zenith angle.
     * @param c_i The cosine of the incident zenith angle.
     * @param s_o The sine of the outgoing zenith angle.
     * @param c_o The cosine of the outgoing zenith angle.
     * @param s_phi The sine of the incident azimuthal angle.
     * @param c_phi The cosine of the incident azimuthal angle.
     * @param s_phi_w The sine of the relative azimuthal angle of the wind.
     * @param c_phi_w The cosine of the relative azimuthal angle of the wind.
     * @return Float The sun glint reflectance.
     */
    Spectrum eval_sun_glint_polarized( ScalarFloat n_real, ScalarFloat n_imag,
                                       Vector3f wi, Vector3f wo,
                                       const Float &s_i, Float c_i,
                                       const Float &s_o, Float c_o,
                                       const Float &s_phi, const Float &c_phi,
                                       const Float &s_phi_w, const Float &c_phi_w) const {

        Spectrum value = dr::zeros<Spectrum>();

        if constexpr (is_polarized_v<Spectrum>){

            dr::masked(c_i, c_i < 1e-6f) = 1e-6f;
            dr::masked(c_o, c_o < 1e-6f) = 1e-6f;

            // Implementation analog to 6SV
            Float z_x = (-s_o * s_phi) / (c_i + c_o);
            Float z_y = (s_i + s_o * c_phi) / (c_i + c_o);

            // Tilt angle (rad)
            const Float tan_tilt = dr::sqrt(z_x * z_x + z_y * z_y);
            const Float tilt     = dr::atan(tan_tilt);

            // Cox-Munk specular probability
            Float specular_prob = eval_cox_munk(s_phi_w, c_phi_w, z_x, z_y);
            auto mask           = Mask(specular_prob < 0.f);
            specular_prob       = dr::select(mask, 0.f, specular_prob);

            Complex2u n_ext(1.f,0.f);
            Complex2u n_water(n_real,n_imag);
            Spectrum fresnel_coeff = fresnel_sunglint_polarized(n_ext, n_water, -wi, wo);

            Spectrum num = dr::Pi<Float> * fresnel_coeff * specular_prob;
            Float denom  = 4.f * c_i * c_o * dr::pow(dr::cos(tilt), 4.f);
            value = num/denom;
        }

        return value;
    }

    /**
     * @brief Evaluate the underwater light reflectance.
     *
     * Evaluates the underwater light reflectance at the given wavelength,
     * index of refraction, and upwelling and downwelling transmittances.
     * The effects of incident and outgoing direction, wind directions,
     * wind speed, chlorinity, and pigmentation are all precomputed in
     * update or modelled by the transmittances. The reflectance is
     * computed taking the product of the upwelling and downwelling
     * transmittance,attenuted by the ratio of upwell to downwelling
     * irradiance.
     *
     * @param wavelength The wavelength at which to evaluate the reflectance.
     * @param n_real The real part of the index of refraction of water.
     * @param n_imag The imaginary part of the index of refraction of water.
     * @param t_d The downwelling transmittance of radiance inside the water body.
     * @param t_u The upwelling transmittance of radiance inside the water body.
     * @return ScalarFloat The underwater light reflectance.
     */
    Float eval_underlight(const ScalarFloat &wavelength,
                          const ScalarFloat &n_real, const ScalarFloat &n_imag,
                          const Float t_d, const Float t_u) const {

        // Analogue to 6SV, we return 0.0 if the wavelength is outside the range
        // of [0.4, 0.7]
        auto outside_range = Mask(wavelength < 400.f || wavelength > 700.f);

        // Compute the underlight term
        Float underlight = (1.f / (dr::sqr(n_real) + dr::sqr(n_imag))) *
                           (m_r_omega * t_u * t_d) /
                           (1.f - m_underlight_alpha * m_r_omega);

        return dr::select(outside_range, 0.f, underlight);
    }

private:

    /**
     * @brief Evaluate the Cox-Munk distribution.
     *
     * Evaluates the Cox-Munk distribution at the given relative wind direction,
     * x and y components of the sensor/emitter direction, and wind speed. The
     * distribution is computed using the Cox-Munk model.
     *
     * @param s_phi_w The sine of the relative wind direction.
     * @param c_phi_w The cosine of the relative wind direction.
     * @param z_x The x component of the sensor/emitter direction.
     * @param z_y The y component of the sensor/emitter direction.
     * @return Float The probability of the Cox-Munk distribution.
     */
    Float eval_cox_munk(const Float &s_phi_w, const Float &c_phi_w,
                        const Float &z_x, const Float &z_y) const {

        const Float xe = (c_phi_w * z_x + s_phi_w * z_y) / m_sigma_c;
        const Float xn = (-s_phi_w * z_x + c_phi_w * z_y) / m_sigma_u;

        const Float xe2 = xe * xe;
        const Float xn2 = xn * xn;

        Float coef = 1.f - (m_c_21 / 2.f) * (xe2 - 1.f) * xn -
                     (m_c_03 / 6.f) * (xn2 - 3.f) * xn;
        coef = coef + (m_c_40 / 24.f) * (xe2 * xe2 - 6.f * xe2 + 3.f);
        coef = coef + (m_c_04 / 24.f) * (xn2 * xn2 - 6.f * xn2 + 3.f);
        coef = coef + (m_c_22 / 4.f) * (xe2 - 1.f) * (xn2 - 1.f);

        Float prob = coef * dr::InvTwoPi<Float> / (m_sigma_u * m_sigma_c) *
                     dr::exp(-(xe2 + xn2) * 0.5f);
        return prob;
    }


    /**
     * @brief Compute the ratio of upwelling to downwelling irradiance.
     *
     * Computes the ratio of upwelling to downwelling irradiance at the given
     * wavelength and pigmentation. The ratio is computed by performing an
     * iterative computation.
     *
     * @param wavelength The wavelength at which to evaluate the ratio.
     * @param pigmentation The pigmentation of the water.
     * @return ScalarFloat The ratio of upwelling to downwelling irradiance.
     * @note This function is only defined for wavelengths in the range [400,
     * 700] nm.
     */
    ScalarFloat r_omega(const ScalarFloat &wavelength,
                        const ScalarFloat &pigmentation) const {
        ScalarFloat pigment_log = dr::log(pigmentation) / dr::log(10.f);

        // Backscattering coefficient
        ScalarFloat molecular_scatter_coeff =
            m_ocean_props.molecular_scatter_coeff_6s(wavelength);
        ScalarFloat scattering_coeff = 0.30f * dr::pow(pigmentation, 0.62);
        ScalarFloat backscatter_ratio =
            0.002f +
            0.02f * (0.5f - 0.25f * pigment_log) * (550.f / wavelength);
        ScalarFloat backscatter_coeff = 0.5f * molecular_scatter_coeff +
                                        scattering_coeff * backscatter_ratio;

        // (Diffuse) attenuation coefficient
        ScalarFloat k          = m_ocean_props.attn_k(wavelength);
        ScalarFloat chi        = m_ocean_props.attn_chi(wavelength);
        ScalarFloat e          = m_ocean_props.attn_e(wavelength);
        ScalarFloat attn_coeff = k + chi * dr::pow(pigmentation, e);

        // If any of the coefficients is zero, we return zero
        if (backscatter_coeff == 0.f || attn_coeff == 0.f)
            return 0.f;

        // Iterative computation of the reflectance
        ScalarFloat u       = 0.75f;
        ScalarFloat r_omega = 0.33f * backscatter_coeff / u / attn_coeff;

        bool converged = false;
        while (!converged) {
            // Update u
            u = (0.9f * (1.f - r_omega)) / (1.f + 2.25f * r_omega);

            // Update reflectance
            ScalarFloat r_omega_new =
                0.33f * backscatter_coeff / (u * attn_coeff);

            // Create a mask that marks the converged values
            if (dr::abs((r_omega_new - r_omega) / r_omega_new) < 0.0001f) {
                converged = true;
                break;
            }

            // Update reflectance ONLY for non-converged values
            r_omega = r_omega_new;
        }

        return r_omega;
    }

private:
    // Ocean properties
    OceanProperties<Float, Spectrum> m_ocean_props;

    // Whitecap constants
    const ScalarFloat m_f_eff_base     = 0.4f;
    const ScalarFloat m_monahan_alpha  = 2.95e-06f;
    const ScalarFloat m_monahan_lambda = 3.52f;

    // Cox-Munk distribution constants
    const ScalarFloat m_c_40 = 0.40f;
    const ScalarFloat m_c_22 = 0.12f;
    const ScalarFloat m_c_04 = 0.23f;
    // Cox-Munk distribution variables
    ScalarFloat m_sigma_c = 1.f;
    ScalarFloat m_sigma_u = 1.f;
    ScalarFloat m_sigma_iso_sqr = 1.f;
    ScalarFloat m_c_21    = 0.f;
    ScalarFloat m_c_03    = 0.f;
    ScalarFloat m_r_omega = 0.f;

    bool m_shadowing;

    // Underlight parameters
    const ScalarFloat m_underlight_alpha = 0.485f;
};

/**
 * @brief Evaluate the transmittance of the radiance over all
 * provided angles. For each angle, computes the quadrature of
 * the tranmittance.
 * @param utils Instance of ocean utilities templated on FloatP
 * with FloatP vectorized over the size of a packet.
 * @param theta Incident zenith angle.
 * @param theta Incident azimuth angle.
 * @param n_real Real part of the index of refraction.
 * @param n_imag Imaginary part of the index of refraction.
 * @param wind_direction Azimuthal angle of the wind direction.
 * @param upwelling Flag for computing upwelling transmittance.
 * Will compute theta according to snells law and invert the
 * index of refraction.
 */
template <typename Float, typename OceanUtilitiesP>
Float eval_ocean_transmittance(const OceanUtilitiesP &utils, Float theta,
                               Float phi, dr::scalar_t<Float> n_real,
                               dr::scalar_t<Float> n_imag,
                               dr::scalar_t<Float> wind_direction,
                               bool upwelling) {
    MI_IMPORT_CORE_TYPES()

    // number of quadrature points
    int res = 64;

    if (upwelling)
        theta = dr::asin(dr::sin(theta) / n_real);

    using FloatX          = dr::DynamicArray<dr::scalar_t<Float>>;
    auto [nodes, weights] = quad::gauss_legendre<FloatX>(res);
    Float result          = dr::zeros<Float>(dr::width(theta));

    auto [nodes_x, nodes_y]     = dr::meshgrid(nodes, nodes);
    auto [weights_x, weights_y] = dr::meshgrid(weights, weights);

    // nodes_x   -> [0,π/2], nodes_y   -> [0,2π],
    // weights_x -> [0,π/4], weights_y -> [0, π].
    nodes_x = dr::fmadd(nodes_x, 0.25f * dr::Pi<FloatX>, 0.25f * dr::Pi<FloatX>);
    nodes_y = dr::fmadd(nodes_y, dr::Pi<FloatX>, dr::Pi<FloatX>);
    weights_x                     = 0.25f * dr::Pi<FloatX> * weights_x;
    weights_y                     = dr::Pi<FloatX> * weights_y;
    auto [s_zeniths, c_zeniths]   = dr::sincos(nodes_x);
    auto [s_azimuths, c_azimuths] = dr::sincos(nodes_y);

    auto [s_the, c_the]     = dr::sincos(theta);
    auto [s_phi_w, c_phi_w] = dr::sincos(phi - wind_direction);

    using FloatP        = dr::Packet<dr::scalar_t<Float>>;
    size_t packet_count = dr::width(theta) / FloatP::Size;
    Assert(dr::width(theta) % FloatP::Size == 0);

    // cycle through each packet.
    for (size_t i = 0; i < packet_count; ++i) {
        FloatP s_the_p, c_the_p, s_phi_w_p, c_phi_w_p;
        s_the_p   = dr::load<FloatP>(s_the.data() + i * FloatP::Size);
        c_the_p   = dr::load<FloatP>(c_the.data() + i * FloatP::Size);
        s_phi_w_p = dr::load<FloatP>(s_phi_w.data() + i * FloatP::Size);
        c_phi_w_p = dr::load<FloatP>(c_phi_w.data() + i * FloatP::Size);

        FloatP result_p = 0.f;
        FloatP td = 0.f, summ = 0.f;

        // compute the quadrature for each packet.
        for (size_t j = 0; j < dr::width(nodes_x); ++j) {

            ScalarFloat s_zenith  = s_zeniths[j];
            ScalarFloat c_zenith  = c_zeniths[j];
            ScalarFloat s_azimuth = s_azimuths[j];
            ScalarFloat c_azimuth = c_azimuths[j];

            ScalarFloat weight_x = weights_x[j];
            ScalarFloat weight_y = weights_y[j];

            ScalarFloat geometry       = c_zenith * s_zenith;
            ScalarFloat geometryWeight = geometry * weight_y * weight_x;

            FloatP glint = utils.eval_sun_glint(
                n_real, n_imag, s_the_p, c_the_p, s_zenith, c_zenith, s_azimuth,
                c_azimuth, s_phi_w_p, c_phi_w_p, upwelling);
            td += glint * geometryWeight;
            summ += geometryWeight;
        }

        dr::masked(td, td >= summ) = summ;
        result_p = 1.f - (td / summ);

        dr::store(result.data() + i * FloatP::Size, result_p);
    }

    return result;
}


template <typename Float, typename Spectrum>
class OceanBSDF final : public BSDF<Float, Spectrum> {
public:
    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture)

    /**
     * @brief Construct a new OceanBSDF object.
     *
     * @param props A set of properties to initialize the oceanic BSDF.
     */
    OceanBSDF(const Properties &props) : Base(props) {
        // Retrieve parameters
        m_wavelength     = props.get<ScalarFloat>("wavelength");
        m_wind_speed     = props.get<ScalarFloat>("wind_speed", 0.1f);
        m_wind_direction = props.get<ScalarFloat>("wind_direction", 0.);
        m_chlorinity     = props.get<ScalarFloat>("chlorinity", 19.f);
        m_pigmentation   = props.get<ScalarFloat>("pigmentation", 0.3f);
        m_shininess      = props.get<ScalarFloat>("shininess", 50.f);
        m_component      = props.get<ScalarInt32>("component", 0);
        m_shadowing      = props.get<bool>("shadowing", true);
        m_accel          = props.get<bool>("accel", true);

        // convert from North Left to East Right.
        m_wind_direction = -m_wind_direction + 90.f;
        // m_wind_direction % 360.
        m_wind_direction = m_wind_direction - (360.f * floor(m_wind_direction / 360.f));
        // Degree to radians.
        m_wind_direction = dr::deg_to_rad(m_wind_direction);

        update();

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
        dr::set_attr(this, "flags", m_flags);
    }

    void traverse(TraversalCallback *callback) override {
        callback->put_parameter("wavelength", m_wavelength,
                                +ParamFlags::NonDifferentiable);
        callback->put_parameter("wind_speed", m_wind_speed,
                                +ParamFlags::Differentiable);
        callback->put_parameter("wind_direction", m_wind_direction,
                                +ParamFlags::Differentiable);
        callback->put_parameter("chlorinity", m_chlorinity,
                                +ParamFlags::Differentiable);
        callback->put_parameter("pigmentation", m_pigmentation,
                                +ParamFlags::Differentiable);
        callback->put_parameter("shininess", m_shininess,
                                +ParamFlags::NonDifferentiable);
        callback->put_parameter("coverage", m_coverage,
                                +ParamFlags::NonDifferentiable);
    }

    void parameters_changed(const std::vector<std::string> &/*keys*/) override {
        update();
    }


    /**
     * @brief Update the variables that can be preprocessed.
     * This includes the complex index of refraction, upwelling
     * and downwelling transmittance textures and the ocean utils
     * variables (cox-munk, underlight).
     */
    void update() {

        // compute the index of refraction
        std::tie(m_n_real, m_n_imag) =
            m_ocean_utils.eval_index_refraction(m_wavelength, m_chlorinity);

        {
            // Pre-compute textures for the upwelling and downwelling
            // transmittances of radiance in the water body.
            using FloatX = DynamicBuffer<ScalarFloat>;
            using FloatP = dr::Packet<dr::scalar_t<Float>>;

            FloatX zeniths = dr::maximum(
                0.f, dr::linspace<FloatX>(0.f, dr::Pi<ScalarFloat> * 0.5f,
                                         MI_OCEAN_TRANSMITTANCE_RES));
            FloatX azimuths = dr::maximum(
                0.f, dr::linspace<FloatX>(0.f, 2.f * dr::Pi<ScalarFloat>,
                                         MI_OCEAN_TRANSMITTANCE_RES));
            auto [zeniths_x, azimuths_y] = dr::meshgrid(zeniths, azimuths);

            // Need to initialize a OceanUtilies which can work with
            // packet-sized float vectors i.e. FloatP. Not the best approach,
            // but the simplest that avoids a complete refactor of the plugin.
            OceanUtilities<FloatP, Spectrum> ocean_utils;
            ocean_utils.update(m_wavelength, m_wind_speed, m_pigmentation, m_shadowing);

            FloatX downwelling = eval_ocean_transmittance(
                ocean_utils, zeniths_x, azimuths_y, m_n_real, m_n_imag,
                m_wind_direction, false);
            FloatX upwelling = eval_ocean_transmittance(
                ocean_utils, zeniths_x, azimuths_y, m_n_real, m_n_imag,
                m_wind_direction, true);

            size_t shape[3] = { MI_OCEAN_TRANSMITTANCE_RES,
                                MI_OCEAN_TRANSMITTANCE_RES, 1 };

            m_downwelling_transmittance =
                Texture2f(TensorXf(downwelling.data(), 3, shape), m_accel,
                        m_accel, dr::FilterMode::Linear, dr::WrapMode::Clamp);
            m_upwelling_transmittance =
                Texture2f(TensorXf(upwelling.data(), 3, shape), m_accel,
                        m_accel, dr::FilterMode::Linear, dr::WrapMode::Clamp);

            dr::eval(m_downwelling_transmittance);
            dr::eval(m_upwelling_transmittance);
        }

        // Update other useful values
        m_ocean_utils.update(m_wavelength, m_wind_speed, m_pigmentation, m_shadowing);

        // Pre-compute whitecap coverage
        m_coverage = m_ocean_utils.eval_whitecap_coverage(m_wind_speed);
    }

    /**
     * @brief Evaluate the whitecap reflectance.
     *
     * Evaluates the whitecap reflectance at the provided wavelength and wind
     * speed.
     *
     * @return Float The whitecap reflectance.
     */
    Float eval_whitecaps() const {
        return m_ocean_utils.eval_whitecaps(m_wavelength, m_wind_speed);
    }

    /**
     * @brief Evaluate the sun glint reflectance.
     *
     * Evaluates the sun glint reflectance at the provided wavelength, incident
     * and outgoing directions, wind direction, wind speed, and chlorinity.
     *
     * @param wi The incident direction of the light in graphics convention.
     * @param wo The outgoing direction of the light in graphics convention.
     * @return Float The sun glint reflectance.
     */
    Spectrum eval_glint(const Vector3f &wi, const Vector3f &wo) const {
        return m_ocean_utils.eval_glint(m_n_real, m_n_imag, wi, wo,
                                        m_wind_direction);
    }

    /**
     * @brief Evaluate the underwater light reflectance.
     *
     * Evaluates the underwater light reflectance at the provided wavelength,
     * incident and outgoing directions, wind direction, wind speed, chlorinity,
     * and pigmentation.
     *
     * @param wi The incident direction of the light.
     * @param wo The outgoing direction of the light.
     * @return Float The underwater light reflectance.
     */
    Float eval_underlight(const Vector3f &wi, const Vector3f &wo) const {
        Vector2f uv;
        Float t_d, t_u;

        uv.x() = (dr::acos(wi.z()) * (dr::InvPi<Float> * 2.f));
        uv.y() = dr::atan2(wi.y(), wi.x()) * dr::InvTwoPi<Float>;
        uv.y() = uv.y() -
                 (1.f * dr::floor(uv.y() / 1.f)); // equivalent to uv.y % 1.f

        m_downwelling_transmittance.eval(uv, &t_d);

        uv.x() = (dr::acos(wo.z()) * (dr::InvPi<Float> * 2.));
        uv.y() = dr::atan2(wi.y(), wi.x()) * dr::InvTwoPi<Float>;
        uv.y() = uv.y() -
                 (1.f * dr::floor(uv.y() / 1.f)); // equivalent to uv.y % 1.f

        m_upwelling_transmittance.eval(uv, &t_u);

        return m_ocean_utils.eval_underlight(m_wavelength, m_n_real, m_n_imag,
                                             t_d, t_u);
    }

    std::pair<BSDFSample3f, Spectrum>
    sample(const BSDFContext &ctx, const SurfaceInteraction3f &si,
           Float sample1, const Point2f &sample2, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);

        bool has_diffuse  = ctx.is_enabled(BSDFFlags::DiffuseReflection, 0),
             has_specular = ctx.is_enabled(BSDFFlags::GlossyReflection, 1);

        Float cos_theta_i = Frame3f::cos_theta(si.wi);
        active &= cos_theta_i > 0.f;

        BSDFSample3f bs   = dr::zeros<BSDFSample3f>();
        Spectrum result(0.f);

        if (unlikely(dr::none_or<false>(active)) || (!has_diffuse && !has_specular))
            return { bs, 0.f };

        Float prob_diff = m_coverage;

        Mask sample_diffuse = active & (sample1 < prob_diff),
             sample_glint   = active && !sample_diffuse;

        if (dr::any_or<true>(sample_diffuse)) {
            // In the case of sampling the diffuse component, the outgoing
            // direction is sampled from a cosine-weighted hemisphere.
            dr::masked(bs.wo, sample_diffuse) =
                warp::square_to_cosine_hemisphere(sample2);
            dr::masked(bs.sampled_component, sample_diffuse) = 0;
            dr::masked(bs.sampled_type, sample_diffuse) =
                +BSDFFlags::DiffuseReflection;
        }

        if (dr::any_or<true>(sample_glint)) {
            if (m_shininess >= 0.f) {
                // For Blinn-Phong, we need to sample the half-vector
                Float ksi_1 = sample2.x(), ksi_2 = sample2.y();

                Float phi_h   = dr::TwoPi<Float> * ksi_1;
                Float theta_h = dr::acos(dr::pow(ksi_2, 1.f / (m_shininess + 2.f)));

                Vector3f half = dr::normalize(
                    Vector3f(dr::sin(theta_h) * dr::cos(phi_h),
                             dr::sin(theta_h) * dr::sin(phi_h),
                             dr::cos(theta_h)));

                Vector3f wo = 2.f * dr::dot(si.wi, half) * half - si.wi;

                // In the case of sampling the glint component, the outgoing
                // direction is sampled using the Blinn-Phong distribution.
                dr::masked(bs.wo, sample_glint)                = wo;
                dr::masked(bs.sampled_component, sample_glint) = 1;
                dr::masked(bs.sampled_type, sample_glint) =
                    +BSDFFlags::GlossyReflection;

            } else {
                dr::masked(bs.wo, sample_glint) =
                    warp::square_to_cosine_hemisphere(sample2);
                dr::masked(bs.sampled_component, sample_glint) = 0;
                dr::masked(bs.sampled_type, sample_glint) =
                    +BSDFFlags::GlossyReflection;
            }
        }

        bs.pdf = pdf(ctx, si, bs.wo, active);
        bs.eta = 1.f;

        active &= bs.pdf > 0.f;
        result = eval(ctx, si, bs.wo, active);

        return { bs, result/bs.pdf & (active && bs.pdf > 0.f) };
    }

    Spectrum eval(const BSDFContext &ctx, const SurfaceInteraction3f &si,
                  const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        bool has_diffuse = ctx.is_enabled(BSDFFlags::DiffuseReflection, 0),
             has_glint   = ctx.is_enabled(BSDFFlags::GlossyReflection, 1);

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        // Ensure incoming and outgoing directions are in the upper hemisphere
        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

        if (unlikely(dr::none_or<false>(active) || !has_glint && !has_diffuse))
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

        if (has_diffuse) {
            whitecap_reflectance   = eval_whitecaps();
            underlight_reflectance = eval_underlight(wo_hat, wi_hat);
            // Diffuse scattering implies no polarization.
            result += depolarizer<Spectrum>(
                whitecap_reflectance + (1.f - whitecap_reflectance) * underlight_reflectance
            );
        }

        if (has_glint) {

            if constexpr( is_polarized_v<Spectrum> ) {
                // If sun glint is enabled, compute the glint reflectance
                glint_reflectance = eval_glint(wo_hat, wi_hat);

                /* The Stokes reference frame vector of this matrix lies in the meridian plane
                spanned by wi and n. */
                Vector3f n(0.f, 0.f, 1.f);
                Vector3f p_axis_in  =
                    dr::normalize(dr::cross(dr::normalize(dr::cross(n,-wo_hat)),-wo_hat));
                Vector3f p_axis_out =
                    dr::normalize(dr::cross(dr::normalize(dr::cross(n,wi_hat)),wi_hat));

                dr::masked( p_axis_in, dr::any(dr::isnan(p_axis_in)) ) =
                    Vector3f(0.f, 1.f, 0.f);
                dr::masked( p_axis_out, dr::any(dr::isnan(p_axis_out)) ) =
                    Vector3f(0.f, 1.f, 0.f);

                // Rotate in/out reference vector of `value` s.t. it aligns with the implicit
                // Stokes bases of -wo_hat & wi_hat. */
                glint_reflectance = mueller::rotate_mueller_basis(glint_reflectance,
                                        -wo_hat, p_axis_in, mueller::stokes_basis(-wo_hat),
                                        wi_hat, p_axis_out, mueller::stokes_basis(wi_hat));

            } else {
                // If sun glint is enabled, compute the glint reflectance
                glint_reflectance = eval_glint(wo_hat, wi_hat);
            }

            result += (1.f - m_coverage) * glint_reflectance;
        }

        dr::masked(result, active) *= cos_theta_o * dr::InvPi<Float>;

        // For debugging purposes, channel indicates which BRDF term to evaluate
        switch (m_component) {
            case 1:
                result[active] = depolarizer<Spectrum>(whitecap_reflectance);
                break;
            case 2:
                result[active] =
                    (1.f - m_coverage) * glint_reflectance;
                break;
            case 3:
                result[active] =
                    depolarizer<Spectrum>((1.f - whitecap_reflectance) * underlight_reflectance);
                break;
            case 4:
                result[active] = depolarizer<Spectrum>(
                    whitecap_reflectance +
                    (1.f - whitecap_reflectance) * underlight_reflectance);
                break;
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

        if (unlikely((!has_diffuse && !has_specular) ||
                     dr::none_or<false>(active)))
            return 0.f;

        Float weight_diffuse = m_coverage, weight_specular = (1.f - m_coverage);

        // Check if the normal has only zeros. If this is the case, use a
        // default normal
        Vector3f normal   = si.n;
        Mask degen_normal = dr::all(dr::eq(normal, Vector3f(0.f)));
        dr::masked(normal, degen_normal) = Vector3f(0.f, 0.f, 1.f);

        Vector3f half    = dr::normalize(si.wi + wo);
        Float projection = dr::dot(half, normal);
        Float D          = (m_shininess + 2.f) / dr::TwoPi<Float> *
                           dr::pow(projection, m_shininess);

        // We multiply the probability of the specular lobe with the pdf of
        // the Blinn-Phong distribution and the probability of the diffuse lobe
        // with the pdf of the cosine-weighted hemisphere.
        Float pdf_diffuse  = warp::square_to_cosine_hemisphere_pdf(wo);
        Float pdf_specular = (m_shininess >= 0.f) ?
                (D * projection) / (4.f * dr::dot(si.wi, half)) :
                warp::square_to_cosine_hemisphere_pdf(wo);

        Float pdf =
            weight_diffuse * pdf_diffuse + weight_specular * pdf_specular;

        // If the outgoing direction is in the lower hemisphere, we return zero
        Float cos_theta_o = Frame3f::cos_theta(wo);

        return dr::select(cos_theta_o > 0.f, pdf, 0.f);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "OceanLegacy[" << std::endl
            << "  component = " << string::indent(m_component) << "," << std::endl
            << "  wavelength = " << string::indent(m_wavelength) << "," << std::endl
            << "  wind_speed = " << string::indent(m_wind_speed) << "," << std::endl
            << "  wind_direction = " << string::indent(m_wind_direction) << "," << std::endl
            << "  chlorinity = " << string::indent(m_chlorinity) << "," << std::endl
            << "  pigmentation = " << string::indent(m_pigmentation) << "," << std::endl
            << "  shininess = " << string::indent(m_shininess) << "," << std::endl
            << "  coverage = " << string::indent(m_coverage) << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS()
private:
    //  User-provided fields
    ScalarInt32 m_component;
    ScalarFloat m_wavelength;
    ScalarFloat m_wind_speed;
    ScalarFloat m_wind_direction;
    ScalarFloat m_chlorinity;
    ScalarFloat m_pigmentation;
    ScalarFloat m_shininess;
    ScalarFloat m_coverage;
    bool m_shadowing;
    bool m_accel;

    // On update fields
    ScalarFloat m_n_real;
    ScalarFloat m_n_imag;

    OceanUtilities<Float, Spectrum> m_ocean_utils;

    Texture2f m_downwelling_transmittance;
    Texture2f m_upwelling_transmittance;
};

MI_IMPLEMENT_CLASS_VARIANT(OceanBSDF, BSDF)
MI_EXPORT_PLUGIN(OceanBSDF, "Ocean material")
NAMESPACE_END(mitsuba)
