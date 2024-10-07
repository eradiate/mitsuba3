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
#include <tuple>

// Transmittance angular resolution
// i.e. texture resolution
#define MI_OCEAN_TRANSMITTANCE_RES 64

NAMESPACE_BEGIN(mitsuba)

/**!

.. _plugin-bsdf-ocean_legacy:

(Legacy 6S) Oceanic reflection model (:monosp:`ocean-legacy`)
---------------------------------------------------------------

.. pluginparameters::

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
     (Default: :monosp:`50.`).

 * - component
   - |int|
   - Debug: specifies which component of the oceanic reflection model to evaluate.
     Default: 0 Component 0 is used to evaluate the total oceanic reflectance.
     Component 1 evaluates the whitecap reflectance. Component 2 evaluates the sun
     glint reflectance. Component 3 evaluates the underlight reflectance.

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

// Header content - THIS CAN BE MOVED IF NECESSARY
#ifndef OCEAN_PROPS
#define OCEAN_PROPS

template <typename Float, typename Spectrum> class OceanProperties {
public:
    MI_IMPORT_TYPES()

    using FloatX = DynamicBuffer<Float>;
    using Value  = dr::Array<ScalarFloat, 1>;

    /**
     * @brief Construct a new Ocean Properties object and initializes the data.
     *
     * Initializes the data for the effective reflectance of whitecaps, the
     * complex index of refraction of water, the water scattering and
     * attenuation coefficients, and the molecular scattering coefficients. The
     * data is taken from various sources in the literature.
     */
    OceanProperties() {
        /*
        * Effective reflectance of whitecaps (Whitlock et al. 1982)
        * Wavelength specified in um.
        * Wavelengths are stored in a regular grid, using range instead.
        constexpr ScalarFloat wc_wavelengths[] = {  
            0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1,
            1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 
            2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1,
            3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0 };
        */
        constexpr ScalarFloat wc_data[] = {
            0.220, 0.220, 0.220, 0.220, 0.220, 0.220, 0.215, 0.210,
            0.200, 0.190, 0.175, 0.155, 0.130, 0.080, 0.100, 0.105,
            0.100, 0.080, 0.045, 0.055, 0.065, 0.060, 0.055, 0.040,
            0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
            0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000
        };

        // Complex index of refraction of water (Hale & Querry 1973)
        // !! Note that the wavelength is in nm to align with the rest of
        // Mitsuba !!
        constexpr ScalarFloat ior_wavelengths[] = {
            200,  225,  250,  275,  300,  325,  345,  375,  400,  425,  445,
            475,  500,  525,  550,  575,  600,  625,  650,  675,  700,  725,
            750,  775,  800,  825,  850,  875,  900,  925,  950,  975,  1000,
            1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2650, 2700, 2750,
            2800, 2850, 2900, 2950, 3000, 3050, 3100, 3150, 3200, 3250, 3300,
            3350, 3400, 3450, 3500, 3600, 3700, 3800, 3900, 4000
        };
        constexpr ScalarFloat ior_real_data[] = {
            1.369, 1.373, 1.362, 1.354, 1.349, 1.346, 1.343, 1.341,
            1.339, 1.338, 1.337, 1.336, 1.335, 1.334, 1.333, 1.333,
            1.332, 1.332, 1.331, 1.331, 1.331, 1.330, 1.330, 1.330,
            1.329, 1.329, 1.329, 1.328, 1.328, 1.328, 1.327, 1.327,
            1.327, 1.324, 1.321, 1.317, 1.312, 1.306, 1.296, 1.279,
            1.242, 1.219, 1.188, 1.157, 1.142, 1.149, 1.201, 1.292,
            1.371, 1.426, 1.467, 1.483, 1.478, 1.467, 1.450, 1.432,
            1.420, 1.410, 1.400, 1.385, 1.374, 1.364, 1.357, 1.351
        };
        constexpr ScalarFloat ior_cplx_data[] = {
            1.10e-07, 4.90e-08, 3.35e-08, 2.35e-08, 1.60e-08, 1.08e-08,
            6.50e-09, 3.50e-09, 1.86e-09, 1.30e-09, 1.02e-09, 9.35e-10,
            1.00e-09, 1.32e-09, 1.96e-09, 3.60e-09, 1.09e-08, 1.39e-08,
            1.64e-08, 2.23e-08, 3.35e-08, 9.15e-08, 1.56e-07, 1.48e-07,
            1.25e-07, 1.82e-07, 2.93e-07, 3.91e-07, 4.86e-07, 1.06e-06,
            2.93e-06, 3.48e-06, 2.89e-06, 9.89e-06, 1.38e-04, 8.55e-05,
            1.15e-04, 1.10e-03, 2.89e-04, 9.56e-04, 3.17e-03, 6.70e-03,
            1.90e-02, 5.90e-02, 1.15e-01, 1.85e-01, 2.68e-01, 2.98e-01,
            2.72e-01, 2.40e-01, 1.92e-01, 1.35e-01, 9.24e-02, 6.10e-02,
            3.68e-02, 2.61e-02, 1.95e-02, 1.32e-02, 9.40e-03, 5.15e-03,
            3.60e-03, 3.40e-03, 3.80e-03, 4.60e-03
        };

        /*
        * Water scattering and attenuation coefficient data (Morel 1988)
        * Wavelength specified in nm.
        * Wavelengths are stored in a regular grid, using range instead.
        constexpr ScalarFloat attn_wavelengths[] = { 
            400, 405, 410, 415, 420, 425, 430, 435, 440, 445, 
            450, 455, 460, 465, 470, 475, 480, 485, 490, 495, 
            500, 505, 510, 515, 520, 525, 530, 535, 540, 545, 
            550, 555, 560, 565, 570, 575, 580, 585, 590, 595, 
            600, 605, 610, 615, 620, 625, 630, 635, 640, 645, 
            650, 655, 660, 665, 670, 675, 680, 685, 690, 695, 
            700 
        };
        */
        constexpr ScalarFloat attn_k[] = {
            0.0209, 0.0200, 0.0196, 0.0189, 0.0183, 0.0182, 0.0171, 0.0170,
            0.0168, 0.0166, 0.0168, 0.0170, 0.0173, 0.0174, 0.0175, 0.0184,
            0.0194, 0.0203, 0.0217, 0.0240, 0.0271, 0.0320, 0.0384, 0.0445,
            0.0490, 0.0505, 0.0518, 0.0543, 0.0568, 0.0615, 0.0640, 0.0640,
            0.0717, 0.0762, 0.0807, 0.0940, 0.1070, 0.1280, 0.1570, 0.2000,
            0.2530, 0.2790, 0.2960, 0.3030, 0.3100, 0.3150, 0.3200, 0.3250,
            0.3300, 0.3400, 0.3500, 0.3700, 0.4050, 0.4180, 0.4300, 0.4400,
            0.4500, 0.4700, 0.5000, 0.5500, 0.6500
        };
        constexpr ScalarFloat attn_chi[] = {
            0.1100, 0.1110, 0.1125, 0.1135, 0.1126, 0.1104, 0.1078, 0.1065,
            0.1041, 0.0996, 0.0971, 0.0939, 0.0896, 0.0859, 0.0823, 0.0788,
            0.0746, 0.0726, 0.0690, 0.0660, 0.0636, 0.0600, 0.0578, 0.0540,
            0.0498, 0.0475, 0.0467, 0.0450, 0.0440, 0.0426, 0.0410, 0.0400,
            0.0390, 0.0375, 0.0360, 0.0340, 0.0330, 0.0328, 0.0325, 0.0330,
            0.0340, 0.0350, 0.0360, 0.0375, 0.0385, 0.0400, 0.0420, 0.0430,
            0.0440, 0.0445, 0.0450, 0.0460, 0.0475, 0.0490, 0.0515, 0.0520,
            0.0505, 0.0440, 0.0390, 0.0340, 0.0300
        };
        constexpr ScalarFloat attn_e[] = {
            0.668, 0.672, 0.680, 0.687, 0.693, 0.701, 0.707, 0.708, 0.707,
            0.704, 0.701, 0.699, 0.700, 0.703, 0.703, 0.703, 0.703, 0.704,
            0.702, 0.700, 0.700, 0.695, 0.690, 0.685, 0.680, 0.675, 0.670,
            0.665, 0.660, 0.655, 0.650, 0.645, 0.640, 0.630, 0.623, 0.615,
            0.610, 0.614, 0.618, 0.622, 0.626, 0.630, 0.634, 0.638, 0.642,
            0.647, 0.653, 0.658, 0.663, 0.667, 0.672, 0.677, 0.682, 0.687,
            0.695, 0.697, 0.693, 0.665, 0.640, 0.620, 0.600
        };
        // IMPORTANT: This table uses the values provided by Morel, which are
        // different than the ones from 6SV
        constexpr ScalarFloat molecular_scatter_coeffs[] = {
            0.00618095, 0.00578095, 0.00547619, 0.00517619, 0.00492222,
            0.0046746,  0.00447143, 0.00426825, 0.00406508, 0.0038619,
            0.00365873, 0.00346667, 0.00331429, 0.0031619,  0.00300952,
            0.00287143, 0.00276984, 0.00265238, 0.0025,     0.00236508,
            0.00226349, 0.0021619,  0.00206032, 0.00195873, 0.00185714,
            0.00177778, 0.00172698, 0.00167619, 0.0016254,  0.0015746,
            0.00152381, 0.00144603, 0.00134444, 0.0013,     0.0013,
            0.00126984, 0.00121905, 0.00116825, 0.00111746, 0.00107,
            0.00102429, 0.00098556, 0.00095,    0.0009181,  0.00088762,
            0.00085714, 0.00082667, 0.00079619, 0.00076571, 0.00073937,
            0.00071397, 0.00069286, 0.00067254, 0.00065222, 0.0006319,
            0.00061159, 0.00059127, 0.00057095, 0.00055063, 0.00053524,
            0.00052
        };

        constexpr ScalarFloat molecular_scatter_coeffs_6s[] = {
            0.0076, 0.0072, 0.0068, 0.0064, 0.0061, 0.0058, 0.0055, 0.0052,
            0.0049, 0.0047, 0.0045, 0.0043, 0.0041, 0.0039, 0.0037, 0.0036,
            0.0034, 0.0033, 0.0031, 0.0030, 0.0029, 0.0027, 0.0026, 0.0025,
            0.0024, 0.0023, 0.0022, 0.0022, 0.0021, 0.0020, 0.0019, 0.0018,
            0.0018, 0.0017, 0.0017, 0.0016, 0.0016, 0.0015, 0.0015, 0.0014,
            0.0014, 0.0013, 0.0013, 0.0012, 0.0012, 0.0011, 0.0011, 0.0010,
            0.0010, 0.0010, 0.0010, 0.0009, 0.0008, 0.0008, 0.0008, 0.0007,
            0.0007, 0.0007, 0.0007, 0.0007, 0.0007
        };

        // Construct distributions from the provided data sets
        // !! Note that unlike the wavelength is in nm to align with the rest of
        // Mitsuba !!
        m_effective_reflectance = ContinuousDistribution<ScalarFloat>(
            ScalarVector2f(200.f, 4000.f), wc_data, std::size(wc_data));

        m_ior_real = IrregularContinuousDistribution<ScalarFloat>(
            ior_wavelengths, ior_real_data, std::size(ior_real_data));

        m_ior_imag = IrregularContinuousDistribution<ScalarFloat>(
            ior_wavelengths, ior_cplx_data, std::size(ior_cplx_data));

        m_attn_k = ContinuousDistribution<ScalarFloat>(
            ScalarVector2f(400.f, 700.f), attn_k, std::size(attn_k));

        m_attn_chi = ContinuousDistribution<ScalarFloat>(
            ScalarVector2f(400.f, 700.f), attn_chi, std::size(attn_chi));

        m_attn_e = ContinuousDistribution<ScalarFloat>(
            ScalarVector2f(400.f, 700.f), attn_e, std::size(attn_e));

        m_molecular_scatter_coeffs = ContinuousDistribution<ScalarFloat>(
            ScalarVector2f(400.f, 700.f), molecular_scatter_coeffs,
            std::size(molecular_scatter_coeffs));

        m_molecular_scatter_coeffs_6s = ContinuousDistribution<ScalarFloat>(
            ScalarVector2f(400.f, 700.f), molecular_scatter_coeffs_6s,
            std::size(molecular_scatter_coeffs_6s));
    }

    /**
     * @brief Evaluate the effective reflectance of whitecaps.
     *
     * Evaluates the effective reflectance of whitecaps at the given
     * wavelength. The value returned already takes into account
     * the base offset of 0.4 as described by (Koepke 1984).
     *
     * @param wavelength The wavelength at which to evaluate the reflectance.
     * @return ScalarFloat The effective reflectance of whitecaps.
     */
    ScalarFloat effective_reflectance(const ScalarFloat &wavelength) const {
        return m_effective_reflectance.eval_pdf(wavelength);
    }

    /**
     * @brief Evaluate the real index of refraction of water.
     *
     * Evaluates the real index of refraction of water at the given
     * wavelength. The value returned is the real part of the complex
     * index of refraction as described by (Hale & Querry 1973).
     *
     * @param wavelength The wavelength at which to evaluate the index of
     * refraction.
     * @return ScalarFloat The real part of the index of refraction.
     */
    ScalarFloat ior_real(const ScalarFloat &wavelength) const {
        return m_ior_real.eval_pdf(wavelength);
    }

    /**
     * @brief Evaluate the complex index of refraction of water.
     *
     * Evaluates the complex index of refraction of water at the given
     * wavelength. The value returned is the imaginary part of the complex
     * index of refraction as described by (Hale & Querry 1973).
     *
     * @param wavelength The wavelength at which to evaluate the index of
     * refraction.
     * @return ScalarFloat The imaginary part of the index of refraction.
     */
    ScalarFloat ior_cplx(const ScalarFloat &wavelength) const {
        return m_ior_imag.eval_pdf(wavelength);
    }

    /**
     * @brief Evaluate the K-term of the attenuation coefficient of water.
     *
     * Evaluates the K-term of the attenuation coefficient of water at the given
     * wavelength. The value returned is the K-term as described by (Morel
     * 1988).
     *
     * @param wavelength The wavelength at which to evaluate the K-term of the
     * attenuation coefficient.
     * @return ScalarFloat The K-term of the attenuation coefficient.
     */
    ScalarFloat attn_k(const ScalarFloat &wavelength) const {
        return m_attn_k.eval_pdf(wavelength);
    }

    /**
     * @brief Evaluate the Chi-term of the attenuation coefficient of water.
     *
     * Evaluates the Chi-term of the attenuation coefficient of water at the
     * given wavelength. The value returned is the Chi-term as described by
     * (Morel 1988).
     *
     * @param wavelength The wavelength at which to evaluate the Chi-term of the
     * attenuation coefficient.
     * @return ScalarFloat The Chi-term of the attenuation coefficient.
     */
    ScalarFloat attn_chi(const ScalarFloat &wavelength) const {
        return m_attn_chi.eval_pdf(wavelength);
    }

    /**
     * @brief Evaluate the E-term of the attenuation coefficient of water.
     *
     * Evaluates the E-term of the attenuation coefficient of water at the given
     * wavelength. The value returned is the E-term as described by (Morel
     * 1988).
     *
     * @param wavelength The wavelength at which to evaluate the E-term of the
     * attenuation coefficient.
     * @return ScalarFloat The E-term of the attenuation coefficient.
     */
    ScalarFloat attn_e(const ScalarFloat &wavelength) const {
        return m_attn_e.eval_pdf(wavelength);
    }

    /**
     * @brief Evaluate the molecular scattering coefficient of water.
     *
     * Evaluates the molecular scattering coefficient of water at the given
     * wavelength. The value returned is the molecular scattering coefficient
     * as described by (Morel 1988).
     *
     * @param wavelength The wavelength at which to evaluate the molecular
     * scattering coefficient.
     * @return ScalarFloat The molecular scattering coefficient.
     */
    ScalarFloat molecular_scatter_coeff(const ScalarFloat &wavelength) const {
        return m_molecular_scatter_coeffs.eval_pdf(wavelength);
    }

    /**
     * @brief Evaluate the molecular scattering coefficient of water.
     *
     * Evaluates the molecular scattering coefficient of water at the given
     * wavelength. The value returned is the molecular scattering coefficient
     * as described by the 6S radiative transfer model.
     *
     * @param wavelength The wavelength at which to evaluate the molecular
     * scattering coefficient.
     * @return ScalarFloat The molecular scattering coefficient.
     */
    ScalarFloat
    molecular_scatter_coeff_6s(const ScalarFloat &wavelength) const {
        return m_molecular_scatter_coeffs_6s.eval_pdf(wavelength);
    }

private:
    // Effective reflectance of whitecaps
    ContinuousDistribution<ScalarFloat> m_effective_reflectance;

    // Real/Complex IOR of water (Hale & Querry 1973)
    IrregularContinuousDistribution<ScalarFloat> m_ior_real;
    IrregularContinuousDistribution<ScalarFloat> m_ior_imag;

    // Water scattering and attenuation coefficients (Morel 1988)
    ContinuousDistribution<ScalarFloat> m_attn_k;
    ContinuousDistribution<ScalarFloat> m_attn_chi;
    ContinuousDistribution<ScalarFloat> m_attn_e;
    ContinuousDistribution<ScalarFloat> m_molecular_scatter_coeffs;
    ContinuousDistribution<ScalarFloat> m_molecular_scatter_coeffs_6s;
};

template <typename Float, typename Spectrum> class OceanUtilities {
public:
    MI_IMPORT_TYPES()

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
                const ScalarFloat &pigmentation) {
        // Update Cox-Munk variables
        m_sigma_c = dr::sqrt(0.003f + 0.00192f * wind_speed);
        m_sigma_u = dr::sqrt(0.00316f * wind_speed);
        m_c_21    = 0.01f - 0.0086f * wind_speed;
        m_c_03    = 0.04f - 0.033f * wind_speed;

        m_r_omega = r_omega(wavelength, pigmentation);
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
        ScalarFloat n_real =
            m_ocean_props.ior_real(wavelength) + friedman_sverdrup(chlorinity);
        ScalarFloat n_imag = m_ocean_props.ior_cplx(wavelength);

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
     * @param wi The incident direction of the light.
     * @param wo The outgoing direction of the light.
     * @param wind_direction The direction of the wind relative to the incident
     * light.
     * @return Float The sun glint reflectance.
     */
    Float eval_glint(const ScalarFloat &n_real, const ScalarFloat &n_imag,
                     const Vector3f &wi, const Vector3f &wo,
                     const ScalarFloat &wind_direction) const {
        // Transform directions into azimuthal and zenithal angles
        Float s_theta_i = dr::sqrt(1 - wi.z() * wi.z());
        Float c_theta_i = wi.z();

        Float s_theta_o = dr::sqrt(1 - wo.z() * wo.z());
        Float c_theta_o = wo.z();

        Float phi_i         = dr::atan2(wi.y(), wi.x());
        Float phi_o         = dr::atan2(wo.y(), wo.x());
        Float phi           = phi_i - phi_o;
        auto [s_phi, c_phi] = dr::sincos(phi);

        Float phi_w             = phi_i - wind_direction;
        auto [s_phi_w, c_phi_w] = dr::sincos(phi_w);

        return eval_sun_glint(n_real, n_imag, s_theta_i, c_theta_i, s_theta_o,
                              c_theta_o, s_phi, c_phi, s_phi_w, c_phi_w);
    }

    /**
     * @brief Evaluate the sun glint reflectance.
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
    Float eval_sun_glint(ScalarFloat n_real, ScalarFloat n_imag,
                         const Float &s_i, Float c_i, const Float &s_o,
                         Float c_o, const Float &s_phi, const Float &c_phi,
                         const Float &s_phi_w, const Float &c_phi_w,
                         const bool invert_real = false) const {

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
        auto mask           = Mask(specular_prob < 0.0f);
        specular_prob       = dr::select(mask, 0.0f, specular_prob);

        Float cos_2_chi = c_o * c_i + s_o * s_i * c_phi;

        cos_2_chi = dr::select(cos_2_chi > 1.0f, 0.999999999f, cos_2_chi);
        cos_2_chi = dr::select(cos_2_chi < -1.0f, -0.999999999f, cos_2_chi);

        Float coschi = dr::sqrt(0.5f * (1.0f + cos_2_chi));
        Float sinchi = dr::sqrt(0.5f * (1.0f - cos_2_chi));

        coschi = dr::select(coschi > 1.0f, 0.999999999f, coschi);
        coschi = dr::select(coschi < -1.0f, -0.999999999f, coschi);
        sinchi = dr::select(sinchi > 1.0f, 0.999999999f, sinchi);
        sinchi = dr::select(sinchi < -1.0f, -0.999999999f, sinchi);

        // Invert the real part of the IOR
        if (invert_real) {
            n_real = 1 / n_real;
            n_imag = 0.0f;
        }

        Float fresnel_coeff = eval_fresnel(n_real, n_imag, coschi, sinchi);
        // Sun glint reflectance
        Float num   = dr::Pi<Float> * fresnel_coeff * specular_prob;
        Float denom = 4.0f * c_i * c_o * dr::pow(dr::cos(tilt), 4.0f);
        return num / denom;
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
        Float underlight = (1.0f / (dr::sqr(n_real) + dr::sqr(n_imag))) *
                           (m_r_omega * t_u * t_d) /
                           (1.0f - m_underlight_alpha * m_r_omega);

        return dr::select(outside_range, 0.0f, underlight);
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

        Float coef = 1.0f - (m_c_21 / 2.0f) * (xe2 - 1.0f) * xn -
                     (m_c_03 / 6.0f) * (xn2 - 3.0f) * xn;
        coef = coef + (m_c_40 / 24.0f) * (xe2 * xe2 - 6.0f * xe2 + 3.0f);
        coef = coef + (m_c_04 / 24.0f) * (xn2 * xn2 - 6.0f * xn2 + 3.0f);
        coef = coef + (m_c_22 / 4.0f) * (xe2 - 1.0f) * (xn2 - 1.0f);

        Float prob = coef * dr::InvTwoPi<Float> / (m_sigma_u * m_sigma_c) *
                     dr::exp(-(xe2 + xn2) * 0.5f);
        return prob;
    }

    /**
     * @brief Evaluate the Fresnel coefficient.
     *
     * Evaluates the Fresnel coefficient at the given real and imaginary parts
     * of the index of refraction, and relevant geometry terms derived from the
     * incoming and outgoing directions. The coefficient is computed using the
     * Fresnel equations.
     *
     * @param n_real The real part of the index of refraction.
     * @param n_imag The imaginary part of the index of refraction.
     * @param coschi The cosine of the geometry term.
     * @param sinchi The sine of the geometry term.
     * @return Float The Fresnel coefficient.
     */
    Float eval_fresnel(const ScalarFloat &n_real, const ScalarFloat &n_imag,
                       const Float &coschi, const Float &sinchi) const {

        const ScalarFloat n_real2 = n_real * n_real;
        const ScalarFloat n_imag2 = n_imag * n_imag;

        const Float s = (n_real2) - (n_imag2) - (sinchi * sinchi);

        const Float a_1 = dr::abs(s);
        const Float a_2 = dr::sqrt(dr::sqr(s) + 4.0f * n_real2 * n_imag2);

        const Float u = dr::sqrt(0.5f * dr::abs(a_1 + a_2));
        const Float v = dr::sqrt(0.5f * dr::abs(a_2 - a_1));

        const Float b_1 = (n_real2 - n_imag2) * coschi;
        const Float b_2 = 2 * n_real * n_imag * coschi;

        const Float right_squared =
            (dr::sqr(coschi - u) + v * v) / (dr::sqr(coschi + u) + v * v);
        const Float left_squared = (dr::sqr(b_1 - u) + dr::sqr(b_2 + v)) /
                                   (dr::sqr(b_1 + u) + dr::sqr(b_2 - v));
        const Float R = (right_squared + left_squared) * 0.5f;

        return R;
    }

    /**
     * @brief Compute the correction to the IOR of water.
     *
     * Computes the correction to the index of refraction of water according to
     * the formulas provided by Friedman (1969) and Sverdrup (1942). The
     * correction is computed as a function of the chlorinity of the water.
     *
     * @param chlorinity The chlorinity of the water.
     * @return ScalarFloat The correction to the index of refraction.
     */
    ScalarFloat friedman_sverdrup(const ScalarFloat &chlorinity) const {
        return 0.00017492711f * (0.03f + 1.805f * chlorinity);
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
        ScalarFloat pigment_log = dr::log(pigmentation) / dr::log(10.0f);

        // Backscattering coefficient
        ScalarFloat molecular_scatter_coeff =
            m_ocean_props.molecular_scatter_coeff_6s(wavelength);
        ScalarFloat scattering_coeff = 0.30f * dr::pow(pigmentation, 0.62);
        ScalarFloat backscatter_ratio =
            0.002f +
            0.02f * (0.5f - 0.25f * pigment_log) * (550.0 / wavelength);
        ScalarFloat backscatter_coeff = 0.5f * molecular_scatter_coeff +
                                        scattering_coeff * backscatter_ratio;

        // (Diffuse) attenuation coefficient
        ScalarFloat k          = m_ocean_props.attn_k(wavelength);
        ScalarFloat chi        = m_ocean_props.attn_chi(wavelength);
        ScalarFloat e          = m_ocean_props.attn_e(wavelength);
        ScalarFloat attn_coeff = k + chi * dr::pow(pigmentation, e);

        // If any of the coefficients is zero, we return zero
        if (backscatter_coeff == 0.0f || attn_coeff == 0.0f)
            return 0.0f;

        // Iterative computation of the reflectance
        ScalarFloat u       = 0.75f;
        ScalarFloat r_omega = 0.33f * backscatter_coeff / u / attn_coeff;

        bool converged = false;
        while (!converged) {
            // Update u
            u = (0.9f * (1.0f - r_omega)) / (1.0f + 2.25f * r_omega);

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
    ScalarFloat m_c_21    = 0.f;
    ScalarFloat m_c_03    = 0.f;
    ScalarFloat m_r_omega = 0.f;

    // Underlight parameters
    const ScalarFloat m_underlight_alpha = 0.485f;
};

#endif // OCEAN_PROPS

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
    nodes_x = dr::fmadd(nodes_x, 0.25 * dr::Pi<FloatX>, 0.25 * dr::Pi<FloatX>);
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
        result_p = 1.0f - (td / summ);

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
        m_accel          = props.get<bool>("accel", true);

        // convert from North Left to East Right.
        m_wind_direction = -m_wind_direction + 90.f;
        // m_wind_direction % 360.
        m_wind_direction = m_wind_direction - ( 360.f*floor(m_wind_direction/360.f) );
        // Degree to radians.
        m_wind_direction = m_wind_direction * dr::Pi<ScalarFloat> / 360.f;

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
            ocean_utils.update(m_wavelength, m_wind_speed, m_pigmentation);

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

        m_ocean_utils.update(m_wavelength, m_wind_speed, m_pigmentation);
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
     * @param wi The incident direction of the light.
     * @param wo The outgoing direction of the light.
     * @return Float The sun glint reflectance.
     */
    Float eval_glint(const Vector3f &wi, const Vector3f &wo) const {
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

    /**
     * @brief Evaluate the oceanic reflectance.
     *
     * Evaluates the oceanic reflectance at the provided wavelength, incident
     * and outgoing directions, wind direction, wind speed, chlorinity, and
     * pigmentation. The reflectance is computed by combining the whitecap, sun
     * glint, and underwater light reflectance, according to the coverage of
     * whitecaps (based on 6S).
     *
     * @param ctx The context of evaluation.
     * @param wi The incident direction of the light.
     * @param wo The outgoing direction of the light.
     * @return Float The oceanic reflectance
     */
    Float eval_ocean(const BSDFContext &ctx, const Vector3f &wi,
                     const Vector3f &wo) const {
        // 6S considers the incident direction to come from the light source
        // `TransportMode` has two states:
        //     - `Radiance`, meaning we trace from the sensor to the light
        //     sources
        //     - `Importance`, meaning we trace from the light sources to the
        //     sensor
        //
        Vector3f wi_hat = ctx.mode == TransportMode::Radiance ? wo : wi,
                 wo_hat = ctx.mode == TransportMode::Radiance ? wi : wo;

        Float coverage = m_ocean_utils.eval_whitecap_coverage(m_wind_speed);
        Float whitecap_reflectance   = eval_whitecaps();
        Float glint_reflectance      = eval_glint(wi_hat, wo_hat);
        Float underlight_reflectance = eval_underlight(wi_hat, wo_hat);

        return (coverage * whitecap_reflectance) +
               (1.f - coverage) * glint_reflectance +
               (1.f - (coverage * whitecap_reflectance)) * underlight_reflectance;
    }

    /**
     * @brief Evaluate the Blinn-Phong BRDF.
     *
     * Evaluates the Blinn-Phong BRDF at the provided incident and outgoing
     * directions, and surface normal. The BRDF is computed using the
     * Blinn-Phong distribution.
     *
     * @param wi The incident direction of the light.
     * @param wo The outgoing direction of the light.
     * @param normal The surface normal.
     * @return Float The Blinn-Phong BRDF.
     * @note This function is not used, but serves as a reference for future
     * implementations.
     */
    Float eval_blinn_phong(const Vector3f &wi, const Vector3f &wo,
                           const Normal3f &normal) const {
        Float coverage = m_ocean_utils.eval_whitecap_coverage(m_wind_speed);
        Float factor   = (m_shininess + 2.f) / (2.f * dr::Pi<Float>);
        Vector3f half  = dr::normalize(wi + wo);
        Float dot      = dr::dot(half, normal);

        // Blinn-phong => clamp dot to above zero
        dot = dr::clamp(dot, 0.f, 1.f);

        Float phong = dr::pow(dot, m_shininess);

        return coverage + (1 - coverage) * factor * phong;
    }

    std::pair<BSDFSample3f, Spectrum>
    sample(const BSDFContext &ctx, const SurfaceInteraction3f &si,
           Float sample1, const Point2f &sample2, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);

        bool has_diffuse  = ctx.is_enabled(BSDFFlags::DiffuseReflection, 0);
        bool has_specular = ctx.is_enabled(BSDFFlags::GlossyReflection, 1);

        Float cos_theta_i = Frame3f::cos_theta(si.wi);
        BSDFSample3f bs   = dr::zeros<BSDFSample3f>();
        active &= cos_theta_i > 0.f;
        if (unlikely(dr::none_or<false>(active)) ||
            (!has_diffuse && !has_specular))
            return { bs, 0.f };

        Float coverage  = m_ocean_utils.eval_whitecap_coverage(m_wind_speed);
        Float prob_diff = coverage;

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
            dr::masked(bs.wo, sample_glint)                = wo;
            dr::masked(bs.sampled_component, sample_glint) = 1;
            dr::masked(bs.sampled_type, sample_glint) =
                +BSDFFlags::GlossyReflection;
        }

        bs.pdf = pdf(ctx, si, bs.wo, active);
        bs.eta = 1.f;

        UnpolarizedSpectrum value =
            eval_ocean(ctx, si.wi, bs.wo) * Frame3f::cos_theta(bs.wo) / bs.pdf;
        return { bs,
                 (depolarizer<Spectrum>(value)) & (active && bs.pdf > 0.f) };
    }

    Spectrum eval(const BSDFContext &ctx, const SurfaceInteraction3f &si,
                  const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        bool has_diffuse = ctx.is_enabled(BSDFFlags::DiffuseReflection, 0);
        bool has_glint   = ctx.is_enabled(BSDFFlags::GlossyReflection, 1);

        if (unlikely(dr::none_or<false>(active) || !has_glint && !has_diffuse))
            return 0.f;

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        // Ensure incoming and outgoing directions are in the upper hemisphere
        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

        // Compute the whitecap reflectance
        UnpolarizedSpectrum result(0.f);
        UnpolarizedSpectrum blinn(0.f);
        UnpolarizedSpectrum whitecap_reflectance(0.f);
        UnpolarizedSpectrum glint_reflectance(0.f);
        UnpolarizedSpectrum underlight_reflectance(0.f);

        // Get the reflected directions
        auto is_reflect =
            Mask(dr::eq(dr::sign(cos_theta_i), dr::sign(cos_theta_o))) &&
            active;

        // 6SV considers the incident direction to come from the light source
        // `TransportMode` has two states:
        //     - `Radiance`, meaning we trace from the sensor to the light
        //     sources
        //     - `Importance`, meaning we trace from the light sources to the
        //     sensor
        Vector3f wi_hat = ctx.mode == TransportMode::Radiance ? wo : si.wi,
                 wo_hat = ctx.mode == TransportMode::Radiance ? si.wi : wo;

        if (has_diffuse) {
            // If whitecaps are enabled, compute the whitecap reflectance
            whitecap_reflectance   = eval_whitecaps();
            underlight_reflectance = eval_underlight(wi_hat, wo_hat);
        }

        if (has_glint)
            // If sun glint is enabled, compute the glint reflectance
            glint_reflectance = eval_glint(wi_hat, wo_hat);

        // Combine the results
        Float coverage = m_ocean_utils.eval_whitecap_coverage(m_wind_speed);

        // For debugging purposes, the channel indicates what term of the BRDF
        // to evaluate
        switch (m_component) {
            case 1:
                result[is_reflect] = (whitecap_reflectance) &active;
                break;
            case 2:
                result[is_reflect] =
                    (1.f - coverage) * glint_reflectance & active;
                break;
            case 3:
                result[is_reflect] =
                    (1.f - (whitecap_reflectance)) * underlight_reflectance &
                    active;
                break;
            case 4:
                result[is_reflect] = coverage;
                break;
            default: {
                result[is_reflect] =
                    whitecap_reflectance + (1.f - coverage) * glint_reflectance +
                    (1.f - (whitecap_reflectance)) * underlight_reflectance;

                // Cosine foreshortening factor - not sure if this is actually
                // needed?
                result[is_reflect] *= cos_theta_o;
            } break;
        }

        return depolarizer<Spectrum>(result) & active;
    }

    Float pdf(const BSDFContext &ctx, const SurfaceInteraction3f &si,
              const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);
        bool has_diffuse  = ctx.is_enabled(BSDFFlags::DiffuseReflection, 0),
             has_specular = ctx.is_enabled(BSDFFlags::GlossyReflection, 1);

        if (unlikely((!has_diffuse && !has_specular) ||
                     dr::none_or<false>(active)))
            return 0.f;

        Float coverage = m_ocean_utils.eval_whitecap_coverage(m_wind_speed);
        Float weight_diffuse = coverage, weight_specular = (1 - coverage);

        // Check if the normal has only zeros. If this is the case, use a
        // default normal
        Vector3f normal   = si.n;
        Mask degen_normal = dr::all(dr::eq(normal, Vector3f(0.f)));
        dr::masked(normal, degen_normal) = Vector3f(0.f, 0.f, 1.f);

        Vector3f half    = dr::normalize(si.wi + wo);
        Float projection = dr::dot(half, normal);
        Float D          = ((m_shininess + 2.0f) /
                   dr::TwoPi<Float>) *dr::pow(projection, m_shininess);

        // We multiply the probability of the specular lobe with the pdf of
        // the Blinn-Phong distribution and the probability of the diffuse lobe
        // with the pdf of the cosine-weighted hemisphere.
        Float pdf_diffuse  = warp::square_to_cosine_hemisphere_pdf(wo),
              pdf_specular = (D * projection) / (4.0f * dr::dot(si.wi, half));

        Float pdf =
            weight_diffuse * pdf_diffuse + weight_specular * pdf_specular;

        // If the outgoing direction is in the lower hemisphere, we return zero
        Float cos_theta_o = Frame3f::cos_theta(wo);

        // return dr::select(cos_theta_o > 0.f, 1.f, 0.f);
        return dr::select(cos_theta_o > 0.f, pdf, 0.f);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "OceanLegacy[" << std::endl
            << "  component = " << string::indent(m_component) << std::endl
            << "  wavelength = " << string::indent(m_wavelength) << std::endl
            << "  wind_speed = " << string::indent(m_wind_speed) << std::endl
            << "  wind_direction = " << string::indent(m_wind_direction)
            << std::endl
            << "  chlorinity = " << string::indent(m_chlorinity) << std::endl
            << "  pigmentation = " << string::indent(m_pigmentation)
            << std::endl
            << "  shininess = " << string::indent(m_shininess) << std::endl
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
