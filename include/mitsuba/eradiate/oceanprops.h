#pragma once

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/distr_1d.h>
#include <mitsuba/core/math.h>
#include <mitsuba/render/fwd.h>
#include <drjit/complex.h>

NAMESPACE_BEGIN(mitsuba)

// Header content - THIS CAN BE MOVED IF NECESSARY
#ifndef OCEAN_PROPS
#define OCEAN_PROPS

template <typename Float, typename Spectrum> class OceanProperties {
public:
    MI_IMPORT_TYPES()

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
    ScalarFloat molecular_scatter_coeff_6s(const ScalarFloat &wavelength) const {
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
template<typename Float>
Float whitecap_coverage_monahan(const Float &wind_speed) {
    static constexpr dr::scalar_t<Float> m_monahan_alpha  = 2.95e-06f;
    static constexpr dr::scalar_t<Float> m_monahan_lambda = 3.52f;
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
 * @param wavelength The wavelength at which to evaluate the reflectance in nm.
 * @param wind_speed The wind speed at which to evaluate the reflectance.
 * @return ScalarFloat The reflectance of whitecaps.
 */
template <typename Float>
Float whitecap_reflectance_frouin(const Float &wavelength,
                                  const Float &wind_speed) {
    // Compute the fractional coverage of whitecaps
    Float coverage = whitecap_coverage_monahan(wind_speed);

    Float eff_reflectance = dr::select(
        wavelength*0.001f >= 0.6f,
        0.22f * dr::exp(-1.75f * dr::pow(wavelength*0.001f - 0.6f, 0.99f)), 0.22f);

    // Compute the whitecap reflectance
    Float whitecap_reflectance = coverage * eff_reflectance;

    return whitecap_reflectance;
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
template<typename Float>
Float friedman_sverdrup(const Float &chlorinity) {
    return 0.00017492711f * (0.03f + 1.805f * chlorinity);
}

/**
 * @brief Evaluate the complex index of refraction of the ocean.
 *
 * @param wavelength The wavelength at which to evaluate the reflectance.
 * @param chlorinity The chlorinity of the water.
 * @return The complex index of refraction of water.
 */
template<typename Float, typename Spectrum, typename ScalarFloat>
std::pair<ScalarFloat, ScalarFloat>
water_ior( const OceanProperties<Float, Spectrum> &ocean_props,
           const ScalarFloat &wavelength,
           const ScalarFloat &chlorinity) {
    ScalarFloat n_real =
        ocean_props.ior_real(wavelength) + friedman_sverdrup<ScalarFloat>(chlorinity);
    ScalarFloat n_imag = ocean_props.ior_cplx(wavelength);

    return { n_real, n_imag };
}

/**
 * @brief Evaluate the Fresnel coefficient as implemented by 6SV
 * (based on Born and Wolf (1975), which we couldn't verify).
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
template<typename Float>
Float fresnel_sunglint_legacy(
    const Float &n_real,
    const Float &n_imag,
    const Float &coschi,
    const Float &sinchi) {

    const Float n_real2 = n_real * n_real;
    const Float n_imag2 = n_imag * n_imag;

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

    return (right_squared + left_squared) * 0.5f;
};

/**
 * @brief Evaluate the polarized Fresnel coefficient as implemented by Mischenko (1997).
 *
 * Evaluates the polarized Fresnel coefficient at given complex indices of
 * refraction and geometry defined by the incoming and outgoing directions.
 * The mueller matrix is computed using the formulation from Mishchenko (1997).
 *
 * @param n_ext The external complex index of refraction.
 * @param n_water The water complex index of refraction.
 * @param wi The incident direction, in the direction of progagation (physics convention).
 * @param wo The outgoing direciton, in the direction of propagation (physics convention).
 * @return MuellerMatrix<UnpolarizedSpectrum> The Fresnel mueller matrix.
 */
template<typename Float, typename UnpolarizedSpectrum>
MuellerMatrix<UnpolarizedSpectrum> fresnel_sunglint_polarized(
    const dr::Complex<UnpolarizedSpectrum> &n_ext,
    const dr::Complex<UnpolarizedSpectrum> &n_water,
    Vector<Float, 3> wi,
    Vector<Float, 3> wo){

    using Mask = dr::mask_t<Float>;
    using Complex2u = dr::Complex<UnpolarizedSpectrum>;
    using Vector3f = Vector<Float, 3>;

    const Complex2u n1 = n_ext;
    const Complex2u n2 = n_water;

    Float mu_i = dr::abs(wi.z());
    Float mu_o = dr::abs(wo.z());

    Float phi_i = -dr::atan2(wi.y(), wi.x());
    Float phi_o = -dr::atan2(wo.y(), wo.x());

    dr::masked(mu_i, mu_i > 0.9999999f) = 0.9999999f;
    dr::masked(mu_o, mu_o > 0.9999999f) = 0.9999999f;

    /* unit vectors of incident and reflected rays */
    Float sitheta_i = dr::sqrt(1.f - mu_i * mu_i);
    Float sitheta_o = dr::sqrt(1.f - mu_o * mu_o);

    wi = Vector3f(sitheta_i * dr::cos(phi_i), sitheta_i * dr::sin (phi_i), -mu_i);
    wo = Vector3f(sitheta_o * dr::cos(phi_o), sitheta_o * dr::sin (phi_o), mu_o);

    // local surface normal k_d
    const Vector3f k_d = wi - wo;
    const Float k_d_norm2 = dr::dot(k_d, k_d);

    // Fresnel reflection coefficition (should that be taken out?)
    // the incident angel wrt local surface normal
    const Float mu_i_l = dr::dot(k_d, wi) / dr::sqrt(k_d_norm2);
    // mu_refr_l, R_r and R_l are all complex. R_r for Fresnel perpendicular and R_l for fresnel parallel
    const Complex2u mu_refr_l = dr::sqrt(1.f - (1.f - mu_i_l * mu_i_l) * n1 * n1 / (n2 * n2));
    const Complex2u R_r = (n1 * mu_i_l - n2 * mu_refr_l) / (n1 * mu_i_l + n2 * mu_refr_l);
    const Complex2u R_l = (n2 * mu_i_l - n1 * mu_refr_l) / (n2 * mu_i_l + n1 * mu_refr_l);

    // theta_v, phi_v, and w represent the polarization frame
    const Vector3f z(0.f, 0.f, 1.f);

    Mask collinear = dr::all(dr::eq(wi, -z));
    const Vector3f phi_v_i = dr::select(collinear, Vector3f(0.f, 1.f, 0.f), dr::normalize(dr::cross(z, wi)));
    const Vector3f theta_v_i = dr::cross(phi_v_i,wi);

    collinear = dr::all(dr::eq(wo, z));
    const Vector3f phi_v_o = dr::select(collinear, Vector3f(0.f,1.f,0.f), dr::normalize(dr::cross(z, wo)));
    const Vector3f theta_v_o = dr::cross(phi_v_o, wo);

    // amplitude scattering matrix
    const Float pi_dot_wo = dr::dot(phi_v_i, wo);
    const Float po_dot_wi = dr::dot(phi_v_o, wi);
    const Float ti_dot_wo = dr::dot(theta_v_i, wo);
    const Float to_dot_wi = dr::dot(theta_v_o, wi);

    const Complex2u f_tt =  pi_dot_wo * po_dot_wi * R_r + ti_dot_wo * to_dot_wi * R_l;
    const Complex2u f_tp = -ti_dot_wo * po_dot_wi * R_r + pi_dot_wo * to_dot_wi * R_l;
    const Complex2u f_pt = -pi_dot_wo * to_dot_wi * R_r + ti_dot_wo * po_dot_wi * R_l;
    const Complex2u f_pp =  ti_dot_wo * to_dot_wi * R_r + pi_dot_wo * po_dot_wi * R_l;

    // stokes transmission matrix
    const Vector3f wi_cross_wo = dr::cross(wi, wo);
    collinear = dr::all(dr::eq(wo, -wi));
    const Float norm2 = dr::select(collinear,
                                    0.000001f,
                                    dr::pow(dr::dot(wi_cross_wo, wi_cross_wo),2));
    Float coeff = 1.f / norm2;

    const UnpolarizedSpectrum n_f_tt = dr::abs(f_tt);
    const UnpolarizedSpectrum n_f_tp = dr::abs(f_tp);
    const UnpolarizedSpectrum n_f_pt = dr::abs(f_pt);
    const UnpolarizedSpectrum n_f_pp = dr::abs(f_pp);

    const UnpolarizedSpectrum M00 = 0.5 * coeff * (n_f_tt * n_f_tt + n_f_tp * n_f_tp + n_f_pt * n_f_pt + n_f_pp * n_f_pp);
    const UnpolarizedSpectrum M01 = 0.5 * coeff * (n_f_tt * n_f_tt - n_f_tp * n_f_tp + n_f_pt * n_f_pt - n_f_pp * n_f_pp);
    const UnpolarizedSpectrum M10 = 0.5 * coeff * (n_f_tt * n_f_tt + n_f_tp * n_f_tp - n_f_pt * n_f_pt - n_f_pp * n_f_pp);
    const UnpolarizedSpectrum M11 = 0.5 * coeff * (n_f_tt * n_f_tt - n_f_tp * n_f_tp - n_f_pt * n_f_pt + n_f_pp * n_f_pp);
    const UnpolarizedSpectrum M02 = -coeff * dr::real(f_tt*dr::conj(f_tp) + f_pt*dr::conj(f_pp));
    const UnpolarizedSpectrum M03 = -coeff * dr::imag(f_tt*dr::conj(f_tp) + f_pt*dr::conj(f_pp));
    const UnpolarizedSpectrum M12 = -coeff * dr::real(f_tt*dr::conj(f_tp) - f_pt*dr::conj(f_pp));
    const UnpolarizedSpectrum M13 = -coeff * dr::imag(f_tt*dr::conj(f_tp) - f_pt*dr::conj(f_pp));
    const UnpolarizedSpectrum M20 = -coeff * dr::real(f_tt*dr::conj(f_pt) + f_tp*dr::conj(f_pp));
    const UnpolarizedSpectrum M21 = -coeff * dr::real(f_tt*dr::conj(f_pt) - f_tp*dr::conj(f_pp));
    const UnpolarizedSpectrum M22 =  coeff * dr::real(f_tt*dr::conj(f_pp) + f_tp*dr::conj(f_pt));
    const UnpolarizedSpectrum M23 =  coeff * dr::imag(f_tt*dr::conj(f_pp) - f_tp*dr::conj(f_pt));
    const UnpolarizedSpectrum M30 =  coeff * dr::imag(f_tt*dr::conj(f_pt) + f_tp*dr::conj(f_pp));
    const UnpolarizedSpectrum M31 =  coeff * dr::imag(f_tt*dr::conj(f_pt) - f_tp*dr::conj(f_pp));
    const UnpolarizedSpectrum M32 = -coeff * dr::imag(f_tt*dr::conj(f_pp) + f_tp*dr::conj(f_pt));
    const UnpolarizedSpectrum M33 =  coeff * dr::real(f_tt*dr::conj(f_pp) - f_tp*dr::conj(f_pt));

    return MuellerMatrix<UnpolarizedSpectrum>(
        M00, M01, M02, M03,
        M10, M11, M12, M13,
        M20, M21, M22, M23,
        M30, M31, M32, M33
    );
};

/**
 * @brief Mean square slope *squared* as described by Cox and Munk (1954).
 *
 * @param wind_speed speed of wind at sea surface (mast height, 10 m) [m/s].
 * @return tuple (Float, Float)
 * Respectively the cross wind and along wind mean square slope squared.
 */
template<typename Float>
std::tuple<Float, Float> cox_munk_crosswind_upwind(const Float& wind_speed) {
    Float sigma_cross_2 = dr::fmadd(wind_speed, 0.00192f, 0.003f);
    Float sigma_along_2 = 0.00316f * wind_speed;

    return { sigma_cross_2, sigma_along_2 };
}

/**
 * @brief Isotropic Mean Squared Slope *squared* as described by Cox and Munk
 * (1954)
 *
 * @param wind_speed speed of wind at sea surface (mast height, 10 m) [m/s].
 * @return Float
 * The isotropic wind mean square slope squared
 */
template<typename Float>
Float cox_munk_msslope_squared(const Float& wind_speed) {
    return dr::fmadd(wind_speed, 0.00512f, 0.003f);
}

/**
 * @brief Evaluates the anisotropic distribution of Cox and Munk (1954), with
 * the Gram Charlier series expansion.
 *
 * @param wind_direction direction in radians, east right convention.
 * @param wind_speed wind_speed speed of wind at sea surface (mast height, 10 m)
 * [m/s].
 * @param sigma_u Upwind root mean slope distribution.
 * @param sigma_c Crosswind root mean slope distribution.
 * @param m Half vector.
 *
 */
template <typename Float>
Float cox_munk_anisotropic_distrib(const Float &wind_direction,
                             const Float &wind_speed, const Float &sigma_u,
                             const Float &sigma_c, const Vector<Float, 3> &m) {

    using Vector3f    = Vector<Float, 3>;
    using ScalarFloat = dr::scalar_t<Float>;

    // Distribution constants
    static constexpr ScalarFloat c_40 = 0.40f;
    static constexpr ScalarFloat c_22 = 0.12f;
    static constexpr ScalarFloat c_04 = 0.23f;

    // Distribution variables
    Float c_21 = 0.01f - 0.0086f * wind_speed;
    Float c_03 = 0.04f - 0.033f * wind_speed;

    auto [s_phi, c_phi] = dr::sincos(wind_direction);

    Vector3f m_p = Vector3f(c_phi * m.x() + s_phi * m.y(),
                            -s_phi * m.x() + c_phi * m.y(), m.z());
    m_p          = dr::normalize(m_p);

    const Float xn  = m_p.x() * dr::rcp(sigma_u * m_p.z());
    const Float xe  = m_p.y() * dr::rcp(sigma_c * m_p.z());
    const Float xe2 = xe * xe;
    const Float xn2 = xn * xn;

    Float coef =
        1.f - (c_21 / 2.f) * (xe2 - 1.f) * xn - (c_03 / 6.f) * (xn2 - 3.f) * xn;
    coef = coef + (c_40 / 24.f) * (xe2 * xe2 - 6.f * xe2 + 3.f);
    coef = coef + (c_04 / 24.f) * (xn2 * xn2 - 6.f * xn2 + 3.f);
    coef = coef + (c_22 / 4.f) * (xe2 - 1.f) * (xn2 - 1.f);

    Float prob = coef * dr::InvTwoPi<Float> * dr::rcp(sigma_u * sigma_c) *
                 dr::exp(-(xe2 + xn2) * 0.5f);
    return prob;
}

/**
 * Evaluates the gram charlier coefficient for the Cox and Munk (1954) only.
 */
template <typename Float>
Float cox_munk_gram_charlier_coef(const Float &wind_direction,
                            const Float &wind_speed, const Float &sigma_u,
                            const Float &sigma_c, const Vector<Float, 3> &m) {

    using Vector3f    = Vector<Float, 3>;
    using ScalarFloat = dr::scalar_t<Float>;

    // Distribution constants
    static ScalarFloat c_40 = 0.40f;
    static ScalarFloat c_22 = 0.12f;
    static ScalarFloat c_04 = 0.23f;

    // Distribution variables
    Float c_21 = 0.01f - 0.0086f * wind_speed;
    Float c_03 = 0.04f - 0.033f * wind_speed;

    auto [s_phi, c_phi] = dr::sincos(wind_direction);

    Vector3f m_p = Vector3f(c_phi * m.x() + s_phi * m.y(),
                            -s_phi * m.x() + c_phi * m.y(), m.z());
    m_p          = dr::normalize(m_p);

    const Float xn = m_p.x() / (sigma_u * m_p.z());
    const Float xe = m_p.y() / (sigma_c * m_p.z());

    const Float xe2 = xe * xe;
    const Float xn2 = xn * xn;

    Float coef =
        1.f - (c_21 / 2.f) * (xe2 - 1.f) * xn - (c_03 / 6.f) * (xn2 - 3.f) * xn;
    coef = coef + (c_40 / 24.f) * (xe2 * xe2 - 6.f * xe2 + 3.f);
    coef = coef + (c_04 / 24.f) * (xn2 * xn2 - 6.f * xn2 + 3.f);
    coef = coef + (c_22 / 4.f) * (xe2 - 1.f) * (xn2 - 1.f);

    return coef;
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
template <typename Float, typename Spectrum, typename ScalarFloat>
ScalarFloat r_omega(const OceanProperties<Float, Spectrum> ocean_props,
                    const ScalarFloat &wavelength,
                    const ScalarFloat &pigmentation) {

    ScalarFloat pigment_log = dr::log(pigmentation) / dr::log(10.f);

    // Backscattering coefficient
    ScalarFloat molecular_scatter_coeff =
        ocean_props.molecular_scatter_coeff_6s(wavelength);
    ScalarFloat scattering_coeff = 0.30f * dr::pow(pigmentation, 0.62);
    ScalarFloat backscatter_ratio =
        0.002f + 0.02f * (0.5f - 0.25f * pigment_log) * (550.f / wavelength);
    ScalarFloat backscatter_coeff =
        0.5f * molecular_scatter_coeff + scattering_coeff * backscatter_ratio;

    // (Diffuse) attenuation coefficient
    ScalarFloat k          = ocean_props.attn_k(wavelength);
    ScalarFloat chi        = ocean_props.attn_chi(wavelength);
    ScalarFloat e          = ocean_props.attn_e(wavelength);
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
        ScalarFloat r_omega_new = 0.33f * backscatter_coeff / (u * attn_coeff);

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

#endif // OCEAN_PROPS

NAMESPACE_END(mitsuba)
