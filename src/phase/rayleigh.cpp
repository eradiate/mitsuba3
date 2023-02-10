#include <mitsuba/core/distr_1d.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/phase.h>

/**!

.. _phase-rayleigh:

Rayleigh phase function (:monosp:`rayleigh`)
-----------------------------------------------

Scattering by particles that are much smaller than the wavelength
of light (e.g. individual molecules in the atmosphere) is well-approximated
by the Rayleigh phase function. This plugin implements an unpolarized
version of this scattering model (*i.e.* the effects of polarization are
ignored). This plugin is useful for simulating scattering in planetary
atmospheres.

This model has no parameters.

.. tabs::
    .. code-tab:: xml

        <phase type="rayleigh" />

    .. code-tab:: python

        'type': 'rayleigh'

*/

NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
class RayleighPhaseFunction final : public PhaseFunction<Float, Spectrum> {
public:
    MI_IMPORT_BASE(PhaseFunction, m_flags)
    MI_IMPORT_TYPES(PhaseFunctionContext)

    RayleighPhaseFunction(const Properties &props) : Base(props) {
        m_depolarization = props.get<ScalarFloat>("depolarization", 0.f);
            
        if (m_depolarization < 0.f || m_depolarization >= 1.f)
            Log(Error, "Depolarization factor must be in [0, 1[");

        m_flags = +PhaseFunctionFlags::Anisotropic;
    }

    MI_INLINE Float eval_rayleigh_pdf(Float cos_theta) const {
        Float rho = (Float) m_depolarization,
              r1 = (1.f - rho) / (1.f + rho / 2.f),
              r2 = (1.f + rho) / (1.f - rho);
        // cos_theta in physics convention
        return (3.f / 16.f) * dr::InvPi<Float> * r1 * (r2 + dr::sqr(cos_theta));
    }

    MI_INLINE Spectrum eval_rayleigh(const PhaseFunctionContext &ctx, 
                                      const MediumInteraction3f &mi,
                                      const Vector3f &wo,
                                      Float cos_theta) const {
        Spectrum phase_val;

        if constexpr (is_polarized_v<Spectrum>) {
            /* cos_theta in physics convention */
            phase_val = mueller::rayleigh_scatter(cos_theta, (Float) m_depolarization);

            /* Due to the coordinate system rotations for polarization-aware
                phase below we need to know the propagation direction of light.
                In the following, light arrives along `-wo_hat` and leaves along
                `+wi_hat`. */
            Vector3f wo_hat = ctx.mode == TransportMode::Radiance ? wo : mi.wi,
                     wi_hat = ctx.mode == TransportMode::Radiance ? mi.wi : wo;

            /* The Stokes reference frame vector of this matrix lies in the 
                scattering plane spanned by wi and wo. */
            Vector3f x_hat = dr::cross(-wo_hat, wi_hat),
                     p_axis_in = dr::normalize(dr::cross(x_hat, -wo_hat)),
                     p_axis_out = dr::normalize(dr::cross(x_hat, wi_hat));

            /* Rotate in/out reference vector of weight s.t. it aligns with the
            implicit Stokes bases of -wo_hat & wi_hat. */
            phase_val = mueller::rotate_mueller_basis(phase_val,
                                                     -wo_hat, p_axis_in, mueller::stokes_basis(-wo_hat),
                                                      wi_hat, p_axis_out, mueller::stokes_basis(wi_hat));

            // If the cross product x_hat is too small, phase_val may be NaN
            dr::masked(phase_val, dr::isnan(phase_val)) = 0.f;

        } else {
            Float rho = (Float) m_depolarization,
                  r1 = (1.f - rho) / (1.f + rho / 2.f),
                  r2 = (1.f + rho) / (1.f - rho);

            phase_val = (3.f / 16.f) * dr::InvPi<Float> * r1 * (r2 + dr::sqr(cos_theta));
        }

        return phase_val;
    }

    std::pair<Vector3f, Spectrum> sample(const PhaseFunctionContext &ctx,
                                         const MediumInteraction3f &mi,
                                         Float /* sample1 */,
                                         const Point2f &sample,
                                         Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::PhaseFunctionSample, active);

        /* The importance sampling below corresponds to the "direct sampling of
        spherical coordinate" scheme described by
        :cite:`Frisvad2011ImportanceSamplingRayleigh` (section 3.B) */ 
        Float z   = 2.f * (2.f * sample.x() - 1.f);
        Float tmp = dr::sqrt(dr::sqr(z) + 1.f);
        Float A   = dr::cbrt(-z + tmp);
        Float B   = dr::cbrt(-z - tmp);
        Float cos_theta = A + B; /* cos_theta in physics convention */
        Float sin_theta = dr::safe_sqrt(1.0f - dr::sqr(cos_theta));
        auto [sin_phi, cos_phi] = dr::sincos(dr::TwoPi<Float> * sample.y());

        /* if θ is the scattering angle in physics convention, and θ' 
        the scattering angle in graphics convention, then θ' = π - θ
        and cos(θ') = -cos(θ) and sin(θ') = sin(θ)*/
        Vector3f wo = Vector3f( sin_theta * cos_phi, sin_theta * sin_phi, -cos_theta );

        /* eval_rayleigh expects cos(θ) in physics convention */
        Spectrum phase_val = eval_rayleigh(ctx, mi, wo, cos_theta);

        return { wo, phase_val };
    }

    Spectrum eval(const PhaseFunctionContext &ctx,
               const MediumInteraction3f &mi, const Vector3f &wo,
               Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::PhaseFunctionEvaluate, active);
        /* if the incident direction is ω in graphics convention, it
        is -ω in physics convention.
        eval_rayleigh expects cos(θ) in physics convention  */
        return eval_rayleigh(ctx, mi, wo, dot(wo, -mi.wi));
    }

    Float pdf(const PhaseFunctionContext &/* ctx */,
              const MediumInteraction3f &mi, const Vector3f &wo,
               Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::PhaseFunctionEvaluate, active);
        /* if the incident direction is ω in graphics convention, it
        is -ω in physics convention.
        eval_rayleigh_pdf expects cos(θ) in physics convention  */
        return eval_rayleigh_pdf(dot(wo, -mi.wi));
    }

    
    std::string to_string() const override { return "RayleighPhaseFunction[]"; }

    MI_DECLARE_CLASS()
private:
    ScalarFloat m_depolarization;
};

MI_IMPLEMENT_CLASS_VARIANT(RayleighPhaseFunction, PhaseFunction)
MI_EXPORT_PLUGIN(RayleighPhaseFunction, "Rayleigh phase function")
NAMESPACE_END(mitsuba)
