#include <mitsuba/core/properties.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/phase.h>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _phase-isotropic:

Isotropic phase function (:monosp:`isotropic`)
----------------------------------------------

This phase function simulates completely uniform scattering,
where all directionality is lost after a single scattering
interaction. It does not have any parameters.

.. tabs::
    .. code-tab:: xml

        <phase type="isotropic" />

    .. code-tab:: python

        'type': 'isotropic'
*/

template <typename Float, typename Spectrum>
class IsotropicPhaseFunction final : public PhaseFunction<Float, Spectrum> {
public:
    MI_IMPORT_BASE(PhaseFunction, m_flags, m_components)
    MI_IMPORT_TYPES(PhaseFunctionContext)

    IsotropicPhaseFunction(const Properties & props) : Base(props) {
        m_flags = +PhaseFunctionFlags::Isotropic;
        dr::set_attr(this, "flags", m_flags);
        m_components.push_back(m_flags);
    }

    std::pair<Vector3f, Spectrum> sample(const PhaseFunctionContext & /* ctx */,
                                         const MediumInteraction3f & /* mi */,
                                         Float /* sample1 */,
                                         const Point2f &sample2,
                                         Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::PhaseFunctionSample, active);

        Vector3f wo  = warp::square_to_uniform_sphere(sample2);
        UnpolarizedSpectrum value = warp::square_to_uniform_sphere_pdf(wo);
        return { wo, depolarizer<Spectrum>(value) };
    }

    Spectrum eval(const PhaseFunctionContext & /* ctx */, const MediumInteraction3f & /* mi */,
                  const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::PhaseFunctionEvaluate, active);
        return UnpolarizedSpectrum(warp::square_to_uniform_sphere_pdf(wo));
    }

    Float pdf(const PhaseFunctionContext & /* ctx */, const MediumInteraction3f & /* mi */,
               const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::PhaseFunctionEvaluate, active);
        return warp::square_to_uniform_sphere_pdf(wo);
    }

    std::string to_string() const override { return "IsotropicPhaseFunction[]"; }

    MI_DECLARE_CLASS()
private:
};

MI_IMPLEMENT_CLASS_VARIANT(IsotropicPhaseFunction, PhaseFunction)
MI_EXPORT_PLUGIN(IsotropicPhaseFunction, "Isotropic phase function")
NAMESPACE_END(mitsuba)
