#include <mitsuba/render/texture.h>
#include <mitsuba/render/interaction.h>
#include <mitsuba/core/properties.h>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _spectrum-delta:

Delta spectrum (:monosp:`delta`)
--------------------------------

.. pluginparameters::

 * - wavelength
   - |float|
   -
   - |exposed|

 * - value
   - |float|
   -
   - |exposed|

 */

template <typename Float, typename Spectrum>
class DeltaSpectrum final : public Texture<Float, Spectrum> {
public:
    MI_IMPORT_TYPES(Texture)

    DeltaSpectrum(const Properties &props) : Texture(props) {
        m_wavelength = dr::opaque<Float>(props.get<ScalarFloat>("wavelength"));
        m_value = dr::opaque<Float>(props.get<ScalarFloat>("value", 1.f));
    }

    void traverse(TraversalCallback *cb) override {
        cb->put("value", m_value, ParamFlags::NonDifferentiable);
        cb->put("wavelength", m_value, ParamFlags::NonDifferentiable);
    }

    // void parameters_changed(const std::vector<std::string> &/*keys*/) override {
    // }

    UnpolarizedSpectrum eval(const SurfaceInteraction3f & /*si*/,
                             Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::TextureEvaluate, active);
        return m_value;
    }

    Float eval_1(const SurfaceInteraction3f & /*it*/, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::TextureEvaluate, active);
        return m_value;
    }

    Color3f eval_3(const SurfaceInteraction3f & /*it*/, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::TextureEvaluate, active);
        return m_value;
    }

    Vector2f eval_1_grad(const SurfaceInteraction3f & /*it*/, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::TextureEvaluate, active);
        return m_value;
    }

    Wavelength pdf_spectrum(const SurfaceInteraction3f & /*si*/, Mask /*active*/) const override {
        NotImplementedError("pdf");
    }

    std::pair<Wavelength, UnpolarizedSpectrum>
    sample_spectrum(const SurfaceInteraction3f & /*si*/,
                    const Wavelength & sample, Mask /*active*/) const override {
        if constexpr (is_spectral_v<Spectrum>) {
            return { m_wavelength, m_value };
        } else {
            DRJIT_MARK_USED(sample);
            return { dr::empty<Wavelength>(), m_value };
        }
    }

    Float mean() const override { return m_value; }

    ScalarVector2f wavelength_range() const override {
        return { dr::slice(m_wavelength), dr::slice(m_wavelength) };
    }

    ScalarFloat spectral_resolution() const override {
        return 0.f;
    }

    ScalarFloat max() const override {
        return dr::slice(dr::max(m_value));
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "RegularSpectrum[" << std::endl
            << "  wavelength = " << string::indent(m_wavelength) << "," << std::endl
            << "  value = " << string::indent(m_value) << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS(DeltaSpectrum)
private:
    Float m_value;
    Float m_wavelength;

    MI_TRAVERSE_CB(Texture, m_value, m_wavelength)
};

MI_EXPORT_PLUGIN(DeltaSpectrum)
NAMESPACE_END(mitsuba)
