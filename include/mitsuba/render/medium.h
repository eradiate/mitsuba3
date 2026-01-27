#pragma once

#include <mitsuba/core/object.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/traits.h>
#include <mitsuba/render/fwd.h>
#include <drjit/call.h>

NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
class MI_EXPORT_LIB Medium : public JitObject<Medium<Float, Spectrum>> {
public:
    MI_IMPORT_TYPES(PhaseFunction, Sampler, Scene, Texture);

    // #ERADIATE_CHANGE_BEGIN: Forward declaration for extremum structure
    // using ExtremumStructure = typename ExtremumStructure<Float, Spectrum>;
    // #ERADIATE_CHANGE_END

    /// Destructor
    ~Medium();

    /// Intersects a ray with the medium's bounding box
    virtual std::tuple<Mask, Float, Float>
    intersect_aabb(const Ray3f &ray) const = 0;

    /// Returns the medium's majorant used for delta tracking
    virtual UnpolarizedSpectrum
    get_majorant(const MediumInteraction3f &mi,
                 Mask active = true) const = 0;

    /// Returns the medium coefficients Sigma_s, Sigma_n and Sigma_t evaluated
    /// at a given MediumInteraction mi
    virtual std::tuple<UnpolarizedSpectrum, UnpolarizedSpectrum,
                       UnpolarizedSpectrum>
    get_scattering_coefficients(const MediumInteraction3f &mi,
                                Mask active = true) const = 0;

    /**
     * \brief Sample a free-flight distance in the medium.
     *
     * This function samples a (tentative) free-flight distance according to an
     * exponential transmittance. It is then up to the integrator to then decide
     * whether the MediumInteraction corresponds to a real or null scattering
     * event.
     *
     * \param ray      Ray, along which a distance should be sampled
     * \param sample   A uniformly distributed random sample
     * \param channel  The channel according to which we will sample the
     * free-flight distance. This argument is only used when rendering in RGB
     * modes.
     *
     * \return         This method returns a MediumInteraction.
     *                 The MediumInteraction will always be valid,
     *                 except if the ray missed the Medium's bounding box.
     */
    MediumInteraction3f sample_interaction(const Ray3f &ray, Float sample,
                                           UInt32 channel, Mask active) const;

    /**
     * \brief Compute the transmittance and PDF
     *
     * This function evaluates the transmittance and PDF of sampling a certain
     * free-flight distance The returned PDF takes into account if a medium
     * interaction occurred (mi.t <= si.t) or the ray left the medium (mi.t >
     * si.t)
     *
     * The evaluated PDF is spectrally varying. This allows to account for the
     * fact that the free-flight distance sampling distribution can depend on
     * the wavelength.
     *
     * \return   This method returns a pair of (Transmittance, PDF).
     *
     */
    std::pair<UnpolarizedSpectrum, UnpolarizedSpectrum>
    transmittance_eval_pdf(const MediumInteraction3f &mi,
                           const SurfaceInteraction3f &si,
                           Mask active) const;

// #RAY_CHANGE_BEGIN, NM 05/06/2024 : add function that calculates the transmittance and pdf  
    virtual std::tuple<MediumInteraction3f, Float, Float> 
    sample_interaction_real(const Ray3f &ray, const SurfaceInteraction3f &si,
                            Float sample, UInt32 channel, Mask active) const;
    
    virtual std::tuple<Float, Float, Mask>
    eval_transmittance_pdf_real(const Ray3f &ray, 
                                    const SurfaceInteraction3f &si,
                                    UInt32 channel, Mask active) const;
// #RAY_CHANGE_END

    /// Return the phase function of this medium
    MI_INLINE const PhaseFunction *phase_function() const {
        return m_phase_function.get();
    }

    /// Returns whether this specific medium instance uses emitter sampling
    MI_INLINE bool use_emitter_sampling() const { return m_sample_emitters; }

    /// Returns whether this medium is homogeneous
    MI_INLINE bool is_homogeneous() const { return m_is_homogeneous; }

    /// Returns whether this medium has a spectrally varying extinction
    MI_INLINE bool has_spectral_extinction() const {
        return m_has_spectral_extinction;
    }

// #ERADIATE_CHANGE_BEGIN: Extremum structure accessor
    /// Returns the extremum structure for local majorant acceleration (nullptr if not used)
    // MI_INLINE const ExtremumStructure* extremum_structure() const {
    //     return m_extremum_structure.get();
    // }

    /// Check if medium uses extremum structure
    // MI_INLINE bool has_extremum_structure() const {
    //     return m_extremum_structure != nullptr;
    // }
// #ERADIATE_CHANGE_END

    void traverse(TraversalCallback *callback) override;

    /// Return a human-readable representation of the Medium
    std::string to_string() const override = 0;

    MI_DECLARE_PLUGIN_BASE_CLASS(Medium)

protected:
    Medium();
    Medium(const Properties &props);

protected:
    ref<PhaseFunction> m_phase_function;
    bool m_sample_emitters;
    bool m_is_homogeneous;
    bool m_has_spectral_extinction;
// #ERADIATE_CHANGE_BEGIN: Extremum structure support
    // ref<ExtremumStructure> m_extremum_structure;
// #ERADIATE_CHANGE_END

    // MI_DECLARE_TRAVERSE_CB(m_phase_function, m_extremum_structure)
    MI_DECLARE_TRAVERSE_CB(m_phase_function)
};

MI_EXTERN_CLASS(Medium)
NAMESPACE_END(mitsuba)

// -----------------------------------------------------------------------
//! @{ \name Enables vectorized method calls on Dr.Jit medium arrays
// -----------------------------------------------------------------------

DRJIT_CALL_TEMPLATE_BEGIN(mitsuba::Medium)
    DRJIT_CALL_GETTER(phase_function)
    DRJIT_CALL_GETTER(use_emitter_sampling)
    DRJIT_CALL_GETTER(is_homogeneous)
    DRJIT_CALL_GETTER(has_spectral_extinction)
    DRJIT_CALL_METHOD(get_majorant)
    DRJIT_CALL_METHOD(intersect_aabb)
    DRJIT_CALL_METHOD(sample_interaction)
    DRJIT_CALL_METHOD(transmittance_eval_pdf)
    DRJIT_CALL_METHOD(sample_interaction_real)
    DRJIT_CALL_METHOD(eval_transmittance_pdf_real)
    DRJIT_CALL_METHOD(get_scattering_coefficients)
DRJIT_CALL_END()

//! @}
// -----------------------------------------------------------------------
