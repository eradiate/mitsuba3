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
    MI_IMPORT_TYPES(PhaseFunction, Sampler, Scene, Texture, ExtremumStructure);

    /// Destructor
    ~Medium();

    /// Intersects a ray with the medium's bounding box
    virtual std::tuple<Mask, Float, Float>
    intersect_aabb(const Ray3f &ray) const = 0;
// #ERADIATE_CHANGE_BEGIN: Overlapping media

    /// Checks if a point is contained by the medium's bounding box
    virtual Mask
    in_aabb(const Point3f &pos) const = 0;
// #ERADIATE_CHANGE_END

    /// Returns the medium's majorant used for delta tracking
    virtual UnpolarizedSpectrum
    get_majorant(const MediumInteraction3f &mi,
                 Mask active = true) const = 0;

    /// Returns the medium's minorant used for residual ratio tracking
     virtual UnpolarizedSpectrum
     get_minorant(const MediumInteraction3f &/*mi*/,
                  Mask /*active*/ = true) const { return UnpolarizedSpectrum(0.f); }

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

// #ERADIATE_CHANGE_BEGIN, NM 05/06/2024 : add function that calculates the transmittance and pdf
    /**
     * \brief Sample a free-flight distance in the medium analytically.
     *
     * This function samples a (tentative) free-flight distance according to an
     * exponential transmittance. It is then up to the integrator to then decide
     * whether the MediumInteraction corresponds to a real or null scattering
     * event.
     *
     * \param ray      Ray, along which a distance should be sampled
     * \param it       The boundary interaction that the sampled distance cannot
     * exceed.
     * \param sample   A uniformly distributed random sample
     * \param channel  The channel according to which we will sample the
     * free-flight distance. This argument is only used when rendering in RGB
     * modes.
     *
     * \return         This method returns a tuple
     *                 (MediumInteraction, Transmittance, PDF).
     *                 The MediumInteraction if an interaction was sampled within
     *                 the medium boudning box and before the bouding iteraction
     *                 \ref it.
     *                 The transmittance and PDF are both computed for all channels
     *                 even if the sampling operation is performed on one channel.
     *
     */
    virtual
    std::tuple<MediumInteraction3f, UnpolarizedSpectrum, UnpolarizedSpectrum>
    sample_interaction_analytical(const Ray3f &ray, const Interaction3f &it,
                            Float sample, UInt32 channel, Mask active) const;

    /** \brief Compute the analytical transmittance along a ray to an interaction.
     *
     * \param ray      Ray, along which to compute the transmittance, use mint
     * \param si       Interaction that marks the end of the segment along which
     *                 to compute the transmittance.
     *
     * \return
     *      The transmittance along a ray
    */
    virtual UnpolarizedSpectrum
    transmittance_eval_analytical(const Ray3f &ray, const Interaction3f &it,
                                  Mask active) const;
// #ERADIATE_CHANGE_END

    /// Return the phase function of this medium
    MI_INLINE const PhaseFunction *phase_function() const {
        return m_phase_function.get();
    }
// #ERADIATE_CHANGE_BEGIN: DDIS
    /// Return the phase function of this medium, non const
    MI_INLINE PhaseFunction *phase_function() {
        return m_phase_function.get();
    }
// #ERADIATE_CHANGE_END

    /// Returns whether this specific medium instance uses emitter sampling
    MI_INLINE bool use_emitter_sampling() const { return m_sample_emitters; }

    /// Returns whether this medium is homogeneous
    MI_INLINE bool is_homogeneous() const { return m_is_homogeneous; }

    /// Returns whether this medium has a spectrally varying extinction
    MI_INLINE bool has_spectral_extinction() const {
        return m_has_spectral_extinction;
    }

// #ERADIATE_CHANGE_BEGIN: Extremum Structure accessor and traversal helpers
    /**
     * \brief Intersects ray with the medium bbox and creates a medium interaction.
     *
     * \param ray   The ray that is used to test the medium bbox.
     *
     * \return
     *      A tuple (mei, mint, maxt): ``mei`` is a  ``MediumInteraction3f``
     *      object initialized with the current ray and medium data. ``mint``
     *      and ``maxt`` represent the minimum and maximum intersection
     *      distances of the ray with the medium's bbox. In case there are no
     *      valid intersection, the range defaults to [0, +Inf].
     */
    std::tuple<MediumInteraction3f, Float, Float>
    prepare_medium_traversal(const Ray3f& ray, Mask active) const;

    /// Returns the extremum structure for local extremum acceleration.
    MI_INLINE const ExtremumStructure *extremum_structure() const {
        return m_extremum_structure.get();
    }

    /// Check if medium uses extremum structure
    MI_INLINE bool has_extremum_structure() const {
        return m_extremum_structure.get() != nullptr;
    }

    MI_INLINE bool use_rrt() const {
        return m_use_rrt;
    }
// #ERADIATE_CHANGE_END

// #ERADIATE_CHANGE_BEGIN: DDIS
    /**
     * \brief Return the ddis phase function of this medium. Can be null for
     * medium that don't specify such phase function.
     */
    MI_INLINE const PhaseFunction *ddis_phase_function() const {
        return m_ddis_phase_function.get();
    }

    MI_INLINE ScalarFloat ddis_threshold() const {
        return m_ddis_threshold;
    }


    /**
     * \brief Rebuild any DDIS-related derived data when the phase function has
     * changed.
     *
     * This method is called by the scene's `parameters_changed()` for every
     * medium whose phase function is marked dirty. The default implementation
     * is a no-op; subclasses that maintain a DDIS phase function override it.
     *
     * This is intentionally separate from `parameters_changed()` so that the
     * scene can drive all media to completion before clearing dirty flags,
     * which is required when multiple media share the same phase function via
     * a scene-level reference.
     */
    virtual void update_ddis_phase_function();

protected:
    /**
     * \brief Create the tabphase irregular plugin used as DDIS phase function.
     *
     * This method is an helper function for child classes.
     */
    virtual ref<PhaseFunction> create_ddis_phase_function();

public:

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
    ref<ExtremumStructure> m_extremum_structure;
    bool m_use_rrt;
// #ERADIATE_CHANGE_END

// #ERADIATE_CHANGE_BEGIN: DDIS
    ref<PhaseFunction> m_ddis_phase_function = nullptr;
    ScalarFloat m_ddis_threshold = 0.f;

    MI_DECLARE_TRAVERSE_CB(m_phase_function, m_extremum_structure,
                           m_ddis_phase_function, m_ddis_threshold)
// #ERADIATE_CHANGE_END
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
// #ERADIATE_CHANGE_BEGIN: Overlapping media
    DRJIT_CALL_METHOD(in_aabb)
// #ERADIATE_CHANGE_END
    DRJIT_CALL_METHOD(sample_interaction)
    DRJIT_CALL_METHOD(transmittance_eval_pdf)
// #ERADIATE_CHANGE_BEGIN: Add function that calculates the transmittance and pdf
    DRJIT_CALL_METHOD(sample_interaction_analytical)
    DRJIT_CALL_METHOD(transmittance_eval_analytical)
// #ERADIATE_CHANGE_END
    DRJIT_CALL_METHOD(get_scattering_coefficients)
// #ERADIATE_CHANGE_BEGIN: Extremum Support && Residual Ratio Tracking && DDIS
    DRJIT_CALL_GETTER(ddis_phase_function)
    DRJIT_CALL_GETTER(ddis_threshold)
    DRJIT_CALL_GETTER(use_rrt)
    DRJIT_CALL_GETTER(extremum_structure)
    DRJIT_CALL_METHOD(prepare_medium_traversal)
// #ERADIATE_CHANGE_END
DRJIT_CALL_END()

//! @}
// -----------------------------------------------------------------------
