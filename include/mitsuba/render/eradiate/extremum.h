#pragma once

#include <mitsuba/core/object.h>
#include <mitsuba/render/interaction.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/eradiate/extremum_segment.h>
#include <mitsuba/render/eradiate/tracking.h>
#include <drjit/call.h>

NAMESPACE_BEGIN(mitsuba)

/**
 * \brief Abstract base class for extremum structures
 *
 * ExtremumStructure provides an interface for spatial data structures that
 * store local extrema (majorant/minorant) of volumetric values.
 * This enables efficient delta tracking with locally-adaptive majorants.
 *
 * To minimize virtual function overhead, the ``traverse_extremum()`` method
 * encapsulates the entire traversal loop internally, requiring only a single
 * virtual call per distance sample.
 */
template <typename Float, typename Spectrum>
class MI_EXPORT_LIB ExtremumStructure : public JitObject<ExtremumStructure<Float, Spectrum>> {
public:
    MI_IMPORT_TYPES(Medium, Sampler, Volume)

    using TrackingStateType    = TrackingState<Float, Spectrum>;
    using TrackingFunctionType = TrackingFunction<Float, Spectrum>;

    /// Destructor
    ~ExtremumStructure();

    /**
     * \brief Build the extremum structure of \c sigma_t over \c domain.
     *
     * Called by the owning medium at construction and parameter update.
     * The default implementation throws.
     *
     * \param domain  Support bound of the medium's extinction in world
     *                space; may be larger than the volume's bounding box
     *                (infinite extents allowed)
     * \param volume  Volume to compute bounds from
     * \param scale   Scale factor applied to the volume values
     */
    virtual void build(const ScalarBoundingBox3f &domain,
                       const Volume *volume, ScalarFloat scale);

    /**
     * \brief Return the segment <tt>[t, t_exit)</tt> containing distance
     * \c t along \c ray, with the structure's local extrema.
     *
     * Stateless query in world-t parameterization. Conventions:
     *
     * - Segments are half-open <tt>[mint, maxt)</tt> and tile exactly:
     *   feeding \c maxt back as the next \c t yields the adjacent segment.
     * - Outside the domain, empty segments are returned — <tt>[t, enter)</tt>
     *   or <tt>[t, +inf)</tt> with value (0, 0): outside the domain the
     *   value is zero by construction.
     */
    virtual ExtremumSegment next_segment(const Ray3f &ray, Float t,
                                         Mask active = true) const;

    /**
     * \brief Traverse the extremum along a ray and applies a callback at each
     * encountered segment.
     *
     * This method traverses the extremum structure segment by segment. At each
     * segment, the callback ``func`` is called to advance the ``state``. This
     * is useful for example to implement Delta Tracking, Ratio Tracking, and
     * Residual Ratio Tracking. The callback is typically defined in the
     * integrator. By default, the traversal is achieved using ``next_segment``.
     * This function can be overriden to achieve better traversal performance.
     *
     * \param ray           Ray along which to sample
     * \param mint          Minimum distance to consider
     * \param maxt          Maximum distance to consider
     * \param channel       Channel from which to sample
     * \param state         Mutable tracking state carried through the traversal loop
     * \param func          Callback function called at every segment.
     * \param active        Mask for active lanes
     *
     * \return
     *      The final tracking state, that includes the medium interaction if
     *      a real scattering event was sampled, and the throughput and pdfs
     *      accumulated throughout the traversal.
     *
     * Note that this function cannot be made abstract because of it would
     * force the requirement for bindings, which are incompatible with
     * function types.
     */
    virtual TrackingStateType traverse_extremum(
        const Ray3f &ray,
        Float mint,
        Float maxt,
        UInt32 channel,
        TrackingStateType state,
        TrackingFunctionType *func,
        Mask active = true
    ) const;

    // Note: this is currently dead code. It is kept in case it is needed in the future.
    /**
     * \brief Evaluate the minorant and majorant at a medium interaction point.
     *
     * This method performs point evaluation at interaction point specified in
     * local space.
     *
     * \param it            Interaction interaction point in local space
     * \param active        Mask for active lanes
     *
     * \return
     *      The minorant and majorant values at the medium interaction point.
     *      Clamped values outside bounds.
     *
     */
    virtual std::tuple<Float, Float> eval_1(
        const Interaction3f & it,
        Mask active = true
    ) const = 0;

    // =============================================================
    //! @{ \name Non-virtual query methods
    // =============================================================

    /// Return the bounding box of the extremum structure
    ScalarBoundingBox3f bbox() const { return m_bbox; }
    //! @}
    // =============================================================

    MI_DECLARE_PLUGIN_BASE_CLASS(ExtremumStructure)

protected:
    ExtremumStructure();
    ExtremumStructure(const Properties &props);

    /**
     * \brief Record the build inputs; to be called by build() implementations.
     *
     * Sets the structure's bounding box to the domain and throws if the
     * structure was already built for a different (domain, volume) pair —
     * one extremum structure belongs to exactly one medium.
     */
    void set_domain(const ScalarBoundingBox3f &domain, const Volume *volume);

protected:
    /// Bounding box of the extremum structure in world space (== the domain)
    ScalarBoundingBox3f m_bbox;

    /// Identity of the volume this structure was built for. Not dereferenced,
    /// only used for identity check.
    const Volume *m_built_volume = nullptr;
};

MI_EXTERN_CLASS(ExtremumStructure)
NAMESPACE_END(mitsuba)

// -----------------------------------------------------------------------
//! @{ \name Enables vectorized method calls on Dr.Jit medium arrays
// -----------------------------------------------------------------------

DRJIT_CALL_TEMPLATE_BEGIN(mitsuba::ExtremumStructure)
    DRJIT_CALL_METHOD(traverse_extremum)
    DRJIT_CALL_METHOD(next_segment)
    DRJIT_CALL_METHOD(eval_1)
DRJIT_CALL_END()

//! @}
// -----------------------------------------------------------------------
