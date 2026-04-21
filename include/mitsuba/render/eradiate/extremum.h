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
 * store local extrema (majorant/minorant) of volumetric extinction coefficients.
 * This enables efficient delta tracking with locally-adaptive majorants.
 *
 * To minimize virtual function overhead, the ``traverse_extremum()`` method 
 * encapsulates the entire traversal loop internally, requiring only a single
 * virtual call per distance sample.
 */
template <typename Float, typename Spectrum>
class MI_EXPORT_LIB ExtremumStructure : public JitObject<ExtremumStructure<Float, Spectrum>> {
public:
    MI_IMPORT_TYPES(Medium, Sampler)

    using TrackingState    = TrackingState<Float, Spectrum>;
    using TrackingFunction = TrackingFunction<Float, Spectrum>;

    /// Destructor
    ~ExtremumStructure();

    /**
     * \brief Traverse the extremum along a ray and applies a callback at each
     * encountered segment.         .
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
    virtual TrackingState traverse_extremum(
        const Ray3f &ray,
        Float mint,
        Float maxt,
        UInt32 channel,
        TrackingState state,
        TrackingFunction *func,
        Mask active = true
    ) const;

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
     * Note: this is currently dead code. It is kept in case it is needed in the
     * future.
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

protected:
    /// Bounding box of the extremum structure in world space
    ScalarBoundingBox3f m_bbox;
};

MI_EXTERN_CLASS(ExtremumStructure)
NAMESPACE_END(mitsuba)

// -----------------------------------------------------------------------
//! @{ \name Enables vectorized method calls on Dr.Jit medium arrays
// -----------------------------------------------------------------------

DRJIT_CALL_TEMPLATE_BEGIN(mitsuba::ExtremumStructure)
    DRJIT_CALL_METHOD(traverse_extremum)
    DRJIT_CALL_METHOD(eval_1)
DRJIT_CALL_END()

//! @}
// -----------------------------------------------------------------------
