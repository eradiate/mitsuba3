#pragma once

#include <mitsuba/core/object.h>
#include <mitsuba/render/interaction.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/eradiate/extremum_segment.h>
#include <drjit/call.h>

NAMESPACE_BEGIN(mitsuba)

/**
 * \brief Abstract base class for extremum structures
 *
 * ExtremumStructure provides an interface for spatial data structures that
 * store local extrema (majorant/minorant) of volumetric extinction coefficients.
 * This enables efficient delta tracking with locally-adaptive majorants.
 *
 * To minimize virtual function overhead, the ``sample_segment()`` method 
 * encapsulates the entire traversal loop internally, requiring only a single
 * virtual call per distance sample.
 */
template <typename Float, typename Spectrum>
class MI_EXPORT_LIB ExtremumStructure : public JitObject<ExtremumStructure<Float, Spectrum>> {
public:
    MI_IMPORT_TYPES(Medium, Sampler)

    /// Destructor
    ~ExtremumStructure();

    /**
     * \brief Sample a segment along a ray with desired optical thickness
     *
     * This method traverses the extremum structure (e.g., via DDA for grids)
     * and returns a segment where the accumulated optical thickness reaches the
     * desired value. The traversal logic is completely encapsulated within
     * this method to minimize virtual call overhead.
     *
     * \param ray           Ray along which to sample
     * \param mint          Minimum distance to consider
     * \param maxt          Maximum distance to consider
     * \param target_ot     Target optical thickness to accumulate
     * \param active        Mask for active lanes
     *
     * \return 
     *      ExtremumSegment containing the sampled distance (mint), segment
     *      bounds, and local majorant/minorant values. If desired_tau cannot
     *      be reached, mint is set to Infinity.
     *      Accumulated optical thickness at segment start.
     */
    virtual std::tuple<ExtremumSegment, Float> sample_segment(
        const Ray3f &ray,
        Float mint,
        Float maxt,
        Float target_ot,
        Mask active = true
    ) const = 0;

    /**
     * \brief Evaluate the minorant and majorant at a medium interaction point.
     *
     * This method performs point evaluation at interaction point specified in 
     * local space.
     *
     * \param it            Interaction interaction point in local space
     * \param active        Mask for active lanes
     *
     * \return The minorant and majorant values at the medium interaction point.
     *         Clamped values outside bounds.
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
    DRJIT_CALL_METHOD(sample_segment)
    DRJIT_CALL_METHOD(eval_1)
DRJIT_CALL_END()

//! @}
// -----------------------------------------------------------------------
