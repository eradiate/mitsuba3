#pragma once

#include <mitsuba/core/object.h>
#include <mitsuba/render/interaction.h>
#include <mitsuba/render/volume.h>
#include <drjit/call.h>

NAMESPACE_BEGIN(mitsuba)

/**
 * \brief Segment along a ray with local extremum values
 *
 * This structure stores the entry/exit distances of a segment along a ray,
 * along with the local majorant and minorant within that segment.
 */
template <typename Float>
struct ExtremumSegment {
    /// Segment entry distance along ray
    Float tmin;
    /// Segment exit distance along ray
    Float tmax;
    /// Local majorant (maximum extinction) in this segment
    Float sigma_maj;
    /// Local minorant (minimum extinction) in this segment
    Float sigma_min;

    DRJIT_STRUCT(ExtremumSegment, tmin, tmax, sigma_maj, sigma_min)
};

/**
 * \brief Abstract base class for extremum structures
 *
 * ExtremumStructure provides an interface for spatial data structures that
 * store local extrema (majorant/minorant) of volumetric extinction coefficients.
 * This enables efficient delta tracking with locally-adaptive majorants.
 *
 * To minimize virtual function overhead, the `sample_segment()` method 
 * encapsulates the entire traversal loop internally, requiring only a single
 * virtual call per distance sample.
 */
template <typename Float, typename Spectrum>
class MI_EXPORT_LIB ExtremumStructure : public JitObject<ExtremumStructure<Float, Spectrum>> {
public:
    MI_IMPORT_TYPES()
    using ExtremumSegmentType = ExtremumSegment<Float>;

    /// Destructor
    ~ExtremumStructure();

    /**
     * \brief Sample a segment along a ray with desired optical depth
     *
     * This method traverses the extremum structure (e.g., via DDA for grids)
     * and returns a segment where the accumulated optical depth reaches the
     * desired value. The traversal logic is completely encapsulated within
     * this method to minimize virtual call overhead.
     *
     * \param ray           Ray along which to sample
     * \param mint          Minimum distance to consider
     * \param maxt          Maximum distance to consider
     * \param desired_tau   Target optical depth to accumulate
     * \param active        Mask for active lanes
     *
     * \return ExtremumSegment containing the sampled distance (tmin), segment
     *         bounds, and local majorant/minorant values. If desired_tau cannot
     *         be reached, tmin is set to Infinity.
     */
    virtual ExtremumSegmentType sample_segment(
        const Ray3f &ray,
        Float mint,
        Float maxt,
        Float desired_tau,
        Mask active
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
DRJIT_CALL_END()

//! @}
// -----------------------------------------------------------------------