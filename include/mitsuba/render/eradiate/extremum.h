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
template <typename Float, typename Spectrum>
struct ExtremumSegment {
    using Mask = dr::mask_t<Float>;

    /// Segment entry distance along ray
    Float tmin;
    /// Segment exit distance along ray
    Float tmax;
    /// Local majorant (maximum extinction) in this segment
    Float majorant;
    /// Local minorant (minimum extinction) in this segment
    Float minorant;

    /** 
     * \brief Create a new invalid extremum segment 
     * 
     * Initializes the minimum and maximum segment distances to \f$\infty\f$ 
     * and \f$-\infty\f$, respectively.
     */
    ExtremumSegment() { reset(); }

    /// Create an extremum segment from its fields.
    ExtremumSegment(
        Float tmin, 
        Float tmax, 
        Float majorant, 
        Float minorant
    ) : tmin(tmin), 
        tmax(tmax), 
        majorant(majorant), 
        minorant(minorant) {  }

    /**
     * This callback method is invoked by dr::zeros<>, and takes care of fields that deviate
     * from the standard zero-initialization convention. In this particular class, 
     * the ``tmin`` and ``tmax`` fields should be set to + and - infinity respectively to 
     * to mark invalid intersection records.
     */
    void zero_(size_t size = 1) {                                                                                                                                                                            
        tmin        = dr::full<Float>(dr::Infinity<Float>, size);
        tmax        = dr::full<Float>(-dr::Infinity<Float>, size);
        minorant   = dr::zeros<Float>(size);
        majorant   = dr::zeros<Float>(size);
    }  

    /**
     * \brief Check whether this is a valid segment
     * 
     * A segment is considered valid when 
     * \code
     * segment.tmin < segment.tmax
     * \endcode
     */
    Mask valid() const {
        return tmin <= tmax;
    }

    /**
     * \brief Mark the extremum segment as invalid.
     *
     * This operation sets segment's minimum
     * and maximum distances to \f$\infty\f$ and \f$-\infty\f$,
     * respectively.
     */
    void reset(){
        tmin = dr::Infinity<Float>;
        tmax = -dr::Infinity<Float>;
    }

    DRJIT_STRUCT_NODEF(ExtremumSegment, tmin, tmax, majorant, minorant)
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
     * \param target_od     Target optical depth to accumulate
     * \param active        Mask for active lanes
     *
     * \return 
     *      ExtremumSegment containing the sampled distance (tmin), segment
     *      bounds, and local majorant/minorant values. If desired_tau cannot
     *      be reached, tmin is set to Infinity.
     *      Accumulated optical thickness at segment start.
     */
    virtual std::tuple<ExtremumSegment, Float> sample_segment(
        const Ray3f &ray,
        Float mint,
        Float maxt,
        Float target_od,
        Mask active
    ) const = 0;

    /**
     * \brief Evaluate the majorant at a medium interaction point.
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
    DRJIT_CALL_METHOD(eval_1)
DRJIT_CALL_END()

//! @}
// -----------------------------------------------------------------------