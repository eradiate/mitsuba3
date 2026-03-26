#pragma once

#include <mitsuba/core/fwd.h>

NAMESPACE_BEGIN(mitsuba)

template<typename Float, typename Spectrum>
struct ExtremumSegment {
    MI_IMPORT_CORE_TYPES()                                                                         \

    /// Segment entry distance along ray
    Float mint;
    /// Segment exit distance along ray
    Float maxt;
    /// Extremum data stored as [minorant, majorant] 
    Vector2f value;

    ExtremumSegment(){ reset(); };

    ExtremumSegment(
        Float mint, 
        Float maxt, 
        Vector2f value
    ) : mint(mint), 
        maxt(maxt), 
        value(value) {}

    ExtremumSegment(
        const Float& mint, 
        const Float& maxt,
        const Float& minorant,
        const Float& majorant
    ) : mint(mint), 
        maxt(maxt), 
        value(Vector2f(minorant, majorant)) {}

    /**
     * This callback method is invoked by dr::zeros<>, and takes care of fields 
     * that deviate from the standard zero-initialization convention. In 
     * ExtremumSegment, the ``mint`` and ``maxt`` fields are set to  + and - 
     * infinity respectively to to mark invalid intersection records.
     */
    void zero_(size_t size = 1) {                                                                                                                                                                            
        mint        = dr::full<Float>(dr::Infinity<Float>, size);
        maxt        = dr::full<Float>(-dr::Infinity<Float>, size);
        value       = dr::zeros<Vector2f>(size);
    }  

    /**
     * \brief Check whether this is a valid segment
     * 
     * A segment is considered valid when 
     * \code
     * segment.mint < segment.maxt
     * \endcode
     */
    Mask valid() const {
        return mint < maxt;
    }

    /**
     * \brief Mark the extremum segment as invalid.
     *
     * This operation sets segment's minimum
     * and maximum distances to \f$\infty\f$ and \f$-\infty\f$,
     * respectively.
     */
    void reset() {
        mint = dr::Infinity<Float>;
        maxt = -dr::Infinity<Float>;
    }

    /// Minorant value over the segment. Accessor to the first element of ``value``.
    Float minorant() {
        return value.x();
    }

    /// Majorant value over the segment. Accessor to the second element of ``value``.
    Float majorant() {
        return value.y();
    }

    DRJIT_STRUCT_NODEF(ExtremumSegment, mint, maxt, value)
};

NAMESPACE_END(mitsuba)