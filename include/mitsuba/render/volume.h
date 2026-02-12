#pragma once

#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/profiler.h>
#include <mitsuba/core/transform.h>
#include <mitsuba/render/interaction.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/render/texture.h>

#include <drjit/texture.h>

NAMESPACE_BEGIN(mitsuba)

/// Abstract base class for 3D volumes.
template <typename Float, typename Spectrum>
class MI_EXPORT_LIB Volume : public JitObject<Volume<Float, Spectrum>> {
public:
    MI_IMPORT_TYPES(Texture)

    // ======================================================================
    //! @{ \name Volume interface
    // ======================================================================

    /// Evaluate the volume at the given surface interaction, with color processing.
    virtual UnpolarizedSpectrum eval(const Interaction3f &it, Mask active = true) const;

    /// Evaluate this volume as a single-channel quantity.
    virtual Float eval_1(const Interaction3f &it, Mask active = true) const;

    /// Evaluate this volume as a three-channel quantity with no color processing (e.g. velocity field).
    virtual Vector3f eval_3(const Interaction3f &it, Mask active = true) const;

   /**
     * Evaluate this volume as a six-channel quantity with no color processing
     * This interface is specifically intended to encode the parameters of an SGGX phase function.
     */
    virtual dr::Array<Float, 6> eval_6(const Interaction3f &it, Mask active = true) const;

    /**
     * \brief Evaluate this volume as a n-channel float quantity
     *
     * This interface is specifically intended to encode a variable number of parameters.
     * Pointer allocation/deallocation must be performed by the caller.
     */
    virtual void eval_n(const Interaction3f &it, Float *out, Mask active = true) const;

    /**
     * Evaluate the volume at the given surface interaction,
     * and compute the gradients of the linear interpolant as well.
     */
    virtual std::pair<UnpolarizedSpectrum, Vector3f> eval_gradient(const Interaction3f &it,
                                                                   Mask active = true) const;

    /// Returns the maximum value of the volume over all dimensions.
    virtual ScalarFloat max() const;

    /**
     * \brief In the case of a multi-channel volume, this function returns
     * the maximum value for each channel.
     *
     * Pointer allocation/deallocation must be performed by the caller.
     */
    virtual void max_per_channel(ScalarFloat *out) const;

// #ERADIATE_CHANGE_BEGIN: Local extremum support
    /**
     * \brief Compute local extrema over a spatial region
     *
     * Returns the maximum value (majorant) and minimum value (minorant)
     * over the specified bounding box region.
     * Only fully implemented for grid-based volumes.
     *
     * \param bounds  Bounding box defining the query region in local space
     * \return (majorant, minorant) pair
     * 
     * The reference to array is a bit of a codesmell, used to avoid the ref-count
     * cost. 
     */
    virtual std::pair<Float, Float>
    extremum(const DynamicBuffer<Float>* array, BoundingBox3f local_bounds) const;

    virtual const ScalarFloat* data() const;
    virtual const DynamicBuffer<Float>* array() const;

    /**
     * \brief Compute local majorant (maximum) over a spatial region
     *
     * Convenience method that returns only the majorant.
     */
    virtual Float majorant(const BoundingBox3f &local_bounds) const {
        return extremum(nullptr, local_bounds).first;
    }

    /**
     * \brief Compute local minorant (minimum) over a spatial region
     *
     * Convenience method that returns only the minorant.
     */
    virtual Float minorant(const BoundingBox3f &local_bounds) const {
        return extremum(nullptr, local_bounds).second;
    }
// #ERADIATE_CHANGE_END

    /// Returns the bounding box of the volume
    ScalarBoundingBox3f bbox() const { return m_bbox; }

// #RAY_CHANGE_BEGIN, NM 24/05/2024 : Add util functions to the volume class
    /// If applicable, returns the dimensions of one grid cell in world space.
    ScalarVector3f voxel_size() const;
// #RAY_CHANGE_END

    /**
     * \brief Returns the resolution of the volume, assuming that it is based
     * on a discrete representation.
     *
     * The default implementation returns <tt>(1, 1, 1)</tt>
     */
    virtual ScalarVector3i resolution() const;

    /**
     * \brief Returns the number of channels stored in the volume
     *
     *  When the channel count is zero, it indicates that the volume
     *  does not support per-channel queries.
     */
    uint32_t channel_count() const { return m_channel_count; }

    //! @}
    // ======================================================================

    /// Returns a human-reable summary
    std::string to_string() const override {
        std::ostringstream oss;
        oss << "Volume[" << std::endl
            << "  to_local = " << m_to_local << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_PLUGIN_BASE_CLASS(Volume)

protected:
    Volume(const Properties &props);

    void update_bbox() {
        ScalarAffineTransform4f to_world = m_to_local.inverse();
        m_bbox = ScalarBoundingBox3f();
        m_bbox.expand(to_world * ScalarPoint3f(0.f, 0.f, 0.f));
        m_bbox.expand(to_world * ScalarPoint3f(0.f, 0.f, 1.f));
        m_bbox.expand(to_world * ScalarPoint3f(0.f, 1.f, 0.f));
        m_bbox.expand(to_world * ScalarPoint3f(0.f, 1.f, 1.f));
        m_bbox.expand(to_world * ScalarPoint3f(1.f, 0.f, 0.f));
        m_bbox.expand(to_world * ScalarPoint3f(1.f, 0.f, 1.f));
        m_bbox.expand(to_world * ScalarPoint3f(1.f, 1.f, 0.f));
        m_bbox.expand(to_world * ScalarPoint3f(1.f, 1.f, 1.f));
    }

protected:
    /// Used to bring points in world coordinates to local coordinates.
    ScalarAffineTransform4f m_to_local;
    /// Bounding box
    ScalarBoundingBox3f m_bbox;
    /// Number of channels stored in the volume
    uint32_t m_channel_count;

    MI_TRAVERSE_CB(Object)
};

MI_EXTERN_CLASS(Volume)
NAMESPACE_END(mitsuba)
