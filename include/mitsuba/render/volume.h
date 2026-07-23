#pragma once

#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/profiler.h>
#include <mitsuba/core/transform.h>
#include <mitsuba/render/interaction.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/render/texture.h>

// #ERADIATE_CHANGE_BEGIN: Local extremum support
#include <mitsuba/core/platform.h>
#include <mitsuba/render/eradiate/volume_utils.h>
// #ERADIATE_CHANGE_END

#include <drjit/texture.h>

NAMESPACE_BEGIN(mitsuba)

/// Abstract base class for 3D volumes.
template <typename Float, typename Spectrum>
class MI_EXPORT_LIB Volume : public JitObject<Volume<Float, Spectrum>> {
public:
    MI_IMPORT_TYPES(Texture, ExtremumStructure)

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

// #ERADIATE_CHANGE_BEGIN: Tracking estimators extension
    /// Returns the minimum value of the volume over all dimensions.
    virtual ScalarFloat min() const;

    /**
     * \brief In the case of a multi-channel volume, this function returns
     * the minimum value for each channel.
     *
     * Pointer allocation/deallocation must be performed by the caller.
     */
    virtual void min_per_channel(ScalarFloat *out) const;
// #ERADIATE_CHANGE_END

// #ERADIATE_CHANGE_BEGIN: Local extremum support
    /**
     * \brief Bounds of \ref eval over a world-space box.
     *
     * Valid for any consumer and any volume: the returned bounds are always
     * correct, possibly loose. The default implementation returns
     * <tt>(0, max())</tt>.
     *
     * \param bbox  Query region in world space
     * \return (minorant, majorant) pair
     */
    virtual std::pair<Float, Float> extremum(const BoundingBox3f &bbox) const;

    /**
     * \brief Bounds of \ref eval over a box in the volume's *native*
     * parameterization, normalized to [0,1]^3.
     *
     * For Cartesian grids the axes are local grid coordinates; for
     * spherical-coordinates volumes they are (r, theta, phi). Callers must
     * first verify the parameterization via the matching frame capability
     * (e.g. \ref spherical_frame()). The default implementation assumes local
     * Cartesian axes and delegates to the world-space query.
     *
     * \param bbox  Query region in the native parameterization, in [0,1]^3
     * \return (minorant, majorant) pair
     */
    virtual std::pair<Float, Float> extremum_local(const BoundingBox3f &bbox) const;

    /**
     * \brief Capability gate for radially parameterized volumes.
     *
     * Returns the spherical frame if this volume's native parameterization
     * is (r, theta, phi). By default, returns an invalid frame.
     */
    virtual SphericalParameters<ScalarFloat> spherical_frame() const {
        return SphericalParameters<ScalarFloat>();
    };

    MI_INLINE bool dirty() const {
        return m_dirty;
    }

    MI_INLINE void set_dirty(bool dirty) {
        m_dirty = dirty;
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

// #ERADIATE_CHANGE_BEGIN: Local extremum support
    /**
     * \brief A Scoped Guard that pins the reference count of the volume.
     *
     * Use for bulk operations in scalar mode.
     */
    struct PinGuard {
        const Volume *volume;
        explicit PinGuard(const Volume* v) : volume(v) { volume->pin_ref_count(); }
        ~PinGuard() { volume->unpin_ref_count(); }

        PinGuard(const PinGuard&) = delete;
        PinGuard& operator=(const PinGuard&) = delete;
    };

    virtual PinGuard pin() const { return PinGuard(this); };
// #ERADIATE_CHANGE_END

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

// #ERADIATE_CHANGE_BEGIN: Local extremum support
    /// Pin the reference count of the data that constitutes the volume e.g. Texture/
    virtual void pin_ref_count() const {};
    /// Unpin the reference count.
    virtual void unpin_ref_count() const {};
// #ERADIATE_CHANGE_END

protected:
    /// Used to bring points in world coordinates to local coordinates.
    ScalarAffineTransform4f m_to_local;
    /// Bounding box
    ScalarBoundingBox3f m_bbox;
    /// Number of channels stored in the volume
    uint32_t m_channel_count;

// #ERADIATE_CHANGE_BEGIN: Local extremum support
    bool m_dirty;
// #ERADIATE_CHANGE_END

    MI_TRAVERSE_CB(Object)
};

MI_EXTERN_CLASS(Volume)
NAMESPACE_END(mitsuba)
