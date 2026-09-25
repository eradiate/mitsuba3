#pragma once

#include <mitsuba/core/object.h>
#include <mitsuba/render/interaction.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/eradiate/extremum_segment.h>
#include <mitsuba/render/eradiate/dda.h>

#include <optional>

NAMESPACE_BEGIN(mitsuba)

/**
 * \brief Abstract base class for extremum structures
 *
 * ExtremumStructure provides an interface for spatial data structures that
 * store local extrema (majorant/minorant) of volumetric extinction coefficients.
 * This enables efficient use of tracking algorithms with locally-adaptive
 * majorants and minorants.
 *
 * Structures are traversed through ``dda_init`` / ``dda_next``, driven by
 * the owning ``Medium``. They are host-side objects: no method is reachable
 * through a Dr.Jit vcall.
 *
 * The extremum structure needs to be built using the ``update_extremum``
 * function, it is **not** called automatically in the constructor. The caller,
 * usually a ``Medium``, passes the volumetric data the extremum is derived from.
 */
template <typename Float, typename Spectrum>
class MI_EXPORT_LIB ExtremumStructure : public JitObject<ExtremumStructure<Float, Spectrum>> {
public:
    MI_IMPORT_TYPES(Medium, Sampler, Volume)

    /// Destructor
    ~ExtremumStructure();

    /// Setter for the bbox over which the structure must be valid.
    MI_INLINE void set_bbox(ScalarBoundingBox3f bbox) { m_bbox = bbox; };

    /// Setter for the scale by which to multiply the extremum values.
    MI_INLINE void set_scale(ScalarFloat scale) { m_scale = scale; }

    /**
     * \brief Update the bbox and scale, and rebuild the structure.
     *
     * The \c bbox parameters indicates the domain over which the extremum
     * structure can be queried. It can be larger or smaller than the underlying
     * volume bbox. It is the extremum's responsibility to be valid over this
     * area. The building implementation is handled in ``build``.
     *
     * \param bbox      The validity bbox of the extremum structure
     * \param volume    The volume from which to derive the extremum structure
     * \param scale     The scale by which to multiply the extremum values
     */
    void update_extremum(const ScalarBoundingBox3f &bbox,
                         const Volume *volume,
                         std::optional<ScalarFloat> scale);

    /**
     * \brief Build the extremum structure of \c volume.
     *
     * Implements the logic that constructs the extremum structure from a
     * \c volume. Called by ``update_extremum`` which is itself called by
     * the owning ``Medium``
     *
     * \param volume  Volume to compute extremum values from
     */
    virtual void build(const Volume *volume) = 0;

    /**
     * \brief Set up a stateful DDA traversal along \c ray.
     *
     * Transforms the ray to structure-local coordinates, locates the entry
     * cell and clips <tt>[mint, maxt]</tt> to the structure's domain. Called
     * once per ray; the returned state is then advanced with ``dda_next``
     * until <tt>state.mint >= state.maxt</tt>.
     *
     * Host-only: reached through a \c Medium, never through a Dr.Jit vcall.
     */
    virtual DDAState dda_init(const Ray3f &ray, Float mint, Float maxt,
                              Mask active = true) const = 0;

    /**
     * \brief Return the segment starting at <tt>state.mint</tt>, together
     * with the state advanced past it.
     *
     * Segments are half-open and tile exactly: <tt>segment.mint ==
     * state.mint</tt> and <tt>segment.maxt == next_state.mint</tt>, exact in
     * \c t.
     *
     * Host-only, like ``dda_init``.
     */
    virtual std::pair<ExtremumSegment, DDAState>
    dda_next(const DDAState &state, Mask active = true) const = 0;

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

    ScalarBoundingBox3f bbox() const { return m_bbox; }
    //! @}
    // =============================================================

    MI_DECLARE_PLUGIN_BASE_CLASS(ExtremumStructure)

protected:
    ExtremumStructure();
    ExtremumStructure(const Properties &props);

protected:
    /// The bbox over which the extremum structure must be valid.
    ScalarBoundingBox3f m_bbox;
    /// Scale by which to multiply the extremum values.
    ScalarFloat m_scale;
};

MI_EXTERN_CLASS(ExtremumStructure)
NAMESPACE_END(mitsuba)

