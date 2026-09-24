#pragma once

#include <mitsuba/core/fwd.h>
#include <mitsuba/core/vector.h>

#include <array>

NAMESPACE_BEGIN(mitsuba)

/**
 * \brief Traversal state carried between successive \c dda_next calls of an
 * \c ExtremumStructure.
 *
 * Holds what a structure cannot cheaply recover per segment: the
 * structure-local ray and the current cell. The cursor \c mint is the entry
 * distance of the next segment; traversal is over once <tt>mint >= maxt</tt>.
 *
 * The ray is kept as a bare origin and direction rather than a \c Ray3f, and
 * anything derivable from them in a few instructions is recomputed by
 * \c dda_next rather than carried, so that the struct stays small and
 * plugin-agnostic.
 */
template <typename Float, typename Spectrum>
struct DDAState {
    MI_IMPORT_CORE_TYPES()

    /// Ray origin in structure-local coordinates.
    Point3f o;
    /// Ray direction in structure-local coordinates. Never renormalised, so
    /// the world ray's \c t parameterization carries over unchanged.
    Vector3f d;
    /// Current cell index. Sentinel values (-1, resolution) mark the exterior.
    Vector3i pi;
    /// Cursor: entry distance of the next segment. Exact in \c t.
    Float mint;
    /// End of the traversal range.
    Float maxt;

    DRJIT_STRUCT(DDAState, o, d, pi, mint, maxt)
};

/// Maximum number of extremum structures a \ref DDAStateList carries. Every
/// entry costs a whole \ref DDAState of loop state.
static constexpr size_t MAX_DDA_OVERLAP = 4;

/**
 * \brief Traversal state of a medium's extremum structures, one
 * \ref DDAState per structure.
 *
 * The entries advance independently: an entry is asked for a new segment only
 * once the shared cursor \c mint reaches the end of the one it last returned.
 * The medium's segment is rebuilt from the cached bounds on every step.
 *
 * Entries past the medium's structure count are never touched.
 */
template <typename Float, typename Spectrum>
struct DDAStateList {
    MI_IMPORT_CORE_TYPES()

    /// Per-structure traversal state. <tt>entries[i].mint</tt> doubles as the
    /// exit distance of entry \c i's cached segment.
    std::array<DDAState<Float, Spectrum>, MAX_DDA_OVERLAP> entries;
    /// Per-structure <tt>[minorant, majorant]</tt> of that cached segment,
    /// zero outside the structure's domain.
    std::array<Vector2f, MAX_DDA_OVERLAP> values;
    /// Cursor: entry distance of the next segment.
    Float mint;
    /// End of the traversal range.
    Float maxt;

    DRJIT_STRUCT(DDAStateList, entries, values, mint, maxt)
};

NAMESPACE_END(mitsuba)
