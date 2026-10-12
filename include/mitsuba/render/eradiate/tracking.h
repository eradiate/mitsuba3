#pragma once

#include <mitsuba/render/fwd.h>

NAMESPACE_BEGIN(mitsuba)

/**
 * \brief Result of one tracking traversal, as returned by
 * ``Medium::delta_track`` and ``Medium::ratio_track``.
 *
 * Every event of the traversal (a free flight, a null or a real collision)
 * has a contribution \c f (a spectrum), a sampling density \c p (a spectrum:
 * the density with which each channel's majorants would have sampled it) and
 * a channel-independent probability \c P of the null / real decision. With
 * \c p_c the density of the hero channel, the record holds the products over
 * all events of <tt>f / (p_c P)</tt>, <tt>f / p_c</tt> and <tt>p / p_c</tt>.
 */
template <typename Float, typename Spectrum>
struct TrackingRecord {
    MI_IMPORT_TYPES()

    /// Interaction of the real collision. Valid only if one occurred.
    MediumInteraction3f mei;
    /// Medium component of the real collision.
    UInt32 component;
    /// Product of <tt>f / (p_c P)</tt>: the path throughput weight.
    UnpolarizedSpectrum throughput;
    /// Product of <tt>f / p_c</tt>: the transmittance estimate.
    UnpolarizedSpectrum transmittance;
    /// Product of <tt>p / p_c</tt>: the densities relative to the hero channel.
    UnpolarizedSpectrum pdf_ratio;

    DRJIT_STRUCT(TrackingRecord, mei, component, throughput, transmittance,
                 pdf_ratio)
};

/// Helper function to index the channel of an ``UnpolarizedSpectrum``.
template< typename Float, typename Spectrum >
MI_INLINE
Float index_spectrum(
    const unpolarized_spectrum_t<Spectrum> &spec,
    const dr::uint32_array_t<Float> &idx
) {
    Float m = spec[0];
    if constexpr (is_rgb_v<Spectrum>) { // Handle RGB rendering
        dr::masked(m, idx == 1u) = spec[1];
        dr::masked(m, idx == 2u) = spec[2];
    } else {
        DRJIT_MARK_USED(idx);
    }
    return m;
}

NAMESPACE_END(mitsuba)
