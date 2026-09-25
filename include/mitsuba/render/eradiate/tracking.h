#pragma once

#include <mitsuba/render/fwd.h>
#include <mitsuba/render/eradiate/extremum_segment.h>
#include <drjit/random.h>

NAMESPACE_BEGIN(mitsuba)

/**
 * \brief State carried through extremum traversal to accumulate the throughput
 * and its PDF.
 *
 * Can be used to for delta tracking, ratio tracking, and residual ratio tracking.
 * Note: Since the number of required dimensions is different for all pixel
 * samples, ``rng`` is used to sample distances and event types.
 *
 */
template< typename Float, typename Spectrum >
struct TrackingState {
    MI_IMPORT_TYPES()

    Ray3f ray;
    dr::PCG32<UInt32> rng;
    MediumInteraction3f mei;
    Float target_ot;
    Mask use_rrt;
    Mask has_spectral_extinction;
    UInt32 sampled_medium_component;

    // Note that ``throughput`` is shared between algorithm and should be accumulated
    // accordingly. If used for volpathmis, new members and data types will need to
    // be introduced.
    UnpolarizedSpectrum throughput;

    DRJIT_STRUCT(TrackingState, ray, rng, mei, target_ot, use_rrt, \
        has_spectral_extinction, sampled_medium_component, throughput)
};

/**
 * \brief Signature of the tracking function callback accepted by
 * ``Medium::dda_track``.
 *
 * One call is one collision attempt within \c segment, not one segment.
 *
 * \param segment
 *      An extremum segment along a ray.
 * \param state
 *      The tracking state that holds interaction information and accumulates
 *      throughput and pdfs.
 * \param channel
 *      The channel to use for sampling.
 * \param active
 *      Represents the active lanes.
 *
 *
 * \return A pair (advance, active):
 *      advance:    If true, tracking has exited the segment and requires a
 *                  new one. If false, repeat the loop with the same segment.
 *      active:     Represent active lanes. Lanes that have sampled a real
 *                  interaction or terminated for other reasons will return
 *                  ``false``, prompting the termination of the traversal.
 */
template <typename Float, typename Spectrum,
          typename TrackState = TrackingState<Float, Spectrum>>
using TrackingFunction =
    std::pair<dr::mask_t<Float>, dr::mask_t<Float>>(
    const ExtremumSegment<Float, Spectrum>& /*segment*/,
    TrackState& /*state*/,
    const dr::uint32_array_t<Float>& /*channel*/,
    dr::mask_t<Float> /*active*/
);

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

/**
 * \brief Delta tracking over one extremum segment.
 *
 * A \ref TrackingFunction: pass it to ``Medium::dda_track``. Terminates a lane on a real
 * scattering event and records the sampled medium component.
 */
template <typename Float, typename Spectrum>
std::pair<dr::mask_t<Float>, dr::mask_t<Float>>
delta_track_segment(const ExtremumSegment<Float, Spectrum> &segment,
                    TrackingState<Float, Spectrum> &state,
                    const dr::uint32_array_t<Float> &channel,
                    dr::mask_t<Float> active) {
    using Mask                = dr::mask_t<Float>;
    using UnpolarizedSpectrum = unpolarized_spectrum_t<Spectrum>;

    UnpolarizedSpectrum &throughput = state.throughput;
    auto &rng                       = state.rng;
    auto &mei                       = state.mei;

    auto medium           = mei.medium;
    Mask act_spectral     = state.has_spectral_extinction && active;
    Mask act_not_spectral = !state.has_spectral_extinction && active;

    Float mint = dr::select(mei.is_valid(),
                            dr::maximum(segment.mint, mei.t), segment.mint);

    Float segment_ot = (segment.maxt - mint) * segment.majorant();
    Mask sampled     = (state.target_ot < segment_ot) && active;
    Float maxt       = segment.maxt;

    if (dr::any_or<true>(sampled))
        dr::masked(maxt, sampled) =
            mint + state.target_ot /
                       dr::maximum(segment.majorant(), dr::Epsilon<Float>);

    Float dt = maxt - mint;

    if (dr::any_or<true>(act_spectral)) {
        UnpolarizedSpectrum tr = dr::exp(-dt * segment.majorant());
        Float pdf              = index_spectrum<Float, Spectrum>(
            dr::select(sampled, tr * segment.majorant(), tr), channel);
        dr::masked(throughput, act_spectral) *= tr / pdf;
    }

    if (dr::any_or<true>(sampled)) {
        mei.t = maxt;
        mei.p = state.ray(maxt);

        auto medium_sample = medium->sample_scattering_properties(
            mei, segment.majorant(), rng.template next_float<Float>(sampled),
            sampled);

        UnpolarizedSpectrum &sigma_s = medium_sample.sigma_s;
        UnpolarizedSpectrum &sigma_n = medium_sample.sigma_n;
        UnpolarizedSpectrum &sigma_t = medium_sample.sigma_t;

        Float null_scatter_prob = dr::mean(sigma_n / segment.majorant());
        Mask null_scatter =
            (rng.template next_float<Float>(sampled) < null_scatter_prob) &&
            sampled;
        Mask real_scatter = !null_scatter && sampled;

        if (dr::any_or<true>(null_scatter && act_spectral))
            dr::masked(throughput, null_scatter && act_spectral) *=
                sigma_n / null_scatter_prob;

        if (dr::any_or<true>(real_scatter)) {
            if (dr::any_or<true>(act_spectral))
                dr::masked(throughput, real_scatter && act_spectral) *=
                    sigma_s / (1.f - null_scatter_prob);

            if (dr::any_or<true>(act_not_spectral))
                dr::masked(throughput, real_scatter && act_not_spectral) *=
                    sigma_s / sigma_t;

            dr::masked(state.sampled_medium_component, real_scatter) =
                medium_sample.sampled_component;

            active &= !real_scatter;
        }

        dr::masked(state.target_ot, sampled) =
            -dr::log(1.f - rng.template next_float<Float>(sampled));
    }

    dr::masked(mei.t, !sampled) = dr::Infinity<Float>;
    dr::masked(state.target_ot, !sampled && active) -= segment_ot;

    return { /*advance=*/!sampled, active };
}

/**
 * \brief Ratio tracking over one extremum segment, residual ratio tracking
 * when <tt>state.use_rrt</tt> is set.
 *
 * A \ref TrackingFunction. Never terminates a lane: a shadow ray runs to the
 * end of its range.
 */
template <typename Float, typename Spectrum>
std::pair<dr::mask_t<Float>, dr::mask_t<Float>>
ratio_track_segment(const ExtremumSegment<Float, Spectrum> &segment,
                    TrackingState<Float, Spectrum> &state,
                    const dr::uint32_array_t<Float> &channel,
                    dr::mask_t<Float> active) {
    using Mask                = dr::mask_t<Float>;
    using UnpolarizedSpectrum = unpolarized_spectrum_t<Spectrum>;

    UnpolarizedSpectrum &throughput = state.throughput;
    auto &rng                       = state.rng;
    auto &mei                       = state.mei;
    Mask use_rrt                    = state.use_rrt;

    auto medium           = mei.medium;
    Mask act_spectral     = state.has_spectral_extinction && active;
    Mask act_not_spectral = !state.has_spectral_extinction && active;

    Float control           = dr::select(use_rrt, segment.minorant(), 0.f);
    Float residual_majorant = segment.majorant() - control;

    Float mint = dr::select(mei.is_valid(),
                            dr::maximum(segment.mint, mei.t), segment.mint);

    Float segment_ot = (segment.maxt - mint) * residual_majorant;
    Mask sampled     = (state.target_ot < segment_ot) && active;
    Float maxt       = segment.maxt;

    if (dr::any_or<true>(sampled))
        dr::masked(maxt, sampled) =
            mint + state.target_ot /
                       dr::maximum(residual_majorant, dr::Epsilon<Float>);

    Float dt = maxt - mint;

    if (dr::any_or<true>(use_rrt))
        dr::masked(throughput, active && use_rrt) *= dr::exp(-dt * control);

    if (dr::any_or<true>(act_spectral)) {
        UnpolarizedSpectrum tr = dr::exp(-dt * residual_majorant);
        Float pdf              = index_spectrum<Float, Spectrum>(
            dr::select(sampled, tr * residual_majorant, tr), channel);
        dr::masked(throughput, act_spectral) *= tr / pdf;
    }

    if (dr::any_or<true>(sampled)) {
        mei.t = maxt;
        mei.p = state.ray(maxt);

        UnpolarizedSpectrum sigma_t, sigma_n;
        std::tie(std::ignore, std::ignore, sigma_t) =
            medium->get_scattering_coefficients(mei, sampled);
        sigma_n = segment.majorant() - sigma_t;

        if (dr::any_or<true>(act_spectral))
            dr::masked(throughput, sampled && act_spectral) *= sigma_n;

        if (dr::any_or<true>(act_not_spectral))
            dr::masked(throughput, sampled && act_not_spectral) *= dr::maximum(
                1.f - (sigma_t - control) / residual_majorant, 0.f);

        dr::masked(state.target_ot, sampled) =
            -dr::log(1.f - rng.template next_float<Float>(sampled));
    }

    dr::masked(mei.t, !sampled) = dr::Infinity<Float>;
    dr::masked(state.target_ot, !sampled && active) -= segment_ot;

    return { /*advance=*/!sampled, active };
}

NAMESPACE_END(mitsuba)
