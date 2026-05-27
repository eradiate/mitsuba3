#include <mitsuba/core/ray.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/core/random.h>
#include <mitsuba/render/emitter.h>
#include <mitsuba/render/integrator.h>
#include <mitsuba/render/records.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/render/eradiate/tracking.h>
#include <mitsuba/render/eradiate/extremum.h>
#include <mitsuba/render/eradiate/extremum_segment.h>


NAMESPACE_BEGIN(mitsuba)

/**
 * \brief This list of flags is used to differentiate the type of paths that
 * are present in the stack.
 */
enum class PathTypeFlag : uint32_t {
    /// Default when standard path tracing with no variance reduction method.
    Standard    = 0x0000,

    /// NLE: natural path from which Clones are spawned.
    Mother      = 0x0001,

    /// NLE: path cloned from a mother path.
    Clone       = 0x0002,

    /// Path that underwent a splitting operation.
    Split       = 0x0004,

    /// Set of path types that constitute the NLE variance reduction method
    NLE         = Mother | Clone,
};
MI_DECLARE_ENUM_OPERATORS(PathTypeFlag)

/**!

.. _integrator-eovolpath:

EO Volumetric path tracer (:monosp:`eovolpath`)
-----------------------------------------------

.. pluginparameters::

 * - max_depth
   - |int|
   - Specifies the longest path depth in the generated output image (where -1 corresponds to
     :math:`\infty`). A value of 1 will only render directly visible light sources. 2 will lead
     to single-bounce (direct-only) illumination, and so on. (Default: -1)

 * - rr_depth
   - |int|
   - Specifies the minimum path depth, after which the implementation will start to use the
     *russian roulette* path termination criterion. (Default: 5)

 * - rr_factor
   - |float|
   - Specifies the maximum probability to keep a path when russian roulette is evaluated.
     (Default: 0.97)

 * - enable_ddis
   - |bool|
   - Activate the DDIS variance reduction method. The `ddis_threshold` used to determine the
     probability of sampling usin the emitter direction is set in the media plugins.
     (Default: false)

 * - hide_emitters
   - |bool|
   - Hide directly visible emitters. (Default: no, i.e. |false|)

 * - enable_pbs
   - |bool|
   - Enable prediction-based path splitting (PBS). At each volumetric event,
     the predicted contribution of the scattered direction is used to determine whether
     to split the path into multiple independent copies or to perform russian rouletter.
     Each split path has a proportionally reduced weight. (Default: |false|)

 * - min_split_threshold
   - |float|
   - Minimum prediction weight required to trigger a split. Only paths whose predicted weight
     exceeds this value are split. Must be greater than 1 for splitting to produce more than one
     copy. (Default: 3.0)

 * - max_split_count
   - |int|
   - Maximum number of path copies created at a single splitting event. The actual count is
     ``min(max_split_count, floor(w_spl))``, where ``w_spl`` is the split prediction weight.
     (Default: 100)

 * - crit_rr_threshold
   - |float|
   - Weight threshold below which Russian Roulette is applied to split paths. Split paths whose
     current prediction weight falls below this value are stochastically terminated to limit the
     cost of low-weight copies. (Default: 0.03)

 * - min_rr_threshold
   - |float|
   - Minimum survival probability for split-path Russian Roulette. The survival probability is
     ``max(w_spl, min_rr_threshold)``, which ensures the kill probability never exceeds
     ``1 - min_rr_threshold`` and that surviving paths are not reweighted above their pre-split
     weight. (Default: 0.02)

 * - enable_nle
   - |bool|
   - Enable the N-tuple Local Estimate (NLE) variance reduction method. Each primary ray is traced as
     a *mother* path. At regular intervals along the mother's trajectory, a *clone* path is forked
     and traced independently to perform additional next-event estimation. (Default: |false|)

 * - first_clone_depth
   - |int|
   - Scatter depth at which the mother path creates its first clone. Clone creation then recurs
     every ``nee_per_clone`` scatters thereafter. (Default: 1)

 * - max_clone_depth
   - |int|
   - Maximum number of scattering events a clone is allowed to trace before it is terminated.
     Controls the amount of next-event estimation work performed per clone. (Default: 10)

 * - nee_per_clone
   - |int|
   - Interval, in scattering events, between successive clone creation events along the mother
     path. A new clone is spawned every ``nee_per_clone`` scatters starting from
     ``first_clone_depth``. (Default: 9)

This plugin provides a volumetric path tracer catered towards Earth Observation simualtions that can
be used to compute approximate solutions of the radiative transfer equation. Its implementation
makes use of multiple importance sampling to combine BSDF and phase function sampling with direct
illumination sampling strategies. On surfaces, it behaves exactly like the standard path tracer.

This integrator has special support for index-matched transmission events (i.e. surface scattering
events that do not change the direction of light). As a consequence, participating media enclosed by
a stencil shape are rendered considerably more efficiently when this shape
has a :ref:`null <bsdf-null>` or :ref:`thin dielectric <bsdf-thindielectric>` BSDF assigned
to it (as compared to, say, a :ref:`dielectric <bsdf-dielectric>` or
:ref:`roughdielectric <bsdf-roughdielectric>` BSDF).

In addition, it implements the VROOM variance reduction methods (Buras and Mayer, 2011) which reduce
noise in the presence of strongly peaked phase functions.

.. note:: This integrator does not implement good sampling strategies to render
    participating media with a spectrally varying extinction coefficient. For these cases,
    it is better to use the more advanced :ref:`volumetric path tracer with
    spectral MIS <integrator-volpathmis>`, which will produce in a significantly less noisy
    rendered image. Note however that it does not implement the EO features added in this integrator.

.. warning:: This integrator does not support forward-mode differentiation.
*/
template <typename Float, typename Spectrum>
class EOVolumetricPathIntegrator : public MonteCarloIntegrator<Float, Spectrum> {

public:
    MI_IMPORT_BASE(MonteCarloIntegrator, m_max_depth, m_rr_depth, m_hide_emitters)
    MI_IMPORT_TYPES(Scene, Sampler, Emitter, EmitterPtr, BSDF, BSDFPtr, Medium,
                    MediumPtr, PhaseFunction, PhaseFunctionContext, ExtremumStructure)

    using TrackingStateType = TrackingState<Float, Spectrum>;

private:

struct PathState {
    Mask active;
    UInt32 depth;
    Ray3f ray;
    Spectrum throughput;
    MediumPtr medium;
    Float eta;
    Interaction3f last_scatter_event;
    Float last_scatter_direction_pdf;
    Mask specular_chain;
    Mask valid_ray;
    Float split_weight_rr;
    UInt32 path_flag;
    UInt32 local_depth;

    DRJIT_STRUCT(PathState, active, depth, ray, throughput, \
        medium, eta, last_scatter_event, last_scatter_direction_pdf, \
        specular_chain, valid_ray, split_weight_rr, \
        path_flag, local_depth)
    };

    static constexpr size_t LS_STACK_SIZE = 4;

    using PathStack = dr::Array<PathState, LS_STACK_SIZE>;
    using CountStack = dr::Array<UInt32, LS_STACK_SIZE>;

struct LoopState {
        PathState current;
        Spectrum result;
        PathStack stack;
        CountStack counts;
        Int32 stack_counter;  // -1 means empty (no saved state)
        Sampler* sampler;
        Mask active;

        DRJIT_STRUCT(LoopState, current, result, stack, counts, stack_counter, sampler, active)
    };

public:

    EOVolumetricPathIntegrator(const Properties &props) : Base(props) {
        m_rr_factor = props.get<ScalarFloat>("rr_factor", 0.97f);

        if (m_rr_factor < 0.f || m_rr_factor > 1.f)
            Throw("`rr_factor` is outside range [0.,1.].");

        m_enable_ddis = props.get<bool>("enable_ddis", false);

        // Split Properties
        m_enable_pbs           = props.get<bool>("enable_pbs", false);
        m_max_split_count       = props.get<ScalarUInt32>("max_split_count", 100);
        m_min_split_threshold   = props.get<ScalarFloat>("min_split_threshold", 3.f);
        m_crit_rr_threshold     = props.get<ScalarFloat>("crit_rr_threshold", 0.03f);
        m_min_rr_threshold      = props.get<ScalarFloat>("min_rr_threshold", 0.02f);

        if (m_enable_pbs && (m_crit_rr_threshold <  0.f || m_crit_rr_threshold >= 1.f))
            Throw("`crit_rr_threshold` must be between 0 and 1.");

        if (m_enable_pbs && (m_min_rr_threshold <  0.f || m_min_rr_threshold >= 1.f))
            Throw("`min_rr_threshold` must be between 0 and 1.");

        if (m_enable_pbs && (m_min_split_threshold < 1.f))
            Throw("`min_split_threshold` must be greater than 1.");

        // NLE Properties
        m_enable_nle          = props.get<bool>("enable_nle", false);
        m_first_clone_depth  = props.get<ScalarUInt32>("first_clone_depth", 1);
        m_max_clone_depth    = props.get<ScalarUInt32>("max_clone_depth", 10);
        m_nee_per_clone      = props.get<ScalarUInt32>("nee_per_clone", 9);

        if (m_enable_nle && (m_max_clone_depth <= 1 || m_nee_per_clone <= 1))
            Throw("`max_clone_depth` and `nee_per_clone` must be larger than one");

        if (m_enable_nle && m_max_clone_depth < m_nee_per_clone )
            Throw("`max_clone_depth` must be larger or equal to `nee_per_clone`.");
    }

    /// Create seed and offsets that can be used to generate a new PCG32 rng.
    MI_INLINE
    std::pair<UInt64, UInt64> new_seed_offset(
            Float sample1, Float sample2) const {
        UInt32 s0 = UInt32(sample1 * 4294967296.f);  // [0,1) -> [0, 2^32)
        UInt32 s1 = UInt32(sample2 * 4294967296.f);
        UInt64 seed, offset;
        seed   = sample_tea_64(s0, s1);
        offset = sample_tea_64(s1, s0);
        return {seed, offset};
    }

    // === Stack helpers ===

    void write(Int32 index, LoopState &ls, const PathState &ps,
               UInt32 count, Mask active = true) const {
        if constexpr (!dr::is_jit_v<Float>) {
            if (active) {
                ls.stack[index]  = ps;
                ls.counts[index] = count;
            }
        } else {
            for (size_t i = 0; i < LS_STACK_SIZE; ++i) {
                Mask cond = (Int32(i) == index) && active;
                dr::masked(ls.stack[i],  cond) = ps;
                dr::masked(ls.counts[i], cond) = count;
            }
        }
    }

    void update_current(Int32 index, LoopState &ls, Mask active = true) const {
        if constexpr (!dr::is_jit_v<Float>) {
            dr::masked(ls.current, active) = ls.stack[int(index)];
        } else {
            for (size_t i = 0; i < LS_STACK_SIZE; ++i)
                ls.current = dr::select((Int32(i) == index) && active, ls.stack[i], ls.current);
        }
    }

    // === Splitting interface ===

    /// Push ps onto the stack at stack_counter+1.
    /// Returns the mask of lanes that successfully pushed (stack was not full).
    Mask push(LoopState &ls, const PathState &ps, UInt32 count, Mask active) const {
        Mask can_push = active && (ls.stack_counter < Int32(LS_STACK_SIZE) - 1);
        Int32 new_sc = ls.stack_counter + 1;
        write(new_sc, ls, ps, count, can_push);
        dr::masked(ls.stack_counter, can_push) += 1;
        return can_push;
    }

    /// Decrement count at the current top and pop if exhausted. Cascade through
    /// the stack if the count is equal to zero.
    /// Returns true for lanes whose stack is now empty (stack_counter < 0).
    Mask pop(LoopState &ls, Mask active) const {
        if constexpr (!dr::is_jit_v<Float>) {
            if (active) {
                for (int i = int(ls.stack_counter); i >= 0; --i) {
                    ls.counts[i]--;
                    if (ls.counts[i] > 0) break;
                    ls.stack_counter--;
                }
            }
        } else {
            for (int i = int(LS_STACK_SIZE) - 1; i >= 0; --i) {
                Mask is_top = active && (Int32(i) == ls.stack_counter);
                dr::masked(ls.counts[i], is_top) -= 1;
                dr::masked(ls.stack_counter, is_top && (ls.counts[i] == 0)) -= 1;
            }
        }
        return active && (ls.stack_counter < Int32(0));
    }

    void terminate_ray(LoopState &ls) const {
        Mask terminated = ls.active && !ls.current.active;
        Mask empty = pop(ls, terminated);
        ls.active &= !empty;

        // Clip to 0 so update_current never receives a -1 index (UB in debug).
        dr::masked(ls.stack_counter, empty) = Int32(0);

        // In scalar mode stack_counter may still be -1 for fully-drained stacks;
        // the clipping above handles it, but the early return is a cheap guard.
        if (dr::none_or<false>(ls.active))
            return;

        update_current(ls.stack_counter, ls, terminated && ls.active);
    }

    std::pair<Spectrum, Mask> sample(const Scene *scene,
                                     Sampler *sampler,
                                     const RayDifferential3f &ray_,
                                     const Medium *initial_medium,
                                     Float * /* aovs */,
                                     Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::SamplingIntegratorSample, active);
        // If there is an environment emitter and emitters are visible: all rays will be valid
        // Otherwise, it will depend on whether a valid interaction is sampled
        Mask valid_ray = !m_hide_emitters && (scene->environment() != nullptr);

        // For now, don't use ray differentials
        Ray3f ray = ray_;

        // Tracks radiance scaling due to index of refraction changes
        Float eta(1.f);

        Spectrum throughput(1.f), result(0.f);
        MediumPtr medium = initial_medium;
        Mask specular_chain = active && !m_hide_emitters;
        UInt32 depth = 0, local_depth = 0;

        UInt32 channel = 0;
        if (is_rgb_v<Spectrum>) {
            uint32_t n_channels = (uint32_t) dr::size_v<Spectrum>;
            channel = (UInt32) dr::minimum(sampler->next_1d(active) * n_channels, n_channels - 1);
        }

        Interaction3f last_scatter_event = dr::zeros<Interaction3f>();
        Float last_scatter_direction_pdf = 1.f;

        // VROOM variables
        UInt32 path_flag =
            dr::select(m_enable_nle, +PathTypeFlag::Mother, +PathTypeFlag::Standard);
        Float split_weight_rr = 1.f;

        /* Set up a Dr.Jit loop (optimizes away to a normal loop in scalar mode,
           generates wavefront or megakernel renderer based on configuration).
           Register everything that changes as part of the loop here */
        PathState ps = {
            active,
            depth,
            ray,
            throughput,
            medium,
            eta,
            last_scatter_event,
            last_scatter_direction_pdf,
            specular_chain,
            valid_ray,
            split_weight_rr,
            path_flag,
            local_depth
        };

        PathStack stack;
        CountStack counts;
        Int32 stack_counter = Int32(-1);

        LoopState ls {
            ps,
            result,
            stack,
            counts,
            stack_counter,
            sampler,
            active
        };

        // Base sentinel: push the initial path state with count 1 so that
        // pop() cascades through it and deactivates the lane when done.
        push(ls, ps, UInt32(1), active);

        dr::tie(ls) = dr::while_loop(dr::make_tuple(ls),
            [](const LoopState& ls) { return ls.active; }, // LoopState active flag
            [this, scene, channel](LoopState& ls) {

            Mask& active = ls.current.active; // PathState active flag
            UInt32& depth = ls.current.depth;
            Ray3f& ray = ls.current.ray;
            Spectrum& throughput = ls.current.throughput;
            Spectrum& result = ls.result;
            MediumPtr& medium = ls.current.medium;
            Float& eta = ls.current.eta;
            Interaction3f& last_scatter_event = ls.current.last_scatter_event;
            Float& last_scatter_direction_pdf = ls.current.last_scatter_direction_pdf;
            Mask& specular_chain = ls.current.specular_chain;
            Mask& valid_ray = ls.current.valid_ray;
            Sampler* sampler = ls.sampler;
            Float& split_weight_rr = ls.current.split_weight_rr;
            UInt32& path_flag    = ls.current.path_flag;
            UInt32& local_depth  = ls.current.local_depth;

            MediumInteraction3f mei;
            SurfaceInteraction3f si;

            Mask is_mother = has_flag(path_flag, PathTypeFlag::Mother);
            Mask is_clone = has_flag(path_flag, PathTypeFlag::Clone);
            Mask is_standard = !is_mother && !is_clone;

            // ----------------- Handle termination of paths ------------------
            active &= dr::any(unpolarized_spectrum(throughput) != 0.f);

            Float q;
            Mask perform_rr;

            // Standard Russian roulette:
            // Try to keep path weights equal to one, while accounting for the
            // solid angle compression at refractive index boundaries. Stop with at least some
            // probability to avoid getting stuck (e.g. due to total internal reflection)
            dr::masked(q, is_standard) =
                dr::minimum(dr::max(unpolarized_spectrum(throughput)) * dr::square(eta), m_rr_factor);
            dr::masked(perform_rr, is_standard) =
                (depth > (uint32_t) m_rr_depth);

            // PBS Russian Roulette:
            // Eliminate path that are predicted to have a small contribution.
            Mask do_pbs_rr = m_enable_pbs && (is_clone || has_flag(path_flag, PathTypeFlag::Split));
            dr::masked(q, do_pbs_rr) =
                dr::maximum(split_weight_rr, m_min_rr_threshold);
            dr::masked(perform_rr, do_pbs_rr) =
                (split_weight_rr < m_crit_rr_threshold);

            // apply russian roulette
            active &= sampler->next_1d(active) < q || !perform_rr;
            dr::masked(throughput, perform_rr) *= dr::rcp(dr::detach(q));

            active &= depth < (uint32_t) m_max_depth;
            if (dr::none_or<false>(active)) {
                terminate_ray(ls);
                return;
            }

            // ----------------------- Sampling the RTE -----------------------
            Mask active_medium  = active && (medium != nullptr);
            Mask active_surface = active && !active_medium;

            Mask act_null_scatter = false, act_medium_scatter = false,
                 escaped_medium = false;

            dr::masked(si, active) = scene->ray_intersect(ray, active);
            dr::masked(ray.maxt, active) = si.t;

            // If the medium does not have a spectrally varying extinction,
            // we can perform a few optimizations to speed up rendering
            Mask is_spectral = active_medium;
            Mask not_spectral = false;
            Float ddis_threshold = -1.f;
            split_weight_rr = 1.f;

            if (dr::any_or<true>(active_medium)) {
                is_spectral &= medium->has_spectral_extinction();
                not_spectral = !is_spectral && active_medium;

                dr::masked(ddis_threshold, active_medium && m_enable_ddis) =
                    medium->ddis_threshold();
            }

            Float mint, maxt;
            if (dr::any_or<true>(active_medium)) {
                std::tie(mei, mint, maxt) =
                    medium->prepare_medium_traversal(ray, active_medium);
                escaped_medium = active_medium && !dr::isfinite(maxt);
                active_medium &= !escaped_medium;
            }

            if (dr::any_or<true>(active_medium)) {
                // Prepare Extremum traversal
                auto extremum = medium->extremum_structure();

                Float sample1 = sampler->next_1d();
                Float sample2 = sampler->next_1d();
                auto [seed, offset] = new_seed_offset(sample1, sample2);
                PCG32<UInt32> rng;
                rng.seed(seed, offset);

                Float target_ot = -dr::log( 1.f - sampler->next_1d(active_medium));

                TrackingStateType state {
                    ray,
                    rng,
                    mei,
                    target_ot,
                    medium->use_rrt(),
                    medium->has_spectral_extinction(),
                    /*throughput=*/UnpolarizedSpectrum(1.f),
                };

                // Traverse extremum segments and perform delta tracking
                state = extremum->traverse_extremum(
                ray, mint, maxt, channel, state,
                [](const ExtremumSegment &segment, TrackingStateType &state,
                    const UInt32 &channel, Mask active) {
                    UnpolarizedSpectrum &throughput = state.throughput;
                    PCG32<UInt32> &rng              = state.rng;
                    MediumInteraction3f &mei        = state.mei;

                    MediumPtr medium = mei.medium;
                    Mask act_spectral = state.has_spectral_extinction && active;
                    Mask act_not_spectral = !state.has_spectral_extinction && active;

                    // Check if the last iteration had a valid interaction
                    // within the segment
                    Float mint = dr::select(
                        mei.is_valid(), dr::maximum(segment.mint, mei.t),
                        segment.mint);

                    Float segment_ot =
                        (segment.maxt - mint) * segment.majorant();
                    Mask sampled = (state.target_ot < segment_ot) && active;
                    Float maxt   = segment.maxt;

                    if (dr::any_or<true>(sampled)) {
                        dr::masked(maxt, sampled) =
                            mint + state.target_ot /
                                        dr::maximum(segment.majorant(),
                                                    dr::Epsilon<Float>);
                    }

                    Float dt = maxt - mint;

                    if (dr::any_or<true>(act_spectral)) {
                        // Accumulate transmittance in the throughput and
                        // pdf (spectral only).
                        UnpolarizedSpectrum tr =
                            dr::exp(-dt * segment.majorant());
                        Float pdf = index_spectrum<Float, Spectrum>(
                            dr::select(sampled, tr * segment.majorant(), tr),
                            channel);
                        dr::masked(throughput, act_spectral) *= tr / pdf;
                    }

                    if (dr::any_or<true>(sampled)) {

                        mei.t = maxt;
                        mei.p = state.ray(maxt);

                        // Retrieve scattering coefficients at position.
                        UnpolarizedSpectrum sigma_s, sigma_n, sigma_t;
                        std::tie(sigma_s, std::ignore, sigma_t) =
                            medium->get_scattering_coefficients(mei, sampled);
                        sigma_n = segment.majorant() - sigma_t;

                        // Sample event type
                        Float null_scatter_prob =
                            dr::mean(sigma_n / segment.majorant());
                        Mask null_scatter =
                            (rng.template next_float<Float>(sampled) < null_scatter_prob)
                            && sampled;
                        Mask real_scatter = !null_scatter && sampled;

                        // Accumulate throughput and pdf given the event
                        // type and is_spectral.
                        if (dr::any_or<true>(null_scatter && act_spectral)) {
                            dr::masked(throughput, null_scatter && act_spectral) *=
                                sigma_n / null_scatter_prob;
                        }

                        if (dr::any_or<true>(real_scatter)) {
                            if (dr::any_or<true>(act_spectral)) {
                                dr::masked(throughput, real_scatter && act_spectral) *=
                                    sigma_s / (1.0f - null_scatter_prob);
                            }

                            if (dr::any_or<true>(act_not_spectral)) {
                                // pdf = sigma_t / sigma_maj, sigma_maj gets
                                // cancelled from sigma_s/sigma_maj
                                dr::masked(throughput, real_scatter && act_not_spectral) *=
                                    sigma_s / sigma_t;
                            }

                            // disable the loop once we encounter a real
                            // scattering interaction
                            active &= !real_scatter;
                        }

                        dr::masked(state.target_ot, sampled) =
                            -dr::log(1.f - state.rng.template next_float<Float>(sampled));
                    }

                    dr::masked(mei.t, !sampled) = dr::Infinity<Float>;
                    dr::masked(state.target_ot, !sampled && active) -=
                        segment_ot;

                    Mask step = !sampled;
                    return std::pair<Mask, Mask>(step, active);
                }, active_medium);

                // Update throughput by the transmittance and pdf weight
                dr::masked(throughput, active_medium) *= state.throughput;
                dr::masked(mei, active_medium) = state.mei;

                escaped_medium |= active_medium && !mei.is_valid();
                active_medium &= mei.is_valid();

                act_medium_scatter = !escaped_medium && active_medium;
                dr::masked(depth, act_medium_scatter) += 1;
                dr::masked(local_depth, act_medium_scatter) += 1;
                dr::masked(last_scatter_event, act_medium_scatter) = mei;
            }

            // Dont estimate lighting if we exceeded number of bounces
            active &= depth < (uint32_t) m_max_depth;
            active &= !(has_flag(path_flag, PathTypeFlag::Clone)
                        && local_depth > (uint32_t)m_max_clone_depth);
            act_medium_scatter &= active;
            if (dr::any_or<true>(act_medium_scatter)) {

                PhaseFunctionContext phase_ctx(sampler);
                auto phase = mei.medium->phase_function();
                auto ddis_phase_function = mei.medium->ddis_phase_function();

                // --------------------- NLE setup ---------------------
                Mask enable_nle = act_medium_scatter && m_enable_nle;

                // Clone creation: mother at/past first_clone_depth, spaced by nee_per_clone
                Mask create_clone = active && is_mother
                    && depth >= (uint32_t)m_first_clone_depth
                    && ((depth - (uint32_t)m_first_clone_depth) % (uint32_t)m_nee_per_clone == 0);

                // --------------------- Emitter sampling ---------------------
                Mask sample_emitters = mei.medium->use_emitter_sampling();
                valid_ray |= act_medium_scatter;
                specular_chain &= !act_medium_scatter;
                specular_chain |= act_medium_scatter && !sample_emitters;

                Mask active_e = act_medium_scatter && sample_emitters;
                Mask perform_ddis = ddis_threshold > 0.f && active_e;
                Mask active_nee = active_e;

                if (dr::any_or<true>(enable_nle)) {
                    // DDIS: restrict to create_clone or existing clone when NLE is on
                    dr::masked(perform_ddis, enable_nle) &=
                        create_clone || is_clone;

                    // NEE gating:
                    //  Mother: before first_clone_depth only
                    //  First clone (creation_depth == first_clone_depth): always NEE
                    //  Subsequent clones: last nee_per_clone scatters only
                    // local_depth is already incremented → 1-indexed at NEE time
                    dr::masked(active_nee, enable_nle) &=
                        (depth <= (uint32_t)m_first_clone_depth)
                        || (is_clone
                        && ( (depth - local_depth) == (uint32_t) m_first_clone_depth
                           || local_depth > (uint32_t) (m_max_clone_depth - m_nee_per_clone)));
                }

                MediumInteraction3f ddis_mei = mei;

                if (dr::any_or<true>(active_e)) {
                    // Externalize direction sampling so ds.d is available for DDIS setup
                    auto [ds, emitter_val] = scene->sample_emitter_direction(
                        mei, sampler->next_2d(active_e), false, active_e);
                    dr::masked(emitter_val, ds.pdf == 0.f) = 0.f;

                    if (dr::any_or<true>(perform_ddis)) {
                        // DDIS incoming direction towards emitter
                        dr::masked(ddis_mei.wi,       perform_ddis) = -ds.d;
                        dr::masked(ddis_mei.sh_frame, perform_ddis) = Frame3f(-ds.d);
                    }

                    if (dr::any_or<true>(active_nee)) {

                        auto [phase_val, phase_pdf] = phase->eval_pdf(phase_ctx, mei, ds.d, active_nee);

                        auto [emitted, ds_out] = sample_emitter(
                            ds, emitter_val, mei, scene, sampler, medium, channel, active_nee);
                        dr::masked(result, active_nee) += throughput * phase_val * emitted *
                            mis_weight(ds.pdf, dr::select(ds.delta, 0.f, phase_pdf));
                    }
                }

                // ------------------ Phase function sampling -----------------
                dr::masked(phase, !act_medium_scatter) = nullptr;

                auto [wo, phase_weight, phase_pdf] = phase->sample(phase_ctx, mei,
                        sampler->next_1d(act_medium_scatter),
                        sampler->next_2d(act_medium_scatter),
                        act_medium_scatter);

                // --------------------- NLE clone creation --------------------
                if (dr::any_or<true>(create_clone)) {
                    PathState mps = ls.current;
                    mps.ray                         = mei.spawn_ray(wo);
                    mps.throughput                  *= phase_weight;
                    mps.last_scatter_direction_pdf  = phase_pdf;
                    // Requires a count of two to account for the continuing clone.
                    write(Int32(0), ls, mps, 2, create_clone);

                    dr::masked(path_flag,   create_clone) = +PathTypeFlag::Clone;
                    dr::masked(local_depth, create_clone) = UInt32(0);

                    is_mother = !create_clone;
                    is_clone = create_clone;
                }

                // ------------------------ DDIS -------------------------------
                if (dr::any_or<true>(perform_ddis)) {
                    Spectrum natural_val = phase_weight * phase_pdf;
                    Float natural_pdf = phase_pdf;

                    Float eps = sampler->next_1d(active_medium);
                    Mask active_ddis = (eps < ddis_threshold) && act_medium_scatter && perform_ddis;

                    if(dr::any_or<true>(active_ddis)) {
                        auto [ddis_wo, ddis_weight, ddis_pdf] =
                            ddis_phase_function->sample(phase_ctx, ddis_mei,
                            sampler->next_1d(active_ddis),
                            sampler->next_2d(active_ddis),
                            active_ddis);

                        dr::masked(wo, active_ddis) = ddis_wo;

                        auto [new_natural_val, new_natural_pdf] =
                            phase->eval_pdf(phase_ctx, mei, wo, active_ddis);
                        dr::masked(natural_val, active_ddis) = new_natural_val;
                        dr::masked(natural_pdf, active_ddis) = new_natural_pdf;
                    }

                    auto [ddis_val, ddis_pdf] =
                        ddis_phase_function->eval_pdf(phase_ctx, ddis_mei, wo, perform_ddis);

                    dr::masked(phase_pdf, perform_ddis) =
                        ((1.f - ddis_threshold) * natural_pdf + ddis_threshold * ddis_pdf);
                    dr::masked(phase_weight, perform_ddis) = natural_val / phase_pdf;

                    // split_weight computed for all perform_ddis lanes; enable_pbs gates PBS below
                    dr::masked(split_weight_rr, perform_ddis) =
                        dr::max(unpolarized_spectrum(ddis_val * throughput));
                }

                act_medium_scatter &= phase_pdf > 0.f;

                // ----------------------- Update Path -------------------------

                Ray3f new_ray  = mei.spawn_ray(wo);
                dr::masked(ray, act_medium_scatter) = new_ray;
                dr::masked(last_scatter_direction_pdf, act_medium_scatter) = phase_pdf;
                dr::masked(throughput, act_medium_scatter) *= phase_weight;
            }

            // --------------------- Surface Interactions ---------------------
            // Interact with the surface only if we haven't interacted with the
            // medium before and there is a valid intersection (already accounted).
            active_surface |= escaped_medium;

            if (dr::any_or<true>(active_surface)) {
                // ---------------------- Hide area emitters ----------------------
                if (m_hide_emitters && dr::any_or<true>(depth == 0u)) {
                    // Are we on the first segment and did we hit an area emitter?
                    // If so, skip all area emitters along this ray
                    Mask skip_emitters = si.is_valid() &&
                                         (si.shape->emitter() != nullptr) &&
                                         (depth == 0);

                    if (dr::any_or<true>(skip_emitters)) {
                        Ray3f ray = si.spawn_ray(ray.d);
                        PreliminaryIntersection3f pi =
                            Base::skip_area_emitters(scene, ray, true, skip_emitters);
                        SurfaceInteraction3f si_after_skip =
                            pi.compute_surface_interaction(ray, +RayFlags::All, skip_emitters);
                        dr::masked(si, skip_emitters) = si_after_skip;
                    }
                }

                // ---------------- Intersection with emitters ----------------
                Mask ray_from_camera = active_surface && (depth == 0u);
                Mask count_direct = ray_from_camera || specular_chain;
                EmitterPtr emitter = si.emitter(scene);
                Mask active_e = active_surface && (emitter != nullptr) &&
                                !((depth == 0u) && m_hide_emitters);

                if (dr::any_or<true>(active_e)) {
                    Float emitter_pdf = 1.0f;
                    if (dr::any_or<true>(active_e && !count_direct)) {
                        // Get the PDF of sampling this emitter using next event estimation
                        DirectionSample3f ds(scene, si, last_scatter_event);
                        emitter_pdf = scene->pdf_emitter_direction(last_scatter_event, ds, active_e);
                    }
                    Spectrum emitted = emitter->eval(si, active_e);
                    Spectrum contrib = dr::select(count_direct, throughput * emitted,
                                                  throughput * mis_weight(last_scatter_direction_pdf, emitter_pdf) * emitted);
                    dr::masked(result, active_e) += contrib;
                }
            }
            active_surface &= si.is_valid();
            active &= !(has_flag(path_flag, PathTypeFlag::Clone)
                        && (local_depth + 1)> (uint32_t)m_max_clone_depth);
            active_surface &= active;

            if (dr::any_or<true>(active_surface)) {

                // --------------------- NLE setup (surface) ---------------------
                Mask enable_nle = active_surface && m_enable_nle;


                // Tentative clone creation (non_null_bsdf not yet known)
                Mask create_clone = enable_nle
                    && is_mother
                    && (depth + 1) >= (uint32_t)m_first_clone_depth
                    && ((depth + 1 - (uint32_t)m_first_clone_depth) % (uint32_t)m_nee_per_clone == 0);

                // --------------------- Emitter sampling ---------------------
                BSDFContext ctx;
                BSDFPtr bsdf  = si.bsdf(ray);

                Mask active_e = active_surface && has_flag(bsdf->flags(), BSDFFlags::Smooth) && (depth + 1 < (uint32_t) m_max_depth);
                Mask perform_ddis = (medium != nullptr) && (ddis_threshold > 0.f) && active_e;
                Mask active_nee = active_e;

                if (dr::any_or<true>(enable_nle)) {
                    dr::masked(perform_ddis, enable_nle) &=
                        create_clone || is_clone;

                    dr::masked(active_nee, enable_nle) &=
                        ((depth + 1) <= (uint32_t)m_first_clone_depth)
                        || (is_clone
                            && ( (depth - local_depth) == (uint32_t) m_first_clone_depth
                            || (local_depth + 1) > (uint32_t) (m_max_clone_depth - m_nee_per_clone)));
                }

                MediumInteraction3f ddis_mei = dr::zeros<MediumInteraction3f>();
                ddis_mei.p = si.p;
                ddis_mei.time = si.time;
                ddis_mei.wavelengths = si.wavelengths;
                ddis_mei.medium = medium;
                PhaseFunctionContext phase_ctx(sampler);

                if (likely(dr::any_or<true>(active_e))) {
                    auto [ds, emitter_val] = scene->sample_emitter_direction(
                        si, sampler->next_2d(active_e), false, active_e);
                    dr::masked(emitter_val, ds.pdf == 0.f) = 0.f;

                    if (dr::any_or<true>(perform_ddis)) {
                        // Set ddis incoming direction towards emitter
                        dr::masked(ddis_mei.wi, perform_ddis) = -ds.d;
                        dr::masked(ddis_mei.sh_frame, perform_ddis) = Frame3f(-ds.d);
                    }

                    // Query the BSDF for that emitter-sampled direction
                    Vector3f wo       = si.to_local(ds.d);
                    Spectrum bsdf_val = bsdf->eval(ctx, si, wo, active_e);
                    bsdf_val = si.to_world_mueller(bsdf_val, -wo, si.wi);

                    // Determine probability of having sampled that same
                    // direction using BSDF sampling.
                    Float bsdf_pdf = bsdf->pdf(ctx, si, wo, active_e);

                    if (dr::any_or<true>(active_nee)) {
                        auto [emitted, ds_out] = sample_emitter(
                            ds, emitter_val, si, scene, sampler, medium, channel, active_nee);
                        dr::masked(result, active_nee) += throughput * bsdf_val *
                            mis_weight(ds.pdf, dr::select(ds.delta, 0.f, bsdf_pdf)) * emitted;
                    }
                }

                // ----------------------- BSDF sampling ----------------------
                auto [bs, bsdf_val] = bsdf->sample(ctx, si, sampler->next_1d(active_surface),
                                                   sampler->next_2d(active_surface), active_surface);

                Mask non_null_bsdf = active_surface && !has_flag(bs.sampled_type, BSDFFlags::Null);

                // Finalize clone creation condition
                create_clone &= non_null_bsdf;

                // --------------------- NLE clone creation --------------------
                if (dr::any_or<true>(create_clone)) {
                    Spectrum bsdf_val_world =
                        si.to_world_mueller(bsdf_val, -bs.wo, si.wi);

                    PathState mps = ls.current;
                    mps.ray                         = si.spawn_ray(si.to_world(bs.wo));
                    mps.throughput                  *= bsdf_val_world;
                    mps.eta                         *= bs.eta;
                    mps.last_scatter_direction_pdf  = bs.pdf;

                    // create_clone already ensures that this is a non Null BSDF
                    mps.depth               += 1;
                    mps.last_scatter_event  = si;
                    mps.valid_ray           = true;

                    mps.specular_chain = specular_chain || has_flag(bs.sampled_type, BSDFFlags::Delta);
                    mps.specular_chain &= !(active_surface && has_flag(bs.sampled_type, BSDFFlags::Smooth));

                    Mask has_medium_trans = create_clone && active_surface && si.is_medium_transition();
                    dr::masked(mps.medium, has_medium_trans) = si.target_medium(si.spawn_ray(si.to_world(bs.wo)).d);

                    write(Int32(0), ls, mps, UInt32(2), create_clone);

                    dr::masked(path_flag,      create_clone) = (uint32_t)PathTypeFlag::Clone;
                    dr::masked(local_depth,    create_clone) = UInt32(0);

                    is_mother = !create_clone;
                    is_clone = create_clone;
                }

                // ------------------------- DDIS ------------------------------
                if (dr::any_or<true>(perform_ddis)) {

                    auto ddis_phase_function = medium->ddis_phase_function();
                    Spectrum natural_val = bsdf_val * bs.pdf;
                    Float natural_pdf = bs.pdf;

                    // DDIS is active if the emitter is in the same lobe as the ray,
                    // and that the BSDF is neither Null or Delta.
                    Float eps = sampler->next_1d(perform_ddis);
                    Mask active_ddis = (eps < ddis_threshold) && perform_ddis;
                    active_ddis &= dr::dot(-ray.d, -ddis_mei.wi) > 0.f;
                    active_ddis &= !has_flag(bs.sampled_type, BSDFFlags::Null)
                                && !has_flag(bs.sampled_type, BSDFFlags::Delta);

                    // Sample from the ddis phase function
                    if (dr::any_or<true>(active_ddis)) {
                        auto [wo, phase_weight, phase_pdf] =
                            ddis_phase_function->sample(phase_ctx, ddis_mei,
                                sampler->next_1d(active_ddis),
                                sampler->next_2d(active_ddis),
                                active_ddis);

                        dr::masked(bs.wo, active_ddis)  = si.to_local(wo);
                        dr::masked(bs.eta, active_ddis) = 1.f;

                        BSDFContext ddis_ctx;
                        ddis_ctx.type_mask = +BSDFFlags::Reflection;

                        // Update the bsdf value with the evalutation from the new outgoing direction
                        auto [new_natural_val, new_natural_pdf] =
                            bsdf->eval_pdf(ddis_ctx, si , bs.wo, active_ddis);
                        dr::masked(natural_val, active_ddis) = new_natural_val;
                        dr::masked(natural_pdf, active_ddis) = new_natural_pdf;
                    }

                    auto [ddis_val, ddis_pdf] =
                        ddis_phase_function->eval_pdf(phase_ctx, ddis_mei,
                                                      si.to_world(bs.wo),
                                                      perform_ddis);

                    dr::masked(bs.pdf, perform_ddis) =
                        ((1.f - ddis_threshold) * natural_pdf + ddis_threshold * ddis_pdf);
                    dr::masked(bsdf_val, perform_ddis) = natural_val / bs.pdf;

                    dr::masked(split_weight_rr, perform_ddis) =
                        dr::max(unpolarized_spectrum(ddis_val * throughput));
                }

                bsdf_val = si.to_world_mueller(bsdf_val, -bs.wo, si.wi);

                dr::masked(throughput, active_surface) *= bsdf_val;
                dr::masked(eta, active_surface) *= bs.eta;

                Ray3f bsdf_ray                  = si.spawn_ray(si.to_world(bs.wo));
                dr::masked(ray, active_surface) = bsdf_ray;

                dr::masked(depth, non_null_bsdf) += 1;
                dr::masked(local_depth, non_null_bsdf) += 1;

                // update the last scatter PDF event if we encountered a non-null scatter event
                dr::masked(last_scatter_event, non_null_bsdf) = si;
                dr::masked(last_scatter_direction_pdf, non_null_bsdf) = bs.pdf;

                valid_ray |= non_null_bsdf;
                specular_chain |= non_null_bsdf && has_flag(bs.sampled_type, BSDFFlags::Delta);
                specular_chain &= !(active_surface && has_flag(bs.sampled_type, BSDFFlags::Smooth));
                act_null_scatter |= active_surface && has_flag(bs.sampled_type, BSDFFlags::Null);
                Mask has_medium_trans                = active_surface && si.is_medium_transition();
                dr::masked(medium, has_medium_trans) = si.target_medium(ray.d);
            }

            // --------------- Prediction Based Splitting ------------------

            Mask enable_pbs = (active_surface || active_medium) && m_enable_pbs && !is_mother;
            Mask split = enable_pbs && split_weight_rr > m_min_split_threshold;
            split &= ls.stack_counter < Int32(LS_STACK_SIZE) - 1;

            if (dr::any_or<true>(split)) {
                UInt32 split_count = dr::minimum(m_max_split_count, UInt32(split_weight_rr));
                dr::masked(throughput, split) /= Float(split_count);
                dr::masked(path_flag, split) |= +PathTypeFlag::Split;
                push(ls, ls.current, split_count, split);
            }

            active &= (active_surface | active_medium);
            if (dr::any_or<true>(!active)) {
                terminate_ray(ls);
            }
        },
        "Volpath integrator");

        return { ls.result, ls.current.valid_ray };
    }

    /// Inner overload: ds and emitter_val pre-sampled; does shadow/transmittance traversal only.
    template <typename Interaction>
    std::tuple<Spectrum, DirectionSample3f>
    sample_emitter(DirectionSample3f ds, Spectrum emitter_val,
                   const Interaction &ref_interaction, const Scene *scene,
                   Sampler *sampler, MediumPtr medium,
                   UInt32 channel, Mask active) const {
        Spectrum transmittance(1.0f);

        dr::masked(emitter_val, ds.pdf == 0.f) = 0.f;
        active &= (ds.pdf != 0.f);

        if (dr::none_or<false>(active)) {
            return { emitter_val, ds };
        }

        Ray3f ray = ref_interaction.spawn_ray_to(ds.p);
        Float max_dist = ray.maxt;

        // Potentially escaping the medium if this is the current medium's boundary
        if constexpr (std::is_convertible_v<Interaction, SurfaceInteraction3f>)
            dr::masked(medium, ref_interaction.is_medium_transition()) =
                ref_interaction.target_medium(ray.d);

        Float total_dist = 0.f;
        SurfaceInteraction3f si = dr::zeros<SurfaceInteraction3f>();
        DirectionSample3f dir_sample = ds;

        struct LoopState {
            Mask active;
            Ray3f ray;
            Float total_dist;
            MediumPtr medium;
            SurfaceInteraction3f si;
            Spectrum transmittance;
            Sampler* sampler;

            DRJIT_STRUCT(LoopState, active, ray, total_dist, \
                medium, si, transmittance, sampler)
        } ls = {
            active,
            ray,
            total_dist,
            medium,
            si,
            transmittance,
            sampler
        };

        dr::tie(ls) = dr::while_loop(dr::make_tuple(ls),
            [](const LoopState& ls) { return dr::detach(ls.active); },
            [this, scene, channel, max_dist](LoopState& ls) {

            Mask& active = ls.active;
            Ray3f& ray = ls.ray;
            Float& total_dist = ls.total_dist;
            MediumPtr& medium = ls.medium;
            SurfaceInteraction3f& si = ls.si;
            Spectrum& transmittance = ls.transmittance;
            Sampler* sampler = ls.sampler;

            Float remaining_dist = max_dist - total_dist;
            ray.maxt = remaining_dist;
            active &= remaining_dist > 0.f;
            if (dr::none_or<false>(active))
                return;

            Mask active_medium  = active && (medium != nullptr);

            dr::masked(si, active) = scene->ray_intersect(ray, active);
            dr::masked(ray.maxt, active) = dr::minimum(si.t, remaining_dist);
            dr::masked(total_dist, active) += ray.maxt;

            MediumInteraction3f mei;
            Float mint, maxt;
            if (dr::any_or<true>(active_medium)) {
                std::tie(mei, mint, maxt) =
                    medium->prepare_medium_traversal(ray, active_medium);
                active_medium &= dr::isfinite(maxt);
            }

            if (dr::any_or<true>(active_medium)) {
                // Prepare extremum traversal
                auto extremum = medium->extremum_structure();

                Float sample1 = sampler->next_1d();
                Float sample2 = sampler->next_1d();
                auto [seed, offset] = new_seed_offset(sample1, sample2);
                PCG32<UInt32> rng;
                rng.seed(seed, offset);

                Float target_ot = -dr::log( 1.f - sampler->next_1d(active_medium));

                TrackingStateType state {
                    ray,
                    rng,
                    mei,
                    target_ot,
                    medium->use_rrt(),
                    medium->has_spectral_extinction(),
                    /*throughput=*/UnpolarizedSpectrum(1.f),
                };

                // Unified Ratio Tracking and Residual Ratio Tracking approach:
                // if `use_rrt` is false, set control to 0., which automatically
                // devolves the algorithm to Ratio Tracking.
                state = extremum->traverse_extremum(
                ray, mint, maxt, channel, state,
                [](const ExtremumSegment& segment, TrackingStateType& state,
                   const UInt32& channel, Mask active) {

                    UnpolarizedSpectrum &throughput = state.throughput;
                    PCG32<UInt32> &rng       = state.rng;
                    MediumInteraction3f& mei = state.mei;
                    Mask use_rrt             = state.use_rrt;
                    MediumPtr medium         = mei.medium;
                    Mask act_spectral =
                        state.has_spectral_extinction && active;
                    Mask act_not_spectral =
                        !state.has_spectral_extinction && active;

                    Float control = dr::select(use_rrt, segment.minorant(), 0.f);
                    Float residual_majorant = segment.majorant() - control;

                    Float mint = dr::select(
                        mei.is_valid(),
                        dr::maximum(segment.mint, mei.t),
                        segment.mint
                    );

                    Float segment_ot = (segment.maxt - mint)*residual_majorant;
                    Mask sampled = (state.target_ot < segment_ot) && active;
                    Float maxt = segment.maxt;

                    if (dr::any_or<true>(sampled)) {
                        dr::masked(maxt, sampled) =
                            mint +
                            state.target_ot / dr::maximum(residual_majorant,
                                                          dr::Epsilon<Float>);
                    }

                    Float dt = maxt - mint;

                    if(dr::any_or<true>(use_rrt)) {
                        // We accumulate control transmittance every interaction
                        // instead of accumulating the optical thickness.
                        // Choice of simplicity vs performance.
                        dr::masked(state.throughput, active && use_rrt) *=
                            dr::exp(-dt * control);
                    }

                    if(dr::any_or<true>(act_spectral)) {
                        // Account for distance sampling weight in spectral mode
                        UnpolarizedSpectrum tr = dr::exp(-dt*residual_majorant);
                        Float pdf = index_spectrum<Float, Spectrum>(
                            dr::select(sampled, tr*residual_majorant, tr),
                            channel);
                        dr::masked(throughput, act_spectral)  *= tr / pdf;
                    }

                    if (dr::any_or<true>(sampled)) {
                        mei.t = maxt;
                        mei.p = state.ray(maxt);

                        // Retrieve scattering coefficients at position.
                        UnpolarizedSpectrum sigma_t, sigma_n;
                        std::tie( std::ignore, std::ignore, sigma_t) =
                            medium->get_scattering_coefficients(mei, sampled);
                        sigma_n = segment.majorant() - sigma_t;

                        if(dr::any_or<true>(act_spectral)) {
                            // Both RT and RRT spectral mode simplify to sigma_n.
                            dr::masked(throughput, sampled && act_spectral) *= sigma_n;
                        }

                        if(dr::any_or<true>(act_not_spectral)) {
                            dr::masked(throughput, sampled && act_not_spectral) *=
                                dr::maximum((1.f - (sigma_t - control)/residual_majorant), 0.f);
                        }

                        dr::masked(state.target_ot, sampled) =
                            -dr::log( 1.f - rng.template next_float<Float>(active) );
                    }
                    dr::masked(mei.t, !sampled) = dr::Infinity<Float>;
                    dr::masked(state.target_ot, !sampled && active) -= segment_ot;

                    // Never disable the loop, continue until maxt
                    Mask step = !sampled;
                    return std::pair<Mask,Mask>(step, active);
                }, active_medium);

                dr::masked(transmittance, active_medium) *= state.throughput;
            }

            // Handle interactions with surfaces
            Mask active_surface = si.is_valid() && active;
            if (dr::any_or<true>(active_surface)) {
                auto bsdf         = si.bsdf(ray);
                Spectrum bsdf_val = bsdf->eval_null_transmission(si, active_surface);
                bsdf_val = si.to_world_mueller(bsdf_val, si.wi, si.wi);
                dr::masked(transmittance, active_surface) *= bsdf_val;
            }

            // Update the ray with new origin & t parameter
            dr::masked(ray, active_surface) = si.spawn_ray(ray.d);
            ray.maxt = remaining_dist;

            // Continue tracing through scene if non-zero weights exist
            active &= active_surface &&
                      dr::any(unpolarized_spectrum(transmittance) != 0.f);

            // If a medium transition is taking place: Update the medium pointer
            Mask has_medium_trans = active_surface && si.is_medium_transition();
            if (dr::any_or<true>(has_medium_trans)) {
                dr::masked(medium, has_medium_trans) = si.target_medium(ray.d);
            }
        },
        "Volpath integrator emitter sampling with Extremum traversal & RRT.");

        return { ls.transmittance * emitter_val, dir_sample };
    }

    /// Outer overload: samples the emitter direction then calls the inner overload.
    template <typename Interaction>
    std::tuple<Spectrum, DirectionSample3f>
    sample_emitter(const Interaction &ref_interaction, const Scene *scene,
                   Sampler *sampler, MediumPtr medium,
                   UInt32 channel, Mask active) const {
        Spectrum transmittance(1.0f);

        auto [ds, emitter_val] = scene->sample_emitter_direction(ref_interaction, sampler->next_2d(active), false, active);
        return sample_emitter(ds, emitter_val, ref_interaction, scene, sampler, medium, channel, active);
    }


    //! @}
    // =============================================================

    std::string to_string() const override {
        return tfm::format("EOVolumetricPathIntegrator[\n"
                           "  max_depth = %i,\n"
                           "  rr_depth = %i\n"
                           "]",
                           m_max_depth, m_rr_depth);
    }

    Float mis_weight(Float pdf_a, Float pdf_b) const {
        pdf_a *= pdf_a;
        pdf_b *= pdf_b;
        Float w = pdf_a / (pdf_a + pdf_b);
        return dr::select(dr::isfinite(w), w, 0.f);
    };

    MI_DECLARE_CLASS(EOVolumetricPathIntegrator)
private:
    ScalarFloat m_rr_factor;

    bool m_enable_ddis;

    bool m_enable_pbs;
    ScalarUInt32 m_max_split_count;
    ScalarFloat m_min_split_threshold;
    ScalarFloat m_crit_rr_threshold;
    ScalarFloat m_min_rr_threshold;

    bool m_enable_nle;
    ScalarUInt32 m_first_clone_depth;
    ScalarUInt32 m_max_clone_depth;
    ScalarUInt32 m_nee_per_clone;
};

MI_EXPORT_PLUGIN(EOVolumetricPathIntegrator)
NAMESPACE_END(mitsuba)
