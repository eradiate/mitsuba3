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

 * - ddis_threshold
   - |float|
   - Specifies the probability to importance sample the phase using the emitter as 
     incident direction. Set to a negative value to disable. (Default: 0.1)

 * - hide_emitters
   - |bool|
   - Hide directly visible emitters. (Default: no, i.e. |false|)

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

In addition, it implements the DDIS variance reduction method (Buras and Mayer, 2011) which reduces 
noise in the presence of strongly peaked phase functions. More variance reduction methods coming up
soon! 

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
                    MediumPtr, PhaseFunctionContext, ExtremumStructure)

    using TrackingState = TrackingState<Float, Spectrum>;

    EOVolumetricPathIntegrator(const Properties &props) : Base(props) {
        m_ddis_threshold = props.get<ScalarFloat>("ddis_threshold", 0.1f);
        m_rr_factor = props.get<ScalarFloat>("rr_factor", 0.97f);

        if (m_ddis_threshold > 1.f)
            Throw("`ddis_threshold` is larger than 1.");

        if (m_rr_factor < 0.f || m_rr_factor > 1.f)
            Throw("`rr_factor` is outside range [0.,1.].");
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
        MediumInteraction3f mei = dr::zeros<MediumInteraction3f>();
        Mask specular_chain = active && !m_hide_emitters;
        UInt32 depth = 0;

        UInt32 channel = 0;
        if (is_rgb_v<Spectrum>) {
            uint32_t n_channels = (uint32_t) dr::size_v<Spectrum>;
            channel = (UInt32) dr::minimum(sampler->next_1d(active) * n_channels, n_channels - 1);
        }

        SurfaceInteraction3f si = dr::zeros<SurfaceInteraction3f>();
        Interaction3f last_scatter_event = dr::zeros<Interaction3f>();
        Float last_scatter_direction_pdf = 1.f;

        /* Set up a Dr.Jit loop (optimizes away to a normal loop in scalar mode,
           generates wavefront or megakernel renderer based on configuration).
           Register everything that changes as part of the loop here */
        struct LoopState {
            Mask active;
            UInt32 depth;
            Ray3f ray;
            Spectrum throughput;
            Spectrum result;
            SurfaceInteraction3f si;
            MediumInteraction3f mei;
            MediumPtr medium;
            Float eta;
            Interaction3f last_scatter_event;
            Float last_scatter_direction_pdf;
            Mask specular_chain;
            Mask valid_ray;
            Sampler* sampler;

            DRJIT_STRUCT(LoopState, active, depth, ray, throughput, result, \
                si, mei, medium, eta, last_scatter_event, \
                last_scatter_direction_pdf, \
                specular_chain, valid_ray, sampler)
        } ls = {
            active,
            depth,
            ray,
            throughput,
            result,
            si,
            mei,
            medium,
            eta,
            last_scatter_event,
            last_scatter_direction_pdf,
            specular_chain,
            valid_ray,
            sampler
        };

        dr::tie(ls) = dr::while_loop(dr::make_tuple(ls),
            [](const LoopState& ls) { return ls.active; },
            [this, scene, channel](LoopState& ls) {

            Mask& active = ls.active;
            UInt32& depth = ls.depth;
            Ray3f& ray = ls.ray;
            Spectrum& throughput = ls.throughput;
            Spectrum& result = ls.result;
            SurfaceInteraction3f& si = ls.si;
            MediumInteraction3f& mei = ls.mei;
            MediumPtr& medium = ls.medium;
            Float& eta = ls.eta;
            Interaction3f& last_scatter_event = ls.last_scatter_event;
            Float& last_scatter_direction_pdf = ls.last_scatter_direction_pdf;
            Mask& specular_chain = ls.specular_chain;
            Mask& valid_ray = ls.valid_ray;
            Sampler* sampler = ls.sampler;

            // ----------------- Handle termination of paths ------------------
            // Russian roulette: try to keep path weights equal to one, while accounting for the
            // solid angle compression at refractive index boundaries. Stop with at least some
            // probability to avoid getting stuck (e.g. due to total internal reflection)
            active &= dr::any(unpolarized_spectrum(throughput) != 0.f);
            Float q = dr::minimum(dr::max(unpolarized_spectrum(throughput)) * dr::square(eta), m_rr_factor);
            Mask perform_rr = (depth > (uint32_t) m_rr_depth);
            active &= sampler->next_1d(active) < q || !perform_rr;
            dr::masked(throughput, perform_rr) *= dr::rcp(dr::detach(q));
            
            active &= depth < (uint32_t) m_max_depth;
            if (dr::none_or<false>(active))
                return;

            // ----------------------- Sampling the RTE -----------------------
            Mask active_medium  = active && (medium != nullptr);
            Mask active_surface = active;

            Mask act_null_scatter = false, act_medium_scatter = false,
                 escaped_medium = false;

            dr::masked(si, active) = scene->ray_intersect(ray, active);
            dr::masked(ray.maxt, active) = si.t; 
            active_surface &= si.is_valid();

            // If the medium does not have a spectrally varying extinction,
            // we can perform a few optimizations to speed up rendering
            Mask is_spectral = active_medium;
            Mask not_spectral = false;
            if (dr::any_or<true>(active_medium)) {
                is_spectral &= medium->has_spectral_extinction();
                not_spectral = !is_spectral && active_medium;
            }

            if (dr::any_or<true>(active_medium)) {
                // Prepare Extremum traversal
                auto extremum = medium->extremum_structure();
                Float mint, maxt;
                std::tie(mei, mint, maxt) = 
                    medium->prepare_medium_traversal(ray, active_medium);

                Float sample1 = sampler->next_1d();
                Float sample2 = sampler->next_1d();
                auto [seed, offset] = new_seed_offset(sample1, sample2);
                PCG32<UInt32> rng;
                rng.seed(seed, offset);

                TrackingState state {
                    ray,
                    rng,
                    mei,
                    /*throughput=*/UnpolarizedSpectrum(1.f),
                };

                // Traverse extremum segments and perform delta tracking
                state = extremum->traverse_extremum(
                ray, mint, maxt, channel, state,
                [](const ExtremumSegment& segment, TrackingState& state, 
                   const UInt32& channel, Mask active) {

                    UnpolarizedSpectrum &throughput = state.throughput;
                    PCG32<UInt32> &rng              = state.rng;
                    MediumInteraction3f& mei        = state.mei; 
                    
                    MediumPtr medium      = mei.medium;
                    Mask act_spectral     = medium->has_spectral_extinction() && active;
                    Mask act_not_spectral = !medium->has_spectral_extinction() && active;

                    // Check if the last iteration had a valid interaction within the segment
                    Float mint = dr::select(
                        mei.is_valid(),
                        dr::maximum(segment.mint, mei.t), 
                        segment.mint
                    );
                    
                    // Sample a new potential interaction
                    Float target_ot = 
                        -dr::log( 1.f - state.rng.template next_float<Float>(active) );
                    Float sampled_t =
                        mint
                        + (target_ot) / dr::maximum(segment.majorant(), dr::Epsilon<Float>);
                    
                    Mask sampled = (sampled_t < segment.maxt) && active;
                    Float dt = dr::select(sampled, sampled_t - mint, segment.maxt - mint);

                    if( dr::any_or<true>(act_spectral) ) {
                        // Accumulate transmittance in the throughput and pdf (spectral only).
                        UnpolarizedSpectrum tr = dr::exp(-dt*segment.majorant());
                        Float pdf = index_spectrum<Float, Spectrum>(
                            dr::select(sampled, tr*segment.majorant(), tr), 
                            channel);
                        dr::masked(throughput, act_spectral)  *= tr/pdf;
                    } 

                    if (dr::any_or<true>(sampled)) {
                        mei.t = sampled_t;
                        mei.p = state.ray(sampled_t);

                        // Retrieve scattering coefficients at position.
                        UnpolarizedSpectrum sigma_s, sigma_n, sigma_t;
                        std::tie(sigma_s, std::ignore, sigma_t) = 
                            medium->get_scattering_coefficients(mei, sampled);
                        sigma_n = segment.majorant() - sigma_t;

                        // Sample event type
                        Float null_scatter_prob = dr::mean(sigma_n / segment.majorant());
                        Mask null_scatter = 
                            (rng.template next_float<Float>(sampled) < null_scatter_prob) && sampled;
                        Mask real_scatter = !null_scatter && sampled;
                        
                        // Accumulate throughput and pdf given the event type and is_spectral.
                        if (dr::any_or<true>(null_scatter && act_spectral)) {
                            dr::masked(throughput, null_scatter && act_spectral) *= 
                                sigma_n / null_scatter_prob;
                        }

                        if (dr::any_or<true>(real_scatter)) {
                            if(dr::any_or<true>(act_spectral)) {
                                dr::masked(throughput, real_scatter && act_spectral) *= 
                                    sigma_s / (1.0f - null_scatter_prob);
                            }

                            if(dr::any_or<true>(act_not_spectral)) {
                                // pdf = sigma_t / sigma_maj, sigma_maj gets cancelled from sigma_s/sigma_maj 
                                dr::masked(throughput, real_scatter && act_not_spectral) *= 
                                    sigma_s / sigma_t;
                            }

                            // disable the loop once we encounter a real scattering interaction
                            active &= !real_scatter;
                        }
                    }

                    dr::masked(mei.t, !sampled) = dr::Infinity<Float>;

                    Mask step = !sampled;
                    return std::pair<Mask, Mask>(step, active);
                });

                // Update throughput by the transmittance and pdf weight
                dr::masked(throughput, active_medium) *= state.throughput;
                dr::masked(mei, active_medium) = state.mei; 

                escaped_medium = active_medium && !mei.is_valid();
                active_medium &= mei.is_valid();

                act_medium_scatter = !escaped_medium && active_medium;
                dr::masked(depth, act_medium_scatter) += 1;
                dr::masked(last_scatter_event, act_medium_scatter) = mei;
            }

            // Dont estimate lighting if we exceeded number of bounces
            active &= depth < (uint32_t) m_max_depth;
            act_medium_scatter &= active;

            if (dr::any_or<true>(act_medium_scatter)) {

                PhaseFunctionContext phase_ctx(sampler);
                auto phase = mei.medium->phase_function();

                // --------------------- Emitter sampling ---------------------
                Mask sample_emitters = mei.medium->use_emitter_sampling();
                valid_ray |= act_medium_scatter;
                specular_chain &= !act_medium_scatter;
                specular_chain |= act_medium_scatter && !sample_emitters;

                Mask active_e = act_medium_scatter && sample_emitters;

                Mask perform_ddis = m_ddis_threshold > 0.f && active_e;
                MediumInteraction3f ddis_mei = mei;

                if (dr::any_or<true>(active_e)) {
                    auto [emitted, ds] = sample_emitter(mei, scene, sampler, medium, channel, active_e);
                    auto [phase_val, phase_pdf] = phase->eval_pdf(phase_ctx, mei, ds.d, active_e);

                    if (dr::any_or<true>(perform_ddis)) {
                        // Set ddis incoming direction towards emitter
                        dr::masked(ddis_mei.wi, perform_ddis) = -ds.d;
                        dr::masked(ddis_mei.sh_frame, perform_ddis) = Frame3f(-ds.d);
                        
                        Float ddis_phase_pdf;
                        std::tie(std::ignore, ddis_phase_pdf) = 
                            phase->eval_pdf(phase_ctx, ddis_mei, ds.d, perform_ddis);

                        // Compute mixture pdf accounting for DDIS probability
                        dr::masked(phase_pdf, perform_ddis) =
                            (1.f - m_ddis_threshold) * phase_pdf + m_ddis_threshold * ddis_phase_pdf;
                    }

                    dr::masked(result, active_e) += throughput * phase_val * emitted *
                                                    mis_weight(ds.pdf, dr::select(ds.delta, 0.f, phase_pdf));
                }

                // ------------------ Phase function sampling -----------------
                dr::masked(phase, !act_medium_scatter) = nullptr;

                Float eps = sampler->next_1d(active_medium);
                Mask active_ddis = (eps < m_ddis_threshold) && act_medium_scatter && perform_ddis;
                MediumInteraction3f sample_mei = dr::select(active_ddis, ddis_mei, mei);

                // ddis off -> phase_weight: p(mu)/pdf(mu); phase_pdf: pdf(mu);
                // ddis on  -> phase_weight: p(mu')/pdf(mu'); phase_pdf: pdf(mu');
                auto [wo, phase_weight, phase_pdf] = phase->sample(phase_ctx, sample_mei,
                    sampler->next_1d(act_medium_scatter),
                    sampler->next_2d(act_medium_scatter),
                    act_medium_scatter);
                act_medium_scatter &= phase_pdf > 0.f;

                if (dr::any_or<true>(perform_ddis)) {
                    auto [natural_val, natural_pdf] = phase->eval_pdf(phase_ctx, mei, wo, perform_ddis);
                    auto [ddis_val, ddis_pdf] = phase->eval_pdf(phase_ctx, ddis_mei, wo, perform_ddis);

                    dr::masked(phase_pdf, perform_ddis) =
                        ((1.f - m_ddis_threshold) * natural_pdf + m_ddis_threshold * ddis_pdf);
                    dr::masked(phase_weight, perform_ddis) = natural_val / phase_pdf;
                }

                Ray3f new_ray  = mei.spawn_ray(wo);
                dr::masked(ray, act_medium_scatter) = new_ray;
                dr::masked(last_scatter_direction_pdf, act_medium_scatter) = phase_pdf;
                dr::masked(throughput, act_medium_scatter) *= phase_weight;
            }

            // --------------------- Surface Interactions ---------------------
            // Interact with the surface only if we haven't interacted with the 
            // medium before and there is a valid intersection (already accounted).
            active_surface &= !act_medium_scatter;
            
            if (dr::any_or<true>(active_surface)) {
                // ---------------------- Hide area emitters ----------------------
                if (m_hide_emitters && dr::any_or<true>(ls.depth == 0u)) {
                    // Are we on the first segment and did we hit an area emitter?
                    // If so, skip all area emitters along this ray
                    Mask skip_emitters = si.is_valid() &&
                                         (si.shape->emitter() != nullptr) &&
                                         (ls.depth == 0);

                    if (dr::any_or<true>(skip_emitters)) {
                        Ray3f ray = si.spawn_ray(ls.ray.d);
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
            if (dr::any_or<true>(active_surface)) {
                // --------------------- Emitter sampling ---------------------
                BSDFContext ctx;
                BSDFPtr bsdf  = si.bsdf(ray);
                Mask active_e = active_surface && has_flag(bsdf->flags(), BSDFFlags::Smooth) && (depth + 1 < (uint32_t) m_max_depth);

                if (likely(dr::any_or<true>(active_e))) {
                    // auto [emitted, ds] = sample_emitter(si, scene, sampler, medium, channel, active_e);
                    auto [emitted, ds] = sample_emitter(si, scene, sampler, medium, channel, active_e);

                    // Query the BSDF for that emitter-sampled direction
                    Vector3f wo       = si.to_local(ds.d);
                    Spectrum bsdf_val = bsdf->eval(ctx, si, wo, active_e);
                    bsdf_val = si.to_world_mueller(bsdf_val, -wo, si.wi);

                    // Determine probability of having sampled that same
                    // direction using BSDF sampling.
                    Float bsdf_pdf = bsdf->pdf(ctx, si, wo, active_e);
                    result[active_e] += throughput * bsdf_val * mis_weight(ds.pdf, dr::select(ds.delta, 0.f, bsdf_pdf)) * emitted;
                }

                // ----------------------- BSDF sampling ----------------------
                auto [bs, bsdf_val] = bsdf->sample(ctx, si, sampler->next_1d(active_surface),
                                                   sampler->next_2d(active_surface), active_surface);
                bsdf_val = si.to_world_mueller(bsdf_val, -bs.wo, si.wi);

                dr::masked(throughput, active_surface) *= bsdf_val;
                dr::masked(eta, active_surface) *= bs.eta;

                Ray3f bsdf_ray                  = si.spawn_ray(si.to_world(bs.wo));
                dr::masked(ray, active_surface) = bsdf_ray;

                Mask non_null_bsdf = active_surface && !has_flag(bs.sampled_type, BSDFFlags::Null);
                dr::masked(depth, non_null_bsdf) += 1;

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
            active &= (active_surface | active_medium);
        },
        "Volpath integrator");

        return { ls.result, ls.valid_ray };
    }

    /// Samples an emitter in the scene and evaluates its attenuated contribution
    template <typename Interaction>
    std::tuple<Spectrum, DirectionSample3f>
    sample_emitter(const Interaction &ref_interaction, const Scene *scene,
                   Sampler *sampler, MediumPtr medium,
                   UInt32 channel, Mask active) const {
        Spectrum transmittance(1.0f);

        auto [ds, emitter_val] = scene->sample_emitter_direction(ref_interaction, sampler->next_2d(active), false, active);
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
            
            if (dr::any_or<true>(active_medium)) {
                // Prepare extremum traversal
                auto extremum = medium->extremum_structure();
                auto [mei, mint, maxt] = 
                    medium->prepare_medium_traversal(ray, active_medium);

                Float sample1 = sampler->next_1d();
                Float sample2 = sampler->next_1d();
                auto [seed, offset] = new_seed_offset(sample1, sample2);
                PCG32<UInt32> rng;
                rng.seed(seed, offset);

                TrackingState state {
                    ray,
                    rng,
                    mei,
                    /*throughput=*/UnpolarizedSpectrum(1.f),
                };

                // Unified Ratio Tracking and Residual Ratio Tracking approach:
                // if `use_rrt` is false, set control to 0., which automatically
                // devolves the algorithm to Ratio Tracking.
                state = extremum->traverse_extremum(
                ray, mint, maxt, channel, state,
                [](const ExtremumSegment& segment, TrackingState& state, 
                   const UInt32& channel, Mask active) {

                    UnpolarizedSpectrum &throughput = state.throughput;
                    PCG32<UInt32> &rng       = state.rng;
                    MediumInteraction3f& mei = state.mei; 
                    Mask use_rrt             = mei.medium->use_rrt();
                    MediumPtr medium         = mei.medium;
                    Mask act_spectral     = medium->has_spectral_extinction() && active;
                    Mask act_not_spectral = !medium->has_spectral_extinction() && active;
                    
                    Float control = dr::select(use_rrt, segment.minorant(), 0.f);
                    Float residual_majorant = segment.majorant() - control;

                    Float mint = dr::select(
                        mei.is_valid(),
                        dr::maximum(segment.mint, mei.t), 
                        segment.mint
                    );

                    Float target_ot = 
                        -dr::log( 1.f - rng.template next_float<Float>(active) );
                    Float sampled_t =
                        mint
                        + (target_ot) / dr::maximum(residual_majorant, dr::Epsilon<Float>);
                    
                    Mask sampled = (sampled_t < segment.maxt) && active;
                    Float dt = dr::select(sampled, sampled_t - mint, segment.maxt - mint);

                    if(dr::any_or<true>(use_rrt)) {
                        // We accumulate control transmittance every interaction
                        // instead of accumulating the optical thickness.
                        // Choice of simplicity vs performance.
                        dr::masked(throughput, active && use_rrt) *= dr::exp( -dt * control );
                    }

                    if(dr::any_or<true>(act_spectral)) {
                        // Account for distance sampling weight in spectral mode
                        UnpolarizedSpectrum tr = dr::exp(-dt*residual_majorant);
                        // dr::masked(throughput, act_spectral)  *= tr;
                        Float pdf = index_spectrum<Float, Spectrum>(
                            dr::select(sampled, tr*residual_majorant, tr), 
                            channel);
                        dr::masked(throughput, act_spectral)  *= tr / pdf;
                    } 

                    if (dr::any_or<true>(sampled)) {
                        mei.t = sampled_t;
                        mei.p = state.ray(sampled_t);

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
                    }
                    dr::masked(mei.t, !sampled) = dr::Infinity<Float>;
                    
                    // Never disable the loop, continue until maxt
                    Mask step = !sampled;
                    return std::pair<Mask,Mask>(step, active);
                });

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
        "Volpath integrator emitter sampling special edition");

        return { ls.transmittance * emitter_val, dir_sample };
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
    ScalarFloat m_ddis_threshold;
};

MI_EXPORT_PLUGIN(EOVolumetricPathIntegrator)
NAMESPACE_END(mitsuba)
