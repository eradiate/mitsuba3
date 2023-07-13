#include <mitsuba/core/plugin.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/texture.h>

NAMESPACE_BEGIN(mitsuba)

MI_VARIANT Medium<Float, Spectrum>::Medium() : m_is_homogeneous(false), m_has_spectral_extinction(true) {}

MI_VARIANT Medium<Float, Spectrum>::Medium(const Properties &props)
: m_majorant_resolution_factor{ 0, 0, 0 }, m_majorant_grid(nullptr),
  m_majorant_factor(1.f), m_id(props.id()) {

    for (auto &[name, obj] : props.objects(false)) {
        auto *phase = dynamic_cast<PhaseFunction *>(obj.get());
        if (phase) {
            if (m_phase_function)
                Throw("Only a single phase function can be specified per medium");
            m_phase_function = phase;
            props.mark_queried(name);
        }
    }
    if (!m_phase_function) {
        // Create a default isotropic phase function
        m_phase_function =
            PluginManager::instance()->create_object<PhaseFunction>(Properties("isotropic"));
    }

    m_sample_emitters = props.get<bool>("sample_emitters", true);
    dr::set_attr(this, "use_emitter_sampling", m_sample_emitters);
    dr::set_attr(this, "phase_function", m_phase_function.get());
}

MI_VARIANT Medium<Float, Spectrum>::~Medium() {}

MI_VARIANT void Medium<Float, Spectrum>::traverse(TraversalCallback *callback) {
    callback->put_object("phase_function", m_phase_function.get(), +ParamFlags::Differentiable);

    if (m_majorant_grid)
        callback->put_object("majorant_grid", m_majorant_grid.get(), +ParamFlags::Differentiable);
}

MI_VARIANT
MI_INLINE
std::tuple<typename Medium<Float, Spectrum>::MediumInteraction3f, Float, Float,
           typename Medium<Float, Spectrum>::Mask>
Medium<Float, Spectrum>::prepare_interaction_sampling(const Ray3f &ray,
                                                      Mask active) const {
    // Initialize basic medium interaction fields
    MediumInteraction3f mei = dr::zeros<MediumInteraction3f>();
    mei.wi                  = -ray.d;
    mei.sh_frame            = Frame3f(mei.wi);
    mei.time                = ray.time;
    mei.wavelengths         = ray.wavelengths;
    mei.medium              = this;

    auto [aabb_its, mint, maxt] = intersect_aabb(ray);
    aabb_its &= (dr::isfinite(mint) || dr::isfinite(maxt));
    active &= aabb_its;
    dr::masked(mint, !active) = 0.f;
    dr::masked(maxt, !active) = dr::Infinity<Float>;

    mint = dr::maximum(0.f, mint);
    maxt = dr::minimum(ray.maxt, maxt);
    mei.mint   = mint;

    return std::make_tuple(mei, mint, maxt, active);
}

MI_VARIANT
typename Medium<Float, Spectrum>::MediumInteraction3f
Medium<Float, Spectrum>::sample_interaction(const Ray3f &ray, Float sample,
                                            UInt32 channel, Mask _active) const {
    MI_MASKED_FUNCTION(ProfilerPhase::MediumSample, _active);

    // initialize basic medium interaction fields
    auto [mei, mint, maxt, active] = prepare_interaction_sampling(ray, _active);

    const Float desired_tau = -dr::log(1.f - sample);
    Float sampled_t;
    if (m_majorant_grid) {
        // --- Spatially-variying majorant (supergrid).
        // 1. Prepare for DDA traversal
        // Adapted from:
        // https://github.com/francisengelmann/fast_voxel_traversal/blob/9664f0bde1943e69dbd1942f95efc31901fbbd42/main.cpp
        // TODO: allow precomputing all this (but be careful when ray origin is
        // updated)
        auto [dda_t, dda_tmax, dda_tdelta] =
            m_majorant_grid->prepare_majorant_grid_traversal(ray, mint, maxt, active);

        // 2. Traverse the medium with DDA until we reach the desired
        // optical depth.
        Mask active_dda = active;
        Mask reached    = false;
        Float tau_acc   = 0.f;
        dr::Loop<Mask> dda_loop("Medium::sample_interaction_dda");
        dda_loop.put(active_dda, reached, dda_t, dda_tmax, tau_acc, mei);
        dda_loop.init();
        while (dda_loop(dr::detach(active_dda))) {
            // Figure out which axis we hit first.
            // `t_next` is the ray's `t` parameter when hitting that axis.
            Float t_next = dr::min(dda_tmax);
            Vector3f tmax_update;
            Mask got_assigned = false;
            for (size_t k = 0; k < 3; ++k) {
                Mask active_k = dr::eq(dda_tmax[k], t_next);
                tmax_update[k] =
                    dr::select(!got_assigned && active_k, dda_tdelta[k], 0);
                got_assigned |= active_k;
            }

            // Lookup and accumulate majorant in current cell.
            dr::masked(mei.t, active_dda) = 0.5f * (dda_t + t_next);
            dr::masked(mei.p, active_dda) = ray(mei.t);
            // TODO: avoid this vcall, could lookup directly from the array
            // of floats (but we still need to account for the bbox, etc).
            Float majorant = m_majorant_grid->eval_1(mei, active_dda);
            Float tau_next = tau_acc + majorant * (t_next - dda_t);

            // For rays that will stop within this cell, figure out
            // the precise `t` parameter where `desired_tau` is reached.
            Float t_precise = dda_t + (desired_tau - tau_acc) / majorant;
            reached |= active_dda && (majorant > 0) && (t_precise < maxt) &&
                       (tau_next >= desired_tau);
            dr::masked(dda_t, active_dda) =
                dr::select(reached, t_precise, t_next);

            // Prepare for next iteration
            active_dda &= !reached && (t_next < maxt);
            dr::masked(dda_tmax, active_dda) = dda_tmax + tmax_update;
            dr::masked(tau_acc, active_dda)  = tau_next;
        }
        // Adopt the stopping location, making sure to convert to the main
        // ray's parametrization.
        sampled_t = dr::select(reached, dda_t, dr::Infinity<Float>);
    } else {
        // --- A single majorant for the whole volume.
        mei.combined_extinction = dr::detach(get_majorant(mei, active));
        Float m   = extract_channel(mei.combined_extinction, channel);
        sampled_t = mint + (desired_tau / m);
    }

    Mask valid_mei = active && (sampled_t <= maxt);
    mei.t          = dr::select(valid_mei, sampled_t, dr::Infinity<Float>);
    mei.p          = ray(sampled_t);

    if (m_majorant_grid) {
        // Otherwise it was already looked up above
        mei.combined_extinction =
            dr::detach(m_majorant_grid->eval_1(mei, valid_mei));
    }

    std::tie(mei.sigma_s, mei.sigma_n, mei.sigma_t) =
        get_scattering_coefficients(mei, valid_mei);
    return mei;
}

MI_VARIANT
std::pair<typename Medium<Float, Spectrum>::UnpolarizedSpectrum,
          typename Medium<Float, Spectrum>::UnpolarizedSpectrum>
Medium<Float, Spectrum>::eval_tr_and_pdf(const MediumInteraction3f &mi,
                                         const SurfaceInteraction3f &si,
                                         Mask active) const {
    MI_MASKED_FUNCTION(ProfilerPhase::MediumEvaluate, active);

    Float t      = dr::minimum(mi.t, si.t) - mi.mint;
    UnpolarizedSpectrum tr  = dr::exp(-t * mi.combined_extinction);
    UnpolarizedSpectrum pdf = dr::select(si.t < mi.t, tr, tr * mi.combined_extinction);
    return { tr, pdf };
}

MI_VARIANT
MI_INLINE Float
Medium<Float, Spectrum>::extract_channel(Spectrum value, UInt32 channel) {
    Float result = value[0];
    if constexpr (is_rgb_v<Spectrum>) { // Handle RGB rendering
        dr::masked(result, dr::eq(channel, 1u)) = value[1];
        dr::masked(result, dr::eq(channel, 2u)) = value[2];
    } else {
        DRJIT_MARK_USED(channel);
    }
    return result;
}

MI_IMPLEMENT_CLASS_VARIANT(Medium, Object, "medium")
MI_INSTANTIATE_CLASS(Medium)
NAMESPACE_END(mitsuba)
