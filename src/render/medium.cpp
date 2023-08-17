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

    // Initialize basic medium interaction fields
    auto [mei, mint, maxt, active] = prepare_interaction_sampling(ray, _active);

    const Float desired_tau = -dr::log(1.f - sample);
    Log(Info, "desired_tau = %s", desired_tau);
    Float sampled_t;
    if (m_majorant_grid) {
        // --- Medium has a majorant supergrid
        sampled_t = m_majorant_grid->traverse_majorant_grid(
            desired_tau, ray, mint, maxt, active);
    } else {
        // --- Medium only has global majorant for the whole volume
        mei.combined_extinction = dr::detach(get_majorant(mei, active));
        Float m   = extract_channel(mei.combined_extinction, channel);
        sampled_t = mint + (desired_tau / m);
    }
    Log(Info, "sampled_t = %s", sampled_t);

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
