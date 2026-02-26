#include <mitsuba/core/plugin.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/render/eradiate/extremum.h>

NAMESPACE_BEGIN(mitsuba)

MI_VARIANT Medium<Float, Spectrum>::Medium()
    : JitObject<Medium>(""),
      m_is_homogeneous(false),
      m_has_spectral_extinction(true),
      m_has_local_extremum(false) {
}

MI_VARIANT Medium<Float, Spectrum>::Medium(const Properties &props)
    : JitObject<Medium>(props.id()) {

// #ERADIATE_CHANGE_BEGIN: Initialize extremum structure from properties
    for (auto &prop : props.objects()) {
        if (PhaseFunction *phase = prop.try_get<PhaseFunction>()) {
            if (m_phase_function)
                Throw("Only a single phase function can be specified per medium");
            m_phase_function = phase;
        }
        if (auto *extremum = prop.try_get<ExtremumStructure>()) {
            if (m_extremum_structure)
                Throw("Only a single extremum structure can be specified per medium");
            m_extremum_structure = extremum;
        }
    }

    m_has_local_extremum = m_extremum_structure != nullptr;
// #ERADIATE_CHANGE_END

    if (!m_phase_function) {
        // Create a default isotropic phase function
        m_phase_function =
            PluginManager::instance()->create_object<PhaseFunction>(Properties("isotropic"));
    }

    m_sample_emitters = props.get<bool>("sample_emitters", true);
}

MI_VARIANT Medium<Float, Spectrum>::~Medium() { }

MI_VARIANT void Medium<Float, Spectrum>::traverse(TraversalCallback *cb) {
    cb->put("phase_function", m_phase_function, ParamFlags::Differentiable);
// #ERADIATE_CHANGE_BEGIN: Traverse extremum structure
    cb->put("extremum_structure", m_extremum_structure, ParamFlags::NonDifferentiable);
// #ERADIATE_CHANGE_END
}

// #ERADIATE_CHANGE_BEGIN: Refactored for extremum structure support
MI_VARIANT
typename Medium<Float, Spectrum>::MediumInteraction3f
Medium<Float, Spectrum>::sample_interaction(const Ray3f &ray, Float sample,
                                            UInt32 channel, Mask active) const {
    MI_MASKED_FUNCTION(ProfilerPhase::MediumSample, active);

    // Initialize basic medium interaction fields
    MediumInteraction3f mei = dr::zeros<MediumInteraction3f>();
    mei.wi          = -ray.d;
    mei.sh_frame    = Frame3f(mei.wi);
    mei.time        = ray.time;
    mei.wavelengths = ray.wavelengths;
    mei.medium      = this;

    // Intersect AABB
    auto [aabb_its, mint, maxt] = intersect_aabb(ray);
    aabb_its &= (dr::isfinite(mint) || dr::isfinite(maxt));
    active &= aabb_its;
    dr::masked(mint, !active) = 0.f;
    dr::masked(maxt, !active) = dr::Infinity<Float>;

    mint = dr::maximum(0.f, mint);
    maxt = dr::minimum(ray.maxt, maxt);
    mei.mint = mint;

    Float target_od = -dr::log(1.f - sample);
    Float sampled_t;
    UnpolarizedSpectrum combined_extinction;

    if (m_has_local_extremum) {
        // Use extremum structure with local majorants
        auto [segment, tau_acc] = m_extremum_structure->sample_segment(ray, mint, maxt, target_od, active);
        sampled_t = segment.tmin +
                (target_od - tau_acc) / dr::maximum(segment.majorant, dr::Epsilon<Float>);
                
        Log(Debug, "Valid segment: %f, sampled_t: %f, maxt: %f", segment.valid(), sampled_t, maxt);
                
        // Store local majorant in combined_extinction
        combined_extinction[0] = segment.majorant;
        if constexpr (is_rgb_v<Spectrum>) {
            combined_extinction[1] = segment.majorant;
            combined_extinction[2] = segment.majorant;
        } else {
            DRJIT_MARK_USED(channel);
        }
    } else {
        // Traditional global majorant sampling
        combined_extinction = get_majorant(mei, active);
        Float m = combined_extinction[0];
        if constexpr (is_rgb_v<Spectrum>) {
            dr::masked(m, channel == 1u) = combined_extinction[1];
            dr::masked(m, channel == 2u) = combined_extinction[2];
        } else {
            DRJIT_MARK_USED(channel);
        }
        sampled_t = mint + (target_od / m);
    }

    // Finalize interaction
    Mask valid_mi = active && (sampled_t <= maxt);
    mei.t = dr::select(valid_mi, sampled_t, dr::Infinity<Float>);
    mei.p = ray(sampled_t);

    // TODO: with local extremum structures, this triggers a redundant evaluation
    // of the extremum grid. To fix this we probably need to change the medium interface.
    std::tie(mei.sigma_s, mei.sigma_n, mei.sigma_t) =
        get_scattering_coefficients(mei, valid_mi);
    mei.combined_extinction = combined_extinction;
    
    return mei;
}
// #ERADIATE_CHANGE_END

MI_VARIANT
std::pair<typename Medium<Float, Spectrum>::UnpolarizedSpectrum,
          typename Medium<Float, Spectrum>::UnpolarizedSpectrum>
Medium<Float, Spectrum>::transmittance_eval_pdf(const MediumInteraction3f &mi,
                                                const SurfaceInteraction3f &si,
                                                Mask active) const {
    MI_MASKED_FUNCTION(ProfilerPhase::MediumEvaluate, active);

    Float t      = dr::minimum(mi.t, si.t) - mi.mint;
    UnpolarizedSpectrum tr  = dr::exp(-t * mi.combined_extinction);
    UnpolarizedSpectrum pdf = dr::select(si.t < mi.t, tr, tr * mi.combined_extinction);
    return { tr, pdf };
}

// #RAY_CHANGE_BEGIN, NM 05/06/2024 : add function that calculates the transmittance and pdf  
MI_VARIANT
std::tuple<typename Medium<Float, Spectrum>::MediumInteraction3f, Float, Float>
Medium<Float, Spectrum>::sample_interaction_real(const Ray3f &/*ray*/, 
                                            const SurfaceInteraction3f &/*si*/, Float /*sample*/,
                                            UInt32 /*channel*/, Mask /*active*/) const {
    // PiecewiseVolPathIntegrator should only be used with piecewise medium                                                
    NotImplementedError("sample_interaction_real");
    MediumInteraction3f mi = dr::zeros<MediumInteraction3f>(); 
    return {mi,0.f,0.f};
}

MI_VARIANT
std::tuple<Float, Float, typename Medium<Float, Spectrum>::Mask>
Medium<Float, Spectrum>::eval_transmittance_pdf_real(const Ray3f &/*ray*/, 
                                    const SurfaceInteraction3f &/*si*/,
                                    UInt32 /*channel*/, Mask /*active*/) const {
    // PiecewiseVolPathIntegrator should only be used with piecewise medium
    NotImplementedError("eval_transmittance_pdf_real");
    return {0.f, 0.f, false};
}
// #RAY_CHANGE_END

MI_IMPLEMENT_TRAVERSE_CB(Medium, Object)
MI_INSTANTIATE_CLASS(Medium)
NAMESPACE_END(mitsuba)
