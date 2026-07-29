#include <mitsuba/core/properties.h>
#include <mitsuba/render/eradiate/extremum.h>
#include <mitsuba/render/medium.h>
#include <drjit/while_loop.h>

NAMESPACE_BEGIN(mitsuba)

MI_VARIANT ExtremumStructure<Float, Spectrum>::ExtremumStructure()
    : JitObject<ExtremumStructure>(""){
}

MI_VARIANT ExtremumStructure<Float, Spectrum>::ExtremumStructure(const Properties &props)
    : JitObject<ExtremumStructure>(props.id()){
}

MI_VARIANT ExtremumStructure<Float, Spectrum>::~ExtremumStructure() {
}

MI_VARIANT
TrackingState<Float, Spectrum>
ExtremumStructure<Float, Spectrum>::traverse_extremum(
    const Ray3f &ray,
    Float mint,
    Float maxt,
    UInt32 channel,
    TrackingStateType state,
    TrackingFunctionType *func,
    Mask active
) const {
    active &= maxt > mint;

    struct LoopState {
        ExtremumSegment segment;
        TrackingStateType state;
        Mask advance;
        Mask active;
        Float t;

        DRJIT_STRUCT(LoopState, segment, state, advance, active, t)
    } ls = {
        dr::zeros<ExtremumSegment>(),
        state,
        /*advance=*/active,
        active,
        mint
    };

    dr::tie(ls) = dr::while_loop(
        dr::make_tuple(ls),
        [](const LoopState &ls) { return ls.active; },
        [this, func, ray, maxt, channel](LoopState &ls) {
            if (dr::any_or<true>(ls.advance)) {
                ls.segment = next_segment(ray, ls.t, ls.active);
                dr::masked(ls.segment.maxt, ls.active) =
                    dr::minimum(ls.segment.maxt, maxt);
            }

            std::tie(ls.advance, ls.active) =
                func(ls.segment, ls.state, channel, ls.active);

            dr::masked(ls.t, ls.advance) = ls.segment.maxt;
            ls.active &= ls.t < maxt;
        },
        "Generic Extremum Traversal");

    return ls.state;
}

MI_VARIANT void ExtremumStructure<Float, Spectrum>::build(
    const ScalarBoundingBox3f & /*domain*/,
    const Volume * /*volume*/,
    ScalarFloat /*scale*/
) {
    NotImplementedError("build");
}

MI_VARIANT
typename ExtremumStructure<Float, Spectrum>::ExtremumSegment
ExtremumStructure<Float, Spectrum>::next_segment(
    const Ray3f & /*ray*/,
    Float /*t*/,
    Mask /*active*/
) const {
    NotImplementedError("next_segment");
}

MI_VARIANT void ExtremumStructure<Float, Spectrum>::set_domain(
    const ScalarBoundingBox3f &domain, const Volume *sigma_t
) {
    if (m_built_volume &&
        (m_built_volume != sigma_t || dr::any(m_bbox.min != domain.min) ||
         dr::any(m_bbox.max != domain.max)))
        Throw("ExtremumStructure: build() called again with a different "
              "domain or volume — one extremum structure belongs to exactly "
              "one medium.");
    m_built_volume = sigma_t;
    m_bbox = domain;
}

MI_INSTANTIATE_CLASS(ExtremumStructure)
NAMESPACE_END(mitsuba)
