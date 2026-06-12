#include <mitsuba/core/properties.h>
#include <mitsuba/render/eradiate/extremum.h>

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
    const Ray3f &/*ray*/,
    Float /*mint*/,
    Float /*maxt*/,
    UInt32 /*channel*/,
    TrackingStateType /*state*/,
    TrackingFunctionType * /*func*/,
    Mask /*active*/
) const {
    NotImplementedError("traverse_extremum");
}

MI_VARIANT
auto ExtremumStructure<Float, Spectrum>::sample_segment(
    const Ray3f &ray,
    Float mint,
    Float maxt,
    Float target_ot,
    Mask active
) const -> std::pair<ExtremumSegment, Float> {
    using SegmentType = ExtremumSegment;

    // The tracking callback cannot capture, so results are smuggled through
    // repurposed TrackingState fields. This is a testing convenience only;
    // production tracking algorithms live in the integrators. Conventions:
    //   mei.mint / mei.t        -> sampled segment entry / exit distance
    //   mei.sigma_s             -> sampled segment minorant
    //   mei.combined_extinction -> sampled segment majorant
    //   throughput              -> accumulated optical thickness
    TrackingStateType state = dr::zeros<TrackingStateType>();
    state.ray       = ray;
    state.target_ot = target_ot;

    TrackingFunctionType *func =
        [](const SegmentType &segment, TrackingStateType &state,
           const UInt32 & /*channel*/, Mask active) -> std::pair<Mask, Mask> {
        Float ot_acc = state.throughput[0];
        Float ot_seg = (segment.maxt - segment.mint) * segment.majorant();
        Mask reached = active && (ot_acc + ot_seg > state.target_ot);

        dr::masked(state.throughput, active && !reached) += ot_seg;

        dr::masked(state.mei.mint, reached) = segment.mint;
        dr::masked(state.mei.t, reached)    = segment.maxt;
        dr::masked(state.mei.sigma_s, reached) =
            UnpolarizedSpectrum(segment.minorant());
        dr::masked(state.mei.combined_extinction, reached) =
            UnpolarizedSpectrum(segment.majorant());

        return { /*advance=*/Mask(true), active && !reached };
    };

    state = traverse_extremum(ray, mint, maxt, UInt32(0), state, func, active);

    SegmentType segment = dr::zeros<SegmentType>();
    Mask found = active && state.mei.is_valid();
    dr::masked(segment.mint, found)  = state.mei.mint;
    dr::masked(segment.maxt, found)  = state.mei.t;
    dr::masked(segment.value, found) = Vector2f(
        state.mei.sigma_s[0], state.mei.combined_extinction[0]);

    return { segment, state.throughput[0] };
}

MI_INSTANTIATE_CLASS(ExtremumStructure)
NAMESPACE_END(mitsuba)
