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
    TrackingState /*state*/,
    TrackingFunction * /*func*/,
    Mask /*active*/
) const {
    NotImplementedError("traverse_extremum");
}

MI_INSTANTIATE_CLASS(ExtremumStructure)
NAMESPACE_END(mitsuba)
