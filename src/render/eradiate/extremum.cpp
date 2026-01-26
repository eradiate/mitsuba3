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

MI_INSTANTIATE_CLASS(ExtremumStructure)
NAMESPACE_END(mitsuba)