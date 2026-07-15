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
