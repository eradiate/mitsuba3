#include <mitsuba/core/properties.h>
#include <mitsuba/render/eradiate/extremum.h>

NAMESPACE_BEGIN(mitsuba)

MI_VARIANT ExtremumStructure<Float, Spectrum>::ExtremumStructure()
    : JitObject<ExtremumStructure>(""), m_scale(1.f) {
}

MI_VARIANT ExtremumStructure<Float, Spectrum>::ExtremumStructure(const Properties &props)
    : JitObject<ExtremumStructure>(props.id()), m_scale(1.f) {
}

MI_VARIANT ExtremumStructure<Float, Spectrum>::~ExtremumStructure() {
}

MI_VARIANT void ExtremumStructure<Float, Spectrum>::update_extremum(
    const ScalarBoundingBox3f &bbox, const Volume *volume,
    std::optional<ScalarFloat> scale) {
    // set scale if provided
    if (scale)
        set_scale(scale.value());
    // set validity bbox
    set_bbox(bbox);

    if (!m_bbox.valid())
        Throw("ExtremumStructure::update_extremum() called with an invalid bbox.");

    // rebuild the extremum structure
    build(volume);
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

MI_INSTANTIATE_CLASS(ExtremumStructure)
NAMESPACE_END(mitsuba)
