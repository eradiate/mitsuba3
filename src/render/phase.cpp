
#include <mitsuba/core/properties.h>
#include <mitsuba/render/phase.h>

NAMESPACE_BEGIN(mitsuba)

// #ERADIATE_CHANGE_BEGIN: DDIS
MI_VARIANT PhaseFunction<Float, Spectrum>::PhaseFunction(const Properties &props)
    : JitObject<PhaseFunction>(props.id()), 
      m_flags(+PhaseFunctionFlags::Empty),
      m_node_count(256) {
}   

MI_VARIANT typename PhaseFunction<Float, Spectrum>::FloatStorage
PhaseFunction<Float, Spectrum>::get_envelope_nodes() const {
    FloatStorage i = dr::arange<FloatStorage>(m_node_count);
    return -1.f + 2.f * i / ScalarFloat(m_node_count - 1);
}

MI_VARIANT void
PhaseFunction<Float, Spectrum>::accumulate_envelope(const FloatStorage &nodes,
                                         FloatStorage &values) const {
    PhaseFunctionContext ctx(nullptr);
    size_t n               = dr::width(nodes);
    MediumInteraction3f mi = dr::zeros<MediumInteraction3f>(n);
    mi.wi                  = Vector3f(0.f, 0.f, 1.f);

    if constexpr (dr::is_jit_v<Float>) {
        Float mu        = nodes;
        Float sin_theta = dr::sqrt(dr::maximum(Float(0), 1.f - mu * mu));
        Vector3f wo(sin_theta, 0.f, -mu);
        auto [val, pdf] = eval_pdf(ctx, mi, wo);
        values          = dr::maximum(values, unpolarized_spectrum(val)[0]);

    } else {
        for (size_t i = 0; i < n; ++i) {
            ScalarFloat mu        = nodes[i];
            ScalarFloat sin_theta = dr::sqrt(dr::maximum(ScalarFloat(0), 1.f - mu * mu));
            Vector3f wo(sin_theta, 0.f, -mu);
            auto [val, pdf] = eval_pdf(ctx, mi, wo);
            values[i]       = dr::maximum(values[i], (ScalarFloat) unpolarized_spectrum(val)[0]);
        }
    }
}
// #ERADIATE_CHANGE_END

// #ERADIATE_CHANGE_BEGIN: DDIS phase dirty flag
MI_VARIANT void PhaseFunction<Float, Spectrum>::parameters_changed(
        const std::vector<std::string> &/*keys*/) {
    m_dirty = true;
}
// #ERADIATE_CHANGE_END
MI_INSTANTIATE_CLASS(PhaseFunction)
NAMESPACE_END(mitsuba)
