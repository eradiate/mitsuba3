
#include <mitsuba/core/properties.h>
#include <mitsuba/render/phase.h>

NAMESPACE_BEGIN(mitsuba)

// #ERADIATE_CHANGE_BEGIN: DDIS
MI_VARIANT PhaseFunction<Float, Spectrum>::PhaseFunction(const Properties &props)
    : JitObject<PhaseFunction>(props.id()), 
      m_flags(+PhaseFunctionFlags::Empty),
      m_node_count(256) {
}
 
MI_VARIANT void PhaseFunction<Float, Spectrum>::get_nodes(
    std::vector<ScalarFloat> &nodes) const {
    
    nodes.resize(m_node_count);

    for (size_t i = 0; i < m_node_count; ++i) {
        nodes[i] = -1.f + 2.f * ScalarFloat(i) / ScalarFloat(m_node_count - 1);
    }
}

MI_VARIANT void PhaseFunction<Float, Spectrum>::eval_max(
    const std::vector<ScalarFloat> &nodes, 
    std::vector<ScalarFloat> &values) const {

    MI_IMPORT_TYPES(PhaseFunctionContext)
    
    size_t n = nodes.size();
    values.resize(n);

    MediumInteraction3f mi = dr::zeros<MediumInteraction3f>();
    mi.wi                  = Vector3f(0.f, 0.f, 1.f);
    PhaseFunctionContext ctx(nullptr);

    for (size_t i = 0; i < n; ++i) {
        ScalarFloat mu        = nodes[i];
        ScalarFloat sin_theta = dr::sqrt(dr::maximum(ScalarFloat(0), 1.f - mu * mu));
        Vector3f wo(sin_theta, 0.f, -mu);
        auto [val, pdf] = eval_pdf(ctx, mi, wo);
        values[i]       = dr::maximum(values[i], (ScalarFloat) dr::slice(unpolarized_spectrum(val)[0]));
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
