#include <mitsuba/core/properties.h>
#include <mitsuba/render/interaction.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/eradiate/extremum.h>
#include <mitsuba/render/eradiate/phase_utils.h>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _medium-overlap:

Overlapping medium (:monosp:`overlap`)
--------------------------------------

.. pluginparameters::

 * - aabb_min, aabb_max
   - |point|
   - Optional override to the medium bounding box. Defaults to the union of
     the components' domains.

 * - (Nested plugin)
   - |medium|
   - One or more nested medium function that constitutes a medium component.
   - |exposed|, |differentiable|

 * - sample_emitters
   - |bool|
   - Flag to specify whether shadow rays should be cast from inside the volume (Default: |true|)
     If the medium is enclosed in a :ref:`dielectric <bsdf-dielectric>` boundary,
     shadow rays are ineffective and turning them off will significantly reduce
     render time. This can reduce render time up to 50% when rendering objects
     with subsurface scattering.

 * - ddis_threshold
   - |float|
   - Specifies the probability to importance sample the phase using the emitter as
     incident direction. Set to a negative value to disable. (Default: 0.1)


This plugin provides an aggregate medium implementation for overlapping
media. It accepts multiple nested media components, making the assumption
that those media do not have a physical boundary, and that the data is not
defined outside of their bounding box. Scattering coefficients are the sums
of the components' coefficients. A tracking collision picks one component
proportionally to its majorant and evaluates that component alone; a real
collision scatters with that component's phase function.

No combined majorant structure is built: the medium traverses its
components' extremum structures side by side and sums their segments. Each
component must have a single extremum structure, so overlaps do not nest
(list every medium in one :monosp:`overlap` instead); a component may be a
:monosp:`repeat` of a plain medium. At most four components are supported.
*/
template <typename Float, typename Spectrum>
class OverlappingMedium final : public Medium<Float, Spectrum> {
public:
    MI_IMPORT_BASE(Medium, m_is_homogeneous, m_has_spectral_extinction,
                    m_phase_function, m_extrema,
                    m_ddis_phase_function, m_ddis_threshold,
                    create_ddis_phase_function
                )
    MI_IMPORT_TYPES(Scene, Sampler, ExtremumStructure,
                    PhaseFunction, PhaseFunctionPtr)

    using FloatStorage = DynamicBuffer<Float>;

    OverlappingMedium(const Properties &props) : Base(props) {
        m_has_spectral_extinction = false;
        m_is_homogeneous = true;

        for (auto &prop : props.objects()) {
            if (auto *component = prop.try_get<Base>()) {
                m_has_spectral_extinction |= component->has_spectral_extinction();
                m_is_homogeneous &= component->is_homogeneous();

                m_components.push_back(component);
            }
        }

        if(m_components.empty())
            Throw("Must have at least one medium component.");


        ScalarBoundingBox3f domain;
        std::vector<const PhaseFunction *> phase_functions;
        for (auto &component : m_components) {
            if (component->component_count() != 1)
                Throw("overlapping medium: component %s must have exactly one "
                      "extremum structure, has %zu (nested overlaps are not "
                      "supported)", component->to_string(),
                      component->component_count());
            m_extrema.push_back(component->extrema()[0]);
            domain.expand(m_extrema.back()->bbox());
            phase_functions.push_back(component->phase_function());
        }

        if (m_extrema.size() > MI_MAX_OVERLAPPING_MEDIA)
            Throw("overlapping medium: too many components (%zu > %zu); "
                  "raise MI_MAX_OVERLAPPING_MEDIA at build time",
                  m_extrema.size(), MI_MAX_OVERLAPPING_MEDIA);

        m_phase_functions_dr = dr::load<DynamicBuffer<PhaseFunctionPtr>>(
            phase_functions.data(), phase_functions.size());
        dr::eval(m_phase_functions_dr);

        m_phase_function = m_components[0]->phase_function();

        m_ddis_threshold = props.get<ScalarFloat>("ddis_threshold", 0.1f);

        if (m_ddis_threshold > 0.f) {
            m_ddis_phase_function = static_cast<PhaseFunction*>(create_ddis_phase_function());
        }

        // Medium bounding box: union of the components' domains by default
        if (props.has_property("aabb_min") && props.has_property("aabb_max")) {
            ScalarPoint3f aabb_min = props.get<ScalarPoint3f>("aabb_min");
            ScalarPoint3f aabb_max = props.get<ScalarPoint3f>("aabb_max");
            m_aabb = ScalarBoundingBox3f(aabb_min, aabb_max);
        } else {
            m_aabb = domain;
        }
    }

    void traverse(TraversalCallback *cb) override {
        for (size_t i = 0; i < m_components.size(); ++i) {
            cb->put("component" + std::to_string(i), m_components[i], ParamFlags::Differentiable);
        }
        if (m_ddis_phase_function != nullptr)
            cb->put("ddis_phase_function", m_ddis_phase_function, ParamFlags::Differentiable);
        Base::traverse(cb);
    }

    UnpolarizedSpectrum
    get_majorant(const MediumInteraction3f &mi, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::MediumEvaluate, active);
        UnpolarizedSpectrum majorant = m_components[0]->get_majorant(mi, active);
        for (size_t i = 1; i < m_components.size(); ++i){
            majorant += m_components[i]->get_majorant(mi, active);
        }
        return majorant;
    }

    std::tuple<UnpolarizedSpectrum, UnpolarizedSpectrum, UnpolarizedSpectrum>
    get_scattering_coefficients(const MediumInteraction3f &mi,
                                Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::MediumEvaluate, active);

        UnpolarizedSpectrum sigmas(0.f), sigmat(0.f), majorant(0.f);

        for (size_t i = 0; i < m_components.size(); ++i) {
            Mask accumulate = active && m_components[i]->in_aabb(mi.p);
            if (dr::any_or<true>(accumulate)) {
                auto [c_sigmas, c_sigman, c_sigmat] =
                    m_components[i]->get_scattering_coefficients(mi, accumulate);
                sigmas += c_sigmas;
                sigmat += c_sigmat;
                majorant += m_components[i]->get_majorant(mi, active);
            }
        }

        // Each component's own sigma_n is relative to *its* local majorant,
        // not the composite's combined one, so it cannot be summed directly.
        // Re-derive sigma_n from the combined majorant instead.
        UnpolarizedSpectrum sigman = majorant - sigmat;

        return { sigmas, sigman, sigmat };
    }

    std::pair<UnpolarizedSpectrum, UnpolarizedSpectrum>
    eval_component(const MediumInteraction3f &mei, const UInt32 &component,
                   Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::MediumEvaluate, active);

        UnpolarizedSpectrum sigma_s(0.f), sigma_t(0.f);
        for (size_t k = 0; k < m_components.size(); ++k) {
            Mask owned = active && component == (uint32_t) k &&
                         m_components[k]->in_aabb(mei.p);
            if (dr::any_or<true>(owned)) {
                auto [c_sigma_s, c_sigma_t] =
                    m_components[k]->eval_component(mei, 0u, owned);
                dr::masked(sigma_s, owned) = c_sigma_s;
                dr::masked(sigma_t, owned) = c_sigma_t;
            }
        }
        return { sigma_s, sigma_t };
    }

    std::pair<PhaseFunctionPtr, MediumInteraction3f>
    phase_function(const MediumInteraction3f &mei, const UInt32 &component,
                   Mask active) const override {
        MediumInteraction3f lookup = mei;
        for (size_t k = 0; k < m_components.size(); ++k) {
            Mask owned = active && component == (uint32_t) k;
            if (dr::any_or<true>(owned))
                dr::masked(lookup, owned) =
                    m_components[k]->phase_function(mei, 0u, owned).second;
        }
        return { dr::gather<PhaseFunctionPtr>(m_phase_functions_dr, component,
                                              active),
                 lookup };
    }

    std::tuple<Mask, Float, Float>
    intersect_aabb(const Ray3f &ray) const override {
        return m_aabb.ray_intersect(ray);
    }

    virtual Mask
    in_aabb(const Point3f &pos) const override {
        return m_aabb.contains(pos);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "OverlappingMedium[" << std::endl
            << "  ddis_phase_function   = " << string::indent(m_ddis_phase_function) << "," << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS(OverlappingMedium)


    void update_ddis_phase_function() override {
        using FloatStorage = DynamicBuffer<Float>;

        // DDIS rebuild is driven exclusively by Scene::parameters_changed()
        // via update_ddis_phase_function(), which ensures all media sharing a
        // phase function via ref are updated before dirty flags are cleared.
        if (m_ddis_threshold <= 0.f || m_ddis_phase_function == nullptr)
            return;

        // Compute envelope nodes from all phase function
        std::vector<FloatStorage> nodes_list;
        for (ref<Base> component: m_components) {
            const PhaseFunction *phase = component->phase_function();
            nodes_list.push_back( phase->get_envelope_nodes() );
        }
        FloatStorage nodes = merge_envelope_nodes( nodes_list );

        // Compute envelope values from envelope nodes and phase function
        FloatStorage values = dr::zeros<FloatStorage>(dr::width(nodes));
        for (ref<Base> component: m_components) {
            const PhaseFunction *phase = component->phase_function();
            phase->accumulate_envelope(nodes, values);
        }

        struct ValuesCallback : TraversalCallback {
            FloatStorage *target_nodes = nullptr;
            FloatStorage *target_values = nullptr;
            void put_value(std::string_view name, void *ptr, uint32_t,
                            const std::type_info &) override {
                if (name == "nodes")
                    target_nodes = static_cast<FloatStorage *>(ptr);
                if (name == "values")
                    target_values = static_cast<FloatStorage *>(ptr);
            }
            void put_object(std::string_view, Object *, uint32_t) override {}
        } cb;

        m_ddis_phase_function->traverse(&cb);

        Assert(cb.target_values && cb.target_nodes);

        if (cb.target_values)
            *cb.target_values = values;

        if (cb.target_nodes)
            *cb.target_nodes = nodes;

        if (cb.target_values || cb.target_nodes)
            m_ddis_phase_function->parameters_changed({});
    }

protected:
    ref<PhaseFunction> create_ddis_phase_function() override {
        using FloatStorage = DynamicBuffer<Float>;

        // Compute envelope nodes from all phase function
        std::vector<FloatStorage> nodes_list;
        for (ref<Base> component: m_components) {
            const PhaseFunction *phase = component->phase_function();
            nodes_list.push_back( phase->get_envelope_nodes() );
        }
        FloatStorage nodes = merge_envelope_nodes( nodes_list );

        // Compute envelope values from envelope nodes and phase function
        FloatStorage values = dr::zeros<FloatStorage>(dr::width(nodes));
        for (ref<Base> component: m_components) {
            const PhaseFunction *phase = component->phase_function();
            phase->accumulate_envelope(nodes, values);
        }

        auto pmgr = PluginManager::instance();
        Properties props_ddis("tabphase_irregular");
        size_t shape = nodes.size();
        props_ddis.set_any("nodes", TensorXf(std::move(nodes), 1, &shape));
        props_ddis.set_any("values", TensorXf(std::move(values), 1, &shape));
        return pmgr->create_object<PhaseFunction>(props_ddis);
    }

private:
    std::vector<ref<Base>> m_components;
    DynamicBuffer<PhaseFunctionPtr> m_phase_functions_dr;
    ScalarBoundingBox3f m_aabb;

    MI_TRAVERSE_CB(Base, m_aabb, m_phase_functions_dr)
};

MI_EXPORT_PLUGIN(OverlappingMedium)
NAMESPACE_END(mitsuba)
