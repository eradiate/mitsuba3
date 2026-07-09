#include <mitsuba/core/frame.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/interaction.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/eradiate/extremum.h>
#include <mitsuba/render/eradiate/phase_utils.h>

NAMESPACE_BEGIN(mitsuba)


/**!

.. _medium-multicomponent:

Multi Component medium (:monosp:`multicomponent`)
--------------------------------------------------------------------

.. pluginparameters::

 * - periodic_min, periodic_max
   - |point|
   - Periodic medium bounding box.

 * - aabb_min, aabb_max
   - |point|
   - Medium bounding box

 * - periodic_condition
   - |string|
   - Periodic condition type. Can be {null, repeat}.

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


This plugin provides an aggregate medium implementation. It accepts multiple
nested media components, making the assumption that those media do not have a
physical boundary, and that the data is not defined outside of their bounding
box.
*/
template <typename Float, typename Spectrum>
class MultiComponentMedium final : public Medium<Float, Spectrum> {
public:
    MI_IMPORT_BASE(Medium, m_is_homogeneous, m_has_spectral_extinction,
                    m_ddis_phase_function, m_ddis_threshold,
                    create_ddis_phase_function
                )
    MI_IMPORT_TYPES(Scene, Sampler, MediumPtr, ExtremumStructure,
                    ExtremumStructurePtr, PhaseFunction)
    using FloatStorage = DynamicBuffer<Float>;

    MultiComponentMedium(const Properties &props) : Base(props) {
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

        m_components_dr = dr::load<DynamicBuffer<MediumPtr>>(
            m_components.data(), m_components.size());
        dr::eval(m_components_dr);

        // TODO: Extremum Structure compatibility, lets first make it work
        //       with a more standard volpath.

        m_ddis_threshold = props.get<ScalarFloat>("ddis_threshold", 0.1f);

        if (m_ddis_threshold > 0.f) {
            m_ddis_phase_function = static_cast<PhaseFunction*>(create_ddis_phase_function());
        }

        // Periodic bounding box
        // TODO:
        //  - determine proper default values for periodic and medium bbox.
        //  - determine the behaviour for homogeneous media.
        if (props.has_property("paabb_min") && props.has_property("paabb_max")) {
            ScalarPoint3f paabb_min = props.get<ScalarPoint3f>("paabb_min");
            ScalarPoint3f paabb_max = props.get<ScalarPoint3f>("paabb_max");
            m_paabb = ScalarBoundingBox3f(paabb_min, paabb_max);
        }

        // Medium bounding box
        if (props.has_property("aabb_min") && props.has_property("aabb_max")) {
            ScalarPoint3f aabb_min = props.get<ScalarPoint3f>("aabb_min");
            ScalarPoint3f aabb_max = props.get<ScalarPoint3f>("aabb_max");
            m_aabb = ScalarBoundingBox3f(aabb_min, aabb_max);
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

    void parameters_changed(const std::vector<std::string> &/*keys*/) override {
        // for now empty..
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

    UnpolarizedSpectrum
    get_minorant(const MediumInteraction3f &mi, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::MediumEvaluate, active);
        UnpolarizedSpectrum minorant = m_components[0]->get_minorant(mi, active);
        for (size_t i = 1; i < m_components.size(); ++i) {
            minorant += m_components[i]->get_minorant(mi, active);
        }
        return minorant;
    }

    std::tuple<UnpolarizedSpectrum, UnpolarizedSpectrum, UnpolarizedSpectrum>
    get_scattering_coefficients(const MediumInteraction3f &mi,
                                Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::MediumEvaluate, active);

        UnpolarizedSpectrum sigmas(0.f), sigmat(0.f);

        for (size_t i = 0; i < m_components.size(); ++i) {
            Mask accumulate = active && m_components[i]->in_aabb(mi.p);
            if (dr::any_or<true>(accumulate)) {
                auto [c_sigmas, c_sigman, c_sigmat] =
                    m_components[i]->get_scattering_coefficients(mi, accumulate);
                sigmas += c_sigmas;
                sigmat += c_sigmat;
            }
        }

        // Each component's own sigma_n is relative to *its* local majorant,
        // not the composite's combined one, so it cannot be summed directly.
        // Re-derive sigma_n from the combined majorant instead.
        UnpolarizedSpectrum sigman = get_majorant(mi, active) - sigmat;

        return { sigmas, sigman, sigmat };
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
        oss << "MultiComponentMedium[" << std::endl
            << "  ddis_phase_function   = " << string::indent(m_ddis_phase_function) << "," << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS(MultiComponentMedium)


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

        Assert(cb.taget_values && cb.target_nodes);

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
    DynamicBuffer<MediumPtr> m_components_dr;
    ScalarBoundingBox3f m_aabb;
    ScalarBoundingBox3f m_paabb;

    MI_TRAVERSE_CB(Base, m_paabb, m_aabb)
};

MI_EXPORT_PLUGIN(MultiComponentMedium)
NAMESPACE_END(mitsuba)
