#include <mitsuba/core/properties.h>
#include <mitsuba/render/interaction.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/eradiate/extremum.h>
#include <mitsuba/render/eradiate/phase_utils.h>
#include <mitsuba/render/eradiate/volume_utils.h>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _medium-overlap:

Overlapping medium (:monosp:`overlap`)
--------------------------------------------------------------------

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
of the components' coefficients, and the phase function at a scattering
event is sampled among the components proportionally to their scattering
coefficient.

The components' extremum structures are aggregated into an
:ref:`extremum_overlap <extremum-extremum_overlap>` structure, so no
combined majorant structure needs to be built.
*/
template <typename Float, typename Spectrum>
class OverlappingMedium final : public Medium<Float, Spectrum> {
public:
    MI_IMPORT_BASE(Medium, m_is_homogeneous, m_has_spectral_extinction,
                    m_phase_function, m_extremum_structure,
                    m_ddis_phase_function, m_ddis_threshold,
                    create_ddis_phase_function
                )
    MI_IMPORT_TYPES(Scene, Sampler, MediumPtr, ExtremumStructure,
                    ExtremumStructurePtr, PhaseFunction, PhaseFunctionPtr)

    using FloatStorage = DynamicBuffer<Float>;
    using MediumSample = MediumSample<Float, Spectrum>;

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

        if (m_components.size() > MAX_OVERLAPPING_VOLUMES)
            Throw("overlapping medium: too many components (%zu > %zu)", m_components.size(), MAX_OVERLAPPING_VOLUMES);

        m_components_dr = dr::load<DynamicBuffer<MediumPtr>>(
            m_components.data(), m_components.size());
        dr::eval(m_components_dr);

        m_phase_function = m_components[0]->phase_function();

        // Aggregate the components' built structures; the aggregate is
        // assembled rather than built.
        Properties props_overlap("extremum_overlap");
        for (size_t i = 0; i < m_components.size(); ++i) {
            ExtremumStructure *structure = m_components[i]->extremum_structure();
            if (!structure)
                Throw("overlapping medium: component %s has no extremum structure",
                      m_components[i]->to_string());
            props_overlap.set("structure_" + std::to_string(i),
                              (Object *) structure);
        }
        m_extremum_structure =
            PluginManager::instance()->create_object<ExtremumStructure>(
                props_overlap);

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
            m_aabb = m_extremum_structure->bbox();
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

    PhaseFunctionPtr phase_function(const UInt32 &component, Mask active /*= true*/) const override {
        return dr::gather<MediumPtr>(m_components_dr, component, active)->phase_function();
    }

    MediumSample sample_scattering_properties(const MediumInteraction3f &mei,
        UnpolarizedSpectrum sigma_maj, Float sample, Mask active) const override {
        // Sample a component proportionally to its scattering coefficient
        // at the interaction poin
        MediumSample ms = dr::zeros<MediumSample>();
        std::array<Float, MAX_OVERLAPPING_VOLUMES> cdf;
        UnpolarizedSpectrum summed_sigmas = 0.f;
        UnpolarizedSpectrum summed_sigmat = 0.f;

        // compute the total sigma numbers.
        for (size_t i = 0; i < m_components.size(); ++i) {
            Mask accumulate = active && m_components[i]->in_aabb(mei.p);
            if (dr::any_or<true>(accumulate)) {
                auto [c_sigmas, c_sigman, c_sigmat] =
                    m_components[i]->get_scattering_coefficients(mei, accumulate);
                summed_sigmas += dr::select(accumulate, c_sigmas, 0.f);
                summed_sigmat += dr::select(accumulate, c_sigmat, 0.f);
            }
            cdf[i] = dr::mean(summed_sigmas);
        }

        dr::masked(ms.sigma_s, active) = summed_sigmas;
        dr::masked(ms.sigma_n, active) = sigma_maj - summed_sigmat;
        dr::masked(ms.sigma_t, active) = summed_sigmat;

        // Sample medium component from CDF.
        sample *= cdf[m_components.size()-1];
        UInt32 index = 0;
        for (size_t i = 0; i + 1 < m_components.size(); ++i)
            dr::masked(index, sample >= cdf[i]) = (uint32_t) i + 1;
        dr::masked(ms.sampled_component, active) = index;


        return ms;
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
    DynamicBuffer<MediumPtr> m_components_dr;
    ScalarBoundingBox3f m_aabb;

    MI_TRAVERSE_CB(Base, m_aabb, m_components_dr)
};

MI_EXPORT_PLUGIN(OverlappingMedium)
NAMESPACE_END(mitsuba)
