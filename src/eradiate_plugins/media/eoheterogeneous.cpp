#include <mitsuba/core/frame.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/interaction.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/eradiate/extremum.h>

NAMESPACE_BEGIN(mitsuba)


/**!

.. _medium-eoheterogeneous:

Heterogeneous medium (Earth Observation) (:monosp:`eoheterogeneous`)
--------------------------------------------------------------------

.. pluginparameters::

 * - albedo
   - |float|, |spectrum| or |volume|
   - Single-scattering albedo of the medium (Default: 0.75).
   - |exposed|, |differentiable|

 * - sigma_t
   - |float|, |spectrum| or |volume|
   - Extinction coefficient in inverse scene units (Default: 1).
   - |exposed|, |differentiable|

 * - scale
   - |float|
   - Optional scale factor that will be applied to the extinction parameter.
     It is provided for convenience when accommodating data based on different
     units, or to simply tweak the density of the medium. (Default: 1)
   - |exposed|

 * - sample_emitters
   - |bool|
   - Flag to specify whether shadow rays should be cast from inside the volume (Default: |true|)
     If the medium is enclosed in a :ref:`dielectric <bsdf-dielectric>` boundary,
     shadow rays are ineffective and turning them off will significantly reduce
     render time. This can reduce render time up to 50% when rendering objects
     with subsurface scattering.

 * - (Nested plugin)
   - |phase|
   - A nested phase function that describes the directional scattering properties of
     the medium. When none is specified, the renderer will automatically use an instance of
     isotropic.
   - |exposed|, |differentiable|

 * - ddis_threshold
   - |float|
   - Specifies the probability to importance sample the phase using the emitter as
     incident direction. Set to a negative value to disable. (Default: 0.1)

 * - aabb_min, aabb_max
   - |point3f|
   - Optional override to the medium bounding box. Uses the bounding box of the
     `sigma_t` volume by default.

This plugin provides a flexible heterogeneous medium implementation, which acquires its data
from nested volume instances. These can be constant, use a procedural function, or fetch data from
disk, e.g. using a 3D grid.

This plugin is identical to the `heterogeneous` plugin but accepts additional parameters to set the
medium's bounding box. It is possbile to set volumes that span a portion of the medium's bbox,
and thus exploit wrapping mechanism, e.g. periodic boundaries, outside this portion.

The medium is parametrized by the single scattering albedo and the extinction coefficient
:math:`\sigma_t`. The extinction coefficient should be provided in inverse scene units.
For instance, when a world-space distance of 1 unit corresponds to a meter, the
extinction coefficient should have units of inverse meters. For convenience,
the scale parameter can be used to correct the units. For instance, when the scene is in
meters and the coefficients are in inverse millimeters, set scale to 1000.

Both the albedo and the extinction coefficient can either be constant or textured,
and both parameters are allowed to be spectrally varying.
*/
template <typename Float, typename Spectrum>
class EOHeterogeneousMedium final : public Medium<Float, Spectrum> {
public:
    MI_IMPORT_BASE(Medium, m_is_homogeneous, m_has_spectral_extinction,
                    m_phase_function, m_extremum_structure,
                    m_ddis_phase_function, m_ddis_threshold,
                    create_ddis_phase_function
                )
    MI_IMPORT_TYPES(Scene, Sampler, Texture, Volume, ExtremumStructure,
                    ExtremumStructurePtr, PhaseFunction)
    using FloatStorage = DynamicBuffer<Float>;

    EOHeterogeneousMedium(const Properties &props) : Base(props) {
        m_is_homogeneous = false;
        m_albedo = props.get_volume<Volume>("albedo", 0.75f);
        m_sigmat = props.get_volume<Volume>("sigma_t", 1.0f);

        m_scale = props.get<ScalarFloat>("scale", 1.0f);
        m_has_spectral_extinction = props.get<bool>("has_spectral_extinction", true);

        m_max_density = dr::opaque<Float>(m_scale * m_sigmat->max());
        m_min_density = dr::opaque<Float>(m_scale * m_sigmat->min());

        for (auto &prop : props.objects()) {
            if (auto *extremum = prop.try_get<ExtremumStructure>()) {
                if (m_extremum_structure)
                    Throw("Only a single extremum structure can be specified per medium");
                m_extremum_structure = extremum;
            }
        }

        if (!m_extremum_structure) {
            // Create a default global extremum structure.
            Properties props_extr("extremum_global");
            props_extr.set("volume", (Object *) m_sigmat.get());
            props_extr.set("scale", m_scale);
            m_extremum_structure =
                PluginManager::instance()->create_object<ExtremumStructure>(props_extr);
        }

        m_ddis_threshold = props.get<ScalarFloat>("ddis_threshold", 0.1f);

        if (m_ddis_threshold > 0.f) {
            m_ddis_phase_function = static_cast<PhaseFunction*>(create_ddis_phase_function());
        }

        // Optional user-provided bbox override
        if (props.has_property("aabb_min") && props.has_property("aabb_max")) {
            ScalarPoint3f aabb_min = props.get<ScalarPoint3f>("aabb_min");
            ScalarPoint3f aabb_max = props.get<ScalarPoint3f>("aabb_max");
            m_aabb = ScalarBoundingBox3f(aabb_min, aabb_max);
        }
    }

    void traverse(TraversalCallback *cb) override {
        cb->put("scale",   m_scale,  ParamFlags::NonDifferentiable);
        cb->put("albedo",  m_albedo, ParamFlags::Differentiable);
        cb->put("sigma_t", m_sigmat, ParamFlags::Differentiable);
        if (m_ddis_phase_function != nullptr)
            cb->put("ddis_phase_function", m_ddis_phase_function, ParamFlags::Differentiable);
        Base::traverse(cb);
    }

    void parameters_changed(const std::vector<std::string> &/*keys*/) override {
        m_max_density = dr::opaque<Float>(m_scale * m_sigmat->max());
        m_min_density = dr::opaque<Float>(m_scale * m_sigmat->min());
    }

    UnpolarizedSpectrum
    get_majorant(const MediumInteraction3f & /* mi */,
                 Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::MediumEvaluate, active);
        return m_max_density;
    }

    UnpolarizedSpectrum
    get_minorant(const MediumInteraction3f & /* mi */,
                 Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::MediumEvaluate, active);
        return m_min_density;
    }

    std::tuple<UnpolarizedSpectrum, UnpolarizedSpectrum, UnpolarizedSpectrum>
    get_scattering_coefficients(const MediumInteraction3f &mi,
                                Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::MediumEvaluate, active);

        auto sigmat = m_scale * m_sigmat->eval(mi, active);
        if (has_flag(m_phase_function->flags(), PhaseFunctionFlags::Microflake))
            sigmat *= m_phase_function->projected_area(mi, active);

        auto sigmas = sigmat * m_albedo->eval(mi, active);
        auto sigman = m_max_density - sigmat;
        return { sigmas, sigman, sigmat };
    }

    std::tuple<Mask, Float, Float>
    intersect_aabb(const Ray3f &ray) const override {
        if (m_aabb.valid()) {
            return m_aabb.ray_intersect(ray);
        }
        return m_sigmat->bbox().ray_intersect(ray);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "EOHeterogeneousMedium[" << std::endl
            << "  albedo  = " << string::indent(m_albedo) <<  "," << std::endl
            << "  sigma_t = " << string::indent(m_sigmat) << "," << std::endl
            << "  scale   = " << string::indent(m_scale) << "," << std::endl
            << "  ddis_phase_function   = " << string::indent(m_ddis_phase_function) << "," << std::endl
            << "  extremum = "<< string::indent(m_extremum_structure) << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS(EOHeterogeneousMedium)
private:
    ref<Volume> m_sigmat, m_albedo;
    ScalarFloat m_scale;
    Float m_max_density;
    Float m_min_density;
    ScalarBoundingBox3f m_aabb;

    MI_TRAVERSE_CB(Base, m_sigmat, m_albedo, m_scale, m_max_density, m_aabb)
};

MI_EXPORT_PLUGIN(EOHeterogeneousMedium)
NAMESPACE_END(mitsuba)
