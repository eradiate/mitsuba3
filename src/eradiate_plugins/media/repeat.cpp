#include <mitsuba/core/frame.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/interaction.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/eradiate/extremum.h>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _medium-repeat:

Repeating medium (:monosp:`repeat`)
-----------------------------------

.. pluginparameters::

 * - (Nested plugin)
   - |medium|
   - The inner medium to tile. Its domain defines the canonical tile.

 * - lattice
   - |vector|
   - Axis-aligned lattice period along x, y, z. Default: the extents of the
     inner medium's domain.

 * - aabb_min, aabb_max
   - |point|
   - Bounds of the tiled region in world space. Default: infinite.

This plugin periodically tiles an inner medium (which may itself be an
aggregate, e.g. :ref:`multicomponent <medium-multicomponent>`) over a
translation lattice:

.. math::

    \sigma_t(x) = \begin{cases}
        \sigma_t^\mathrm{inner}(\mathrm{fold}(x)) & x \in \mathrm{domain} \\
        0 & \text{otherwise}
    \end{cases}

All folding happens in this plugin; the inner medium and the integrator are
unaware of the tiling. The extremum structure of the inner medium is wrapped
into an :ref:`extremum_repeat <extremum-extremum_repeat>` structure.
*/
template <typename Float, typename Spectrum>
class RepeatMedium final : public Medium<Float, Spectrum> {
public:
    MI_IMPORT_BASE(Medium, m_is_homogeneous, m_has_spectral_extinction,
                    m_phase_function, m_extremum_structure,
                    m_ddis_phase_function, m_ddis_threshold,
                    ddis_phase_function, ddis_threshold
                )
    MI_IMPORT_TYPES(Scene, Sampler, ExtremumStructure, PhaseFunction,
                    PhaseFunctionPtr)

    using MediumContext = MediumContext<Float, Spectrum>;

    RepeatMedium(const Properties &props) : Base(props) {
        m_is_homogeneous = false;

        for (auto &prop : props.objects()) {
            if (auto *inner = prop.try_get<Base>()) {
                if (m_inner)
                    Throw("repeat accepts a single nested medium");
                m_inner = inner;
            }
        }
        if (!m_inner)
            Throw("repeat requires a nested medium");

        m_has_spectral_extinction = m_inner->has_spectral_extinction();
        m_phase_function = m_inner->phase_function();

        ExtremumStructure *inner_structure = m_inner->extremum_structure();
        if (!inner_structure)
            Throw("repeat: the nested medium has no extremum structure");

        ScalarBoundingBox3f tile = inner_structure->bbox();
        if (!tile.valid() ||
            !dr::all(dr::isfinite(tile.min) && dr::isfinite(tile.max)))
            Throw("repeat: the nested medium must have a finite domain (one "
                  "tile), got %s", tile);
        m_origin = tile.min;

        // The lattice frame is stated once, here, and shared with the
        // extremum structure
        Properties props_repeat("extremum_repeat");
        props_repeat.set("structure", (Object *) inner_structure);

        m_lattice = props.get<ScalarVector3f>("lattice", tile.extents());
        m_lattice_rcp = 1.f / m_lattice;
        props_repeat.set("lattice", m_lattice);

        if (props.has_property("aabb_min") && props.has_property("aabb_max")) {
            m_aabb = ScalarBoundingBox3f(props.get<ScalarPoint3f>("aabb_min"),
                                         props.get<ScalarPoint3f>("aabb_max"));
            props_repeat.set("aabb_min", ScalarPoint3f(m_aabb.min));
            props_repeat.set("aabb_max", ScalarPoint3f(m_aabb.max));
        } else {
            m_aabb = ScalarBoundingBox3f(
                ScalarPoint3f(-dr::Infinity<ScalarFloat>),
                ScalarPoint3f(dr::Infinity<ScalarFloat>));
        }

        m_extremum_structure =
            PluginManager::instance()->create_object<ExtremumStructure>(
                props_repeat);

        m_ddis_threshold = m_inner->ddis_threshold();
        m_ddis_phase_function =
            const_cast<PhaseFunction *>(m_inner->ddis_phase_function());
    }

    void traverse(TraversalCallback *cb) override {
        cb->put("inner", m_inner, ParamFlags::Differentiable);
        if (m_ddis_phase_function != nullptr)
            cb->put("ddis_phase_function", m_ddis_phase_function,
                    ParamFlags::Differentiable);
        Base::traverse(cb);
    }

    UnpolarizedSpectrum
    get_majorant(const MediumInteraction3f &mi, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::MediumEvaluate, active);
        return m_inner->get_majorant(fold(mi), active);
    }

    UnpolarizedSpectrum
    get_minorant(const MediumInteraction3f &mi, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::MediumEvaluate, active);
        return m_inner->get_minorant(fold(mi), active);
    }

    std::tuple<UnpolarizedSpectrum, UnpolarizedSpectrum, UnpolarizedSpectrum>
    get_scattering_coefficients(const MediumInteraction3f &mi,
                                Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::MediumEvaluate, active);
        return m_inner->get_scattering_coefficients(fold(mi), active);
    }

    MediumContext get_medium_context(const MediumInteraction3f &mei,
        UnpolarizedSpectrum sigma_maj, Float sample, Mask active) const override {
        return m_inner->get_medium_context(fold(mei), sigma_maj, sample, active);
    }

    PhaseFunctionPtr phase_function(const MediumInteraction3f &mi,
                                    Float sample,
                                    Mask active) const override {
        return m_inner->phase_function(fold(mi), sample, active);
    }

    PhaseFunctionPtr phase_function(const UInt32 &component,
                                    Mask active /*= true*/) const override {
        return m_inner->phase_function(component, active);
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
        oss << "RepeatMedium[" << std::endl
            << "  inner = " << string::indent(m_inner) << "," << std::endl
            << "  lattice = " << m_lattice << "," << std::endl
            << "  aabb = " << m_aabb << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS(RepeatMedium)

private:
    /// Fold the interaction point into the canonical tile
    MediumInteraction3f fold(const MediumInteraction3f &mi) const {
        Vector3f c = (mi.p - m_origin) * m_lattice_rcp;
        MediumInteraction3f folded = mi;
        folded.p = mi.p - m_lattice * dr::floor(c);
        return folded;
    }

private:
    ref<Base> m_inner;
    ScalarVector3f m_lattice = ScalarVector3f(1.f);
    ScalarVector3f m_lattice_rcp = ScalarVector3f(1.f);
    ScalarPoint3f m_origin = 0.f;
    ScalarBoundingBox3f m_aabb;

    MI_TRAVERSE_CB(Base, m_inner)
};

MI_EXPORT_PLUGIN(RepeatMedium)
NAMESPACE_END(mitsuba)
