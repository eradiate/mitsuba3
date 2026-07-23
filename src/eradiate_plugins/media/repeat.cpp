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

 * - lattice_0, lattice_1, lattice_2
   - |vector|
   - Lattice translation vectors. May be rotated (tilted tiling); per-tile
     rotation of the content is not supported. Default: the axis-aligned
     extents of the inner medium's domain.

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
                    create_ddis_phase_function
                )
    MI_IMPORT_TYPES(Scene, Sampler, ExtremumStructure, PhaseFunction,
                    PhaseFunctionPtr)

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

        ScalarVector3f l[3];
        for (size_t i = 0; i < 3; ++i) {
            ScalarVector3f def(0.f);
            def[i] = tile.extents()[i];
            l[i] = props.get<ScalarVector3f>("lattice_" + std::to_string(i),
                                             def);
            props_repeat.set("lattice_" + std::to_string(i), l[i]);
        }
        m_lattice = ScalarMatrix3f(l[0].x(), l[1].x(), l[2].x(),
                                   l[0].y(), l[1].y(), l[2].y(),
                                   l[0].z(), l[1].z(), l[2].z());
        m_lattice_inv = dr::inverse(m_lattice);

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

        m_ddis_threshold = props.get<ScalarFloat>("ddis_threshold", 0.1f);
        if (m_ddis_threshold > 0.f)
            m_ddis_phase_function =
                static_cast<PhaseFunction *>(create_ddis_phase_function());
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

    PhaseFunctionPtr phase_function(const MediumInteraction3f &mi,
                                    Float sample,
                                    Mask active) const override {
        return m_inner->phase_function(fold(mi), sample, active);
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
        Vector3f c = Matrix3f(m_lattice_inv) * (mi.p - m_origin);
        MediumInteraction3f folded = mi;
        folded.p = mi.p - Matrix3f(m_lattice) * dr::floor(c);
        return folded;
    }

private:
    ref<Base> m_inner;
    ScalarMatrix3f m_lattice = dr::identity<ScalarMatrix3f>();
    ScalarMatrix3f m_lattice_inv = dr::identity<ScalarMatrix3f>();
    ScalarPoint3f m_origin = 0.f;
    ScalarBoundingBox3f m_aabb;

    MI_TRAVERSE_CB(Base, m_inner)
};

MI_EXPORT_PLUGIN(RepeatMedium)
NAMESPACE_END(mitsuba)
