#include <mitsuba/core/properties.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/interaction.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/render/sampler.h>
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
   - Bounds of the tiled region in world space. Required.

This plugin periodically tiles an inner medium (which may itself be an
aggregate, e.g. :ref:`multicomponent <medium-multicomponent>`) over a
translation lattice:

.. math::

    \sigma_t(x) = \begin{cases}
        \sigma_t^\mathrm{inner}(\mathrm{fold}(x)) & x \in \mathrm{domain} \\
        0 & \text{otherwise}
    \end{cases}

Point queries are folded into the canonical tile by this plugin. The inner
medium's extremum structures are traversed tile by tile over the lattice, so
the inner medium and the integrator are unaware of the tiling. The canonical
tile must fit inside one lattice period, and the inner medium must not itself
contain a :monosp:`repeat`.
*/
template <typename Float, typename Spectrum>
class RepeatMedium final : public Medium<Float, Spectrum> {
public:
    MI_IMPORT_BASE(Medium, m_is_homogeneous, m_has_spectral_extinction,
                    m_phase_function, m_extrema,
                    m_ddis_phase_function, m_ddis_threshold,
                    ddis_phase_function, ddis_threshold
                )
    MI_IMPORT_TYPES(Scene, Sampler, PhaseFunction, PhaseFunctionPtr)

    using MediumSample = MediumSample<Float, Spectrum>;

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

        m_extrema = m_inner->dda_entries();
        if (m_extrema.empty())
            Throw("repeat: the nested medium has no extremum structure");

        ScalarBoundingBox3f tile;
        for (const auto &entry : m_extrema) {
            if (entry.tiling)
                Throw("repeat: the nested medium is already tiled");
            tile.expand(entry.structure->bbox());
        }
        if (!tile.valid() ||
            !dr::all(dr::isfinite(tile.min) && dr::isfinite(tile.max)))
            Throw("repeat: the nested medium must have a finite domain (one "
                  "tile), got %s", tile);
        m_origin = tile.min;

        m_lattice = props.get<ScalarVector3f>("lattice", tile.extents());
        m_lattice_rcp = 1.f / m_lattice;
        if (dr::any(tile.extents() > m_lattice))
            Throw("repeat: the tile %s does not fit inside a lattice period %s",
                  tile, m_lattice);

        m_aabb = ScalarBoundingBox3f(props.get<ScalarPoint3f>("aabb_min"),
                                     props.get<ScalarPoint3f>("aabb_max"));

        for (auto &entry : m_extrema) {
            ScalarPoint3f cell_min = entry.structure->bbox().min;
            entry.tiling = { ScalarBoundingBox3f(cell_min, cell_min + m_lattice),
                             m_aabb };
        }

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

    MediumSample sample_scattering_properties(const MediumInteraction3f &mei,
        UnpolarizedSpectrum sigma_maj, Float sample, Mask active) const override {
        return m_inner->sample_scattering_properties(fold(mei), sigma_maj, sample, active);
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
