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
   - The inner medium to repeat.

 * - lattice
   - |vector|
   - Axis-aligned lattice period along x, y, z. Default: the extents of the
     inner medium's domain.

 * - aabb_min, aabb_max
   - |point|
   - Bounds of the tiled region in world space. Required.

This plugin periodically tiles an inner medium (which may itself be an
aggregate, e.g. :ref:`overlap <medium-overlap>`) over a
translation lattice:

.. math::

    \sigma_t(x) = \begin{cases}
        \sigma_t^\mathrm{inner}(\mathrm{fold}(x)) & x \in \mathrm{domain},\ \mathrm{fold}(x) \in \mathrm{domain}^\mathrm{inner} \\
        0 & \text{otherwise}
    \end{cases}

Point queries are folded into the lattice period by this plugin. The lattice,
can be larger than the inner medium's bounding box, creating a gap in which
point queries read zero. Inner's bounding box must fit inside one lattice
period, and the inner medium must not itself contain a :monosp:`repeat`.
The inner medium and integrator are unaware of the tiling.
*/
template <typename Float, typename Spectrum>
class RepeatMedium final : public Medium<Float, Spectrum> {
public:
    MI_IMPORT_BASE(Medium, m_is_homogeneous, m_has_spectral_extinction,
                    m_phase_function, m_extrema,
                    m_ddis_phase_function, m_ddis_threshold,
                    ddis_phase_function, ddis_threshold
                )
    MI_IMPORT_TYPES(Scene, Sampler, PhaseFunction, PhaseFunctionPtr,
                    Extremum)

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
        m_phase_function          = m_inner->phase_function();

        const auto &inner_extrema = m_inner->extrema();
        if (inner_extrema.empty())
            Throw("repeat: the nested medium has no extremum structure");


        // Get periodic bbox, use the extremum because the medium does not expose bbox.
        ScalarBoundingBox3f pbox;
        for (const auto &structure : inner_extrema) {
            pbox.expand(structure->bbox());
            if (structure->class_name() == "RepeatExtremum")
                Throw("Cannot nest repeat media!");
        }

        if (!pbox.valid() || !dr::all(dr::isfinite(pbox.min) && dr::isfinite(pbox.max)))
            Throw("The nested medium must have a finite domain (one pbox), got %s", pbox);

        m_origin      = pbox.min;
        m_lattice     = props.get<ScalarVector3f>("lattice", pbox.extents());
        m_lattice_rcp = dr::rcp(m_lattice);

        if (dr::any(pbox.extents() > m_lattice))
            Throw("The nested medium bounds %s do not fit inside a lattice period %s",
                  pbox, m_lattice);

        m_aabb = ScalarBoundingBox3f(props.get<ScalarPoint3f>("aabb_min"),
                                     props.get<ScalarPoint3f>("aabb_max"));

        for (const auto &structure : inner_extrema) {
            Properties props_repeat("repeatextremum");
            props_repeat.set("inner", (Object *) structure.get());
            props_repeat.set("lattice", m_lattice);
            props_repeat.set("bbox_min", m_aabb.min);
            props_repeat.set("bbox_max", m_aabb.max);
            auto extremum = PluginManager::instance()->create_object<Extremum>(props_repeat);
            m_extrema.push_back(extremum);
        }

        m_ddis_threshold = m_inner->ddis_threshold();
        m_ddis_phase_function = const_cast<PhaseFunction *>(m_inner->ddis_phase_function());
    }

    void traverse(TraversalCallback *cb) override {
        cb->put("inner", m_inner, ParamFlags::Differentiable);
        if (m_ddis_phase_function != nullptr)
            cb->put("ddis_phase_function",m_ddis_phase_function, ParamFlags::Differentiable);
        Base::traverse(cb);
    }

    UnpolarizedSpectrum
    get_majorant(const MediumInteraction3f &mei, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::MediumEvaluate, active);
        return m_inner->get_majorant(fold(mei), active);
    }

    std::tuple<UnpolarizedSpectrum, UnpolarizedSpectrum, UnpolarizedSpectrum>
    get_scattering_coefficients(const MediumInteraction3f &mei,
                                Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::MediumEvaluate, active);
        MediumInteraction3f folded = fold(mei);
        Mask inside = m_inner->in_aabb(folded.p);
        auto [sigmas, sigman, sigmat] =
            m_inner->get_scattering_coefficients(folded, active);
        return { dr::select(inside, sigmas, 0.f),
                 dr::select(inside, sigman, sigman + sigmat),
                 dr::select(inside, sigmat, 0.f) };
    }

    std::pair<UnpolarizedSpectrum, UnpolarizedSpectrum>
    eval_component(const MediumInteraction3f &mei, const UInt32 &component,
                   Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::MediumEvaluate, active);
        MediumInteraction3f folded = fold(mei);
        Mask inside = m_inner->in_aabb(folded.p);
        auto [sigma_s, sigma_t] =
            m_inner->eval_component(folded, component, active && inside);
        return { dr::select(inside, sigma_s, 0.f),
                 dr::select(inside, sigma_t, 0.f) };
    }

    std::pair<PhaseFunctionPtr, MediumInteraction3f>
    phase_function(const MediumInteraction3f &mei, const UInt32 &component,
                   Mask active) const override {
        return m_inner->phase_function(fold(mei), component, active);
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
    /// Fold the interaction point into the canonical lattice period
    MediumInteraction3f fold(const MediumInteraction3f &mei) const {
        Vector3f c = (mei.p - m_origin) * m_lattice_rcp;
        MediumInteraction3f folded = mei;
        folded.p = mei.p - m_lattice * dr::floor(c);
        return folded;
    }

private:
    ref<Base> m_inner;
    ScalarVector3f m_lattice;
    ScalarVector3f m_lattice_rcp;
    ScalarPoint3f m_origin;
    ScalarBoundingBox3f m_aabb;

    MI_TRAVERSE_CB(Base, m_inner)
};

MI_EXPORT_PLUGIN(RepeatMedium)
NAMESPACE_END(mitsuba)
