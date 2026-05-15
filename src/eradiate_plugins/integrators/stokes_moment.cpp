#include <mitsuba/render/integrator.h>
#include <mitsuba/render/records.h>
#include <mitsuba/render/mueller.h>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _integrator-stokes_moment:

Stokes-moment integrator (:monosp:`stokes_moment`)
-----------------------------------------------------

.. pluginparameters::

 * - use_stokes
   - |bool|
   - Enable Stokes vector AOVs. Requires a polarized rendering mode. Default: false.

 * - use_moment
   - |bool|
   - Enable second-moment AOVs for variance estimation. Default: false.

 * - meridian_align
   - |bool|
   - When :paramtype:`use_stokes` is true, align the Stokes vector to the meridian
     plane instead of the sensor's x-axis. Default: false.

 * - (Nested plugin)
   - :paramtype:`integrator`
   - Sub-integrator (exactly one must be specified) which will be sampled along
     the integrator.

This integrator combines the :ref:`stokes <integrator-stokes>` and
:ref:`moment <integrator-moment>` integrators into a single plugin, resolving
the issue that :monosp:`stokes` cannot be nested inside :monosp:`moment`
(the sensor pointer would not be set correctly). Enable either or both modes
via :paramtype:`use_stokes` and :paramtype:`use_moment`.

AOV layout mirrors the behaviour of the integrators it combines. When both modes 
are active, the layout is (N = child AOV count):

.. code-block:: text

    [0..11]                  name.S0.R name.S0.G name.S0.B  name.S1.R ... name.S3.B
    [12..12+N-1]             child AOVs
    [12+N..12+N+2]           name.X  name.Y  name.Z
    [12+N+3..2*(12+N+3)-1]   m2_ mirror of all preceding channels

This plugins departs from the :ref:`stokes <integrator-stokes>` integrator by
prepending the nested integrator name to the stokes AOV names.
 */

template <typename Float, typename Spectrum>
class StokesMomentIntegrator final : public SamplingIntegrator<Float, Spectrum> {
public:
    MI_IMPORT_BASE(SamplingIntegrator)
    MI_IMPORT_TYPES(Scene, Sensor, Sampler, Medium)

    StokesMomentIntegrator(const Properties &props) : Base(props) {
        m_use_stokes   = props.get<bool>("use_stokes", false);
        m_use_moment   = props.get<bool>("use_moment", false);
        m_meridian_align = props.get<bool>("meridian_align", false);

        if (!m_use_stokes && !m_use_moment)
            Throw("At least one of 'use_stokes' or 'use_moment' must be true.");

        if (m_use_stokes && !is_polarized_v<Spectrum>)
            Throw("'use_stokes' requires a polarized rendering mode.");

        for (auto &prop : props.objects()) {
            Base *integrator = prop.try_get<Base>();
            if (!integrator)
                Throw("Child objects must be of type 'SamplingIntegrator'.");
            if (m_integrator)
                Throw("Only one nested integrator may be specified.");
            m_integrator      = integrator;
            m_integrator_name = std::string(prop.name());
        }

        if (!m_integrator)
            Throw("A nested integrator must be specified.");

        // Build AOV name list
        std::vector<std::string> child_aovs = m_integrator->aov_names();
        m_child_aov_count = child_aovs.size();

        if (m_use_stokes) {
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 3; ++j)
                    m_aov_names.push_back(
                        m_integrator_name + "." + "S" + std::to_string(i) + "." + std::string(1, "RGB"[j])
                    );
        }

        for (auto &name : child_aovs) {
            // Prefix child AOVs when moment is active, matching standalone moment behaviour
            if (m_use_moment)
                m_aov_names.push_back(m_integrator_name + "." + name);
            else
                m_aov_names.push_back(name);
        }

        if (m_use_moment) {
            m_aov_names.push_back(m_integrator_name + ".X");
            m_aov_names.push_back(m_integrator_name + ".Y");
            m_aov_names.push_back(m_integrator_name + ".Z");
            // Append m2_ mirror of every first-moment AOV
            size_t first_moment_count = m_aov_names.size();
            for (size_t i = 0; i < first_moment_count; ++i)
                m_aov_names.push_back("m2_" + m_aov_names[i]);
        }
    }

    std::pair<Spectrum, Mask> sample(const Scene *scene,
                                     Sampler *sampler,
                                     const RayDifferential3f &ray,
                                     const Medium *medium,
                                     Float *aovs,
                                     Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::SamplingIntegratorSample, active);

        size_t stokes_offset = m_use_stokes ? 12 : 0;

        // Child writes its AOVs after the Stokes block (if any)
        auto [spec, mask] = m_integrator->sample(
            scene, sampler, ray, medium, aovs + stokes_offset, active);

        // spec_for_moment may be updated to the rotated spectrum below
        Spectrum spec_for_moment = spec;

        // --- Stokes processing ---
        if (m_use_stokes) {
            if constexpr (is_polarized_v<Spectrum>) {
                if (!m_sensor)
                    Throw("The `sample()` method for this integrator must "
                          "exclusively be called through the `render()` method!");

                Vector3f current_basis = mueller::stokes_basis(-ray.d);
                Vector3f target_basis  = dr::zeros<Vector3f>();

                if (m_meridian_align) {
                    Vector3f vertical        = Vector3f(0.f, 0.f, 1.f);
                    Vector3f tmp             = dr::cross(vertical, -ray.d);
                    Mask ray_is_vertical     = dr::norm(tmp) < math::RayEpsilon<Float>;
                    target_basis[ray_is_vertical]  = Vector3f(1.f, 0.f, 0.f);
                    target_basis[!ray_is_vertical] = dr::cross(dr::normalize(tmp), -ray.d);
                } else {
                    Vector3f vertical = m_sensor->world_transform() * Vector3f(0.f, 1.f, 0.f);
                    target_basis      = dr::cross(ray.d, vertical);
                }

                auto spec_rot = mueller::rotate_stokes_basis(-ray.d, current_basis, target_basis) * spec;
                spec_for_moment = spec_rot;

                Float *aovs_s = aovs;
                for (int i = 0; i < 4; ++i) {
                    Color3f rgb;
                    if constexpr (is_monochromatic_v<Spectrum>) {
                        rgb = spec_rot.entry(i, 0).x();
                    } else if constexpr (is_rgb_v<Spectrum>) {
                        rgb = spec_rot.entry(i, 0);
                    } else {
                        static_assert(is_spectral_v<Spectrum>);
                        auto pdf = pdf_rgb_spectrum(ray.wavelengths);
                        UnpolarizedSpectrum _spec =
                            spec_rot.entry(i, 0) * dr::select(pdf != 0.f, dr::rcp(pdf), 0.f);
                        rgb = spectrum_to_srgb(_spec, ray.wavelengths, active);
                    }
                    *aovs_s++ = rgb.r(); *aovs_s++ = rgb.g(); *aovs_s++ = rgb.b();
                }
            }
        }

        // --- Moment processing ---
        if (m_use_moment) {
            size_t first_moment_count = stokes_offset + m_child_aov_count + 3;
            Float *aovs_xyz = aovs + stokes_offset + m_child_aov_count;

            UnpolarizedSpectrum spec_u = unpolarized_spectrum(spec_for_moment);

            Color3f xyz;
            if constexpr (is_monochromatic_v<Spectrum>) {
                xyz = spec_u.x();
            } else if constexpr (is_rgb_v<Spectrum>) {
                xyz = srgb_to_xyz(spec_u, active);
            } else {
                static_assert(is_spectral_v<Spectrum>);
                auto pdf = pdf_rgb_spectrum(ray.wavelengths);
                spec_u *= dr::select(pdf != 0.f, dr::rcp(pdf), 0.f);
                xyz = spectrum_to_xyz(spec_u, ray.wavelengths, active);
            }

            *aovs_xyz++ = xyz.x(); *aovs_xyz++ = xyz.y(); *aovs_xyz++ = xyz.z();
            // aovs_xyz now points to the start of the m2_ block

            // Square every first-moment value into its mirror m2_ slot
            Float *m2_start = aovs_xyz;
            for (size_t k = 0; k < first_moment_count; ++k)
                m2_start[k] = dr::square(aovs[k]);
        }

        return { spec, mask };
    }

    TensorXf render(Scene *scene,
                    Sensor *sensor,
                    UInt32 seed    = 0,
                    uint32_t spp   = 0,
                    bool develop   = true,
                    bool evaluate  = true) override {
        m_sensor = sensor;
        TensorXf result = Base::render(scene, sensor, seed, spp, develop, evaluate);
        m_sensor = nullptr;
        return result;
    }

    std::vector<std::string> aov_names() const override {
        return m_aov_names;
    }

    void traverse(TraversalCallback *cb) override {
        cb->put("integrator", m_integrator, ParamFlags::Differentiable);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "StokesMomentIntegrator[" << std::endl
            << "  use_stokes = " << m_use_stokes << "," << std::endl
            << "  use_moment = " << m_use_moment << "," << std::endl
            << "  integrator = " << string::indent(m_integrator, 2) << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS(StokesMomentIntegrator)

private:
    bool m_use_stokes;
    bool m_use_moment;
    bool m_meridian_align;
    Sensor *m_sensor = nullptr;
    ref<Base> m_integrator;
    std::string m_integrator_name;
    std::vector<std::string> m_aov_names;
    size_t m_child_aov_count;

    MI_TRAVERSE_CB(Base, m_integrator)
};

MI_EXPORT_PLUGIN(StokesMomentIntegrator)
NAMESPACE_END(mitsuba)
