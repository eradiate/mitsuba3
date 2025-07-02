#include <mitsuba/core/distr_1d.h>
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

NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
class PiecewiseSphericalMedium final : public Medium<Float, Spectrum> {
public:
    MI_IMPORT_BASE(Medium, m_is_homogeneous, m_has_spectral_extinction,
                   m_phase_function)
    MI_IMPORT_TYPES(Scene, Sampler, Texture, Volume)

    using ScalarIndex       = uint32_t;
    using ScalarSize        = uint32_t;
    using FloatStorage      = DynamicBuffer<Float>;
    using index             = dr::uint32_array_t<Float>;

    PiecewiseSphericalMedium(const Properties &props) : Base(props) {
        m_is_homogeneous = false;
        m_albedo         = props.volume<Volume>("albedo", 0.75f);
        m_sigmat         = props.volume<Volume>("sigma_t", 1.f);

        m_has_spectral_extinction =
            props.get<bool>("has_spectral_extinction", true);

        m_max_density = dr::opaque<Float>(m_sigmat->max());

        precompute_optical_thickness();

        dr::set_attr(this, "is_homogeneous", m_is_homogeneous);
        dr::set_attr(this, "has_spectral_extinction", m_has_spectral_extinction);
    }

    std::tuple<typename Medium<Float, Spectrum>::MediumInteraction3f, Float, Float>
    sample_interaction_real(const Ray3f &ray, const SurfaceInteraction3f &si,
                            Float sample, UInt32 channel,
                            Mask active) const override {

    }

    std::tuple<Float, Float, Mask>
    eval_transmittance_pdf_real(const Ray3f &ray,
                                const SurfaceInteraction3f &si, UInt32 channel,
                                Mask active) const override {
    }

    void traverse(TraversalCallback *callback) override {
    }

    void parameters_changed(const std::vector<std::string> & /*keys*/ = {}) override {
    }

    void precompute() const override {
        std::vector<ScalarVector3f> directions = precompute_directions();
        //precompute_optical_thickness()
    }

    std::vector<ScalarVector3f> precompute_directions() const {
        ScalarPoint3f p = m_sigmat->bbox().min;

        //  Compute angles based on a regularly spaced interval
        int32_t sided_samples = 4;
        int32_t extent = sided_samples * 2 + 1;
        int32_t total_samples = extent * extent;

        //  Storage for the direction vectors
        //  TODO -> Change this to a drjit::DynamicArray?
        std::vector<ScalarVector3f> directions(total_samples);

        for (int32_t off_x = -sided_samples; off_x < sided_samples + 1; off_x++) {
            for (int32_t off_y = -sided_samples; off_y < sided_samples + 1; off_y++) {    
                ScalarPoint3f offset = p + ScalarPoint3f(off_x, off_y, 1.0f);
                ScalarVector3f direction = drjit::normalize(offset - p);

                int32_t idx = (off_x + sided_samples) * extent + (off_y + sided_samples);
                directions[idx] = direction;
            }
        }

        return directions;
    }

    void precompute_optical_thickness() const {
        static constexpr size_t SpectralSize = dr::array_size_v<UnpolarizedSpectrum>;
        using ScalarUnpolarized = dr::Array<dr::scalar_t<UnpolarizedSpectrum>, SpectralSize>;

        // Check that the first two dimensions are equal to 1 and that z is one
        // or more.
        const ScalarVector3i resolution = m_sigmat->resolution();
        const ScalarVector3f voxel_size = m_sigmat->voxel_size();
        if (resolution.x() > 1 || resolution.y() > 1)
            Throw("PiecewiseMedium: x or y resolution bigger than one, assumed "
                  "shape is [1,1,n]");

        MediumInteraction3f mei = dr::zeros<MediumInteraction3f>();
        ScalarPoint3f min       = m_sigmat->bbox().min;
        ScalarVector3f step     = ScalarVector3d(0.f, 0.f, voxel_size.z());

        //  The interaction point from where we calculate directions
        ScalarPoint3f p         = m_sigmat-> bbox().min;

        //  Precompute

        //  Here we assume we start at the top shell
        for (int32_t i = 0; i < (int32_t) resolution.z(); ++i) {

        }
    }

    UnpolarizedSpectrum get_majorant(const MediumInteraction3f &mi,
                                     Mask active) const override {
    }

    std::tuple<UnpolarizedSpectrum, UnpolarizedSpectrum, UnpolarizedSpectrum>
    get_scattering_coefficients(const MediumInteraction3f &mi,
                                Mask active) const override {
    }

    std::tuple<Mask, Float, Float>
    intersect_aabb(const Ray3f &ray) const override {
        return m_sigmat->bbox().ray_intersect(ray);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "PiecewiseSphericalMedium[" << std::endl
            << "  albedo        = " << string::indent(m_albedo) << std::endl
            << "  sigma_t       = " << string::indent(m_sigmat) << std::endl
            << "]";
        return oss.str();
    }

    static Float extract_channel(UnpolarizedSpectrum value, UInt32 channel) {
        Float result = value[0];
        if constexpr (is_rgb_v<Spectrum>) { // Handle RGB rendering
            dr::masked(result, dr::eq(channel, 1u)) = value[1];
            dr::masked(result, dr::eq(channel, 2u)) = value[2];
        } else {
            DRJIT_MARK_USED(channel);
        }

        return result;
    }

    MI_DECLARE_CLASS()

private:
    ref<Volume> m_sigmat, m_albedo;

    Float m_max_density;
    FloatStorage m_cum_opt_thickness;
    FloatStorage m_reverse_cum_opt_thickness;
};

MI_IMPLEMENT_CLASS_VARIANT(PiecewiseSphericalMedium, Medium)
MI_EXPORT_PLUGIN(PiecewiseSphericalMedium, "Piecewise Spherical Medium")
NAMESPACE_END(mitsuba)