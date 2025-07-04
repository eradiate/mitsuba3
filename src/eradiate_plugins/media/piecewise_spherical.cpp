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

        precompute();

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
        std::vector<ScalarFloat> angles = precompute_directions();
        precompute_optical_thickness(angles);
    }

    std::vector<ScalarFloat> precompute_directions() const {
        //  Compute angles based on a regularly spaced interval
        int32_t sided_samples = 0;
        int32_t extent = sided_samples * 2 + 1;
        int32_t total_samples = extent;

        //  Storage for the direction vectors
        //  TODO -> Change this to a drjit::DynamicArray?
        std::vector<ScalarFloat> angles(total_samples);
    
        //  Critical angle (define how?)
        ScalarFloat critical_angle = 0.0f;
        ScalarFloat step = 0.0f;

        if (total_samples != 1) {
            critical_angle = 80.0f * (dr::Pi<ScalarFloat> / 180.0f);
            step = critical_angle / sided_samples;
        }

        ScalarFloat to_deg = 180.0f / dr::Pi<ScalarFloat>;

        //  Compute all angles theta for sampling
        for (int32_t idx_x = -sided_samples; idx_x < sided_samples + 1; idx_x++) {
            //  Based on step, compute angle
            ScalarFloat theta = idx_x * step;
            ScalarFloat alpha = theta + (dr::Pi<ScalarFloat> / 2.0f);

            angles[idx_x + sided_samples] = alpha;
        }

        return angles;
    }

    void precompute_optical_thickness(const std::vector<ScalarFloat> angles) const {
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
        ScalarPoint3f p         = m_sigmat-> bbox().max;
        ScalarFloat r_max       = p.z();

        ScalarFloat to_deg = 180.0f / dr::Pi<ScalarFloat>;
    

        //  Here we assume we start at the top shell
        for (int32_t i = 0; i < (int32_t) angles.size(); ++i) {
            //  We assume that the voxel size is the radius of a spherical shell!
            //  The intersection point with the shell from P can then be found by 
            //  using the equation of a circle (x^2 + y^2 = r^2) and of the line
            //  defining the directional through P (y = theta * x + r_max) where 
            //  r_max is defined by the z-coordinate of the max bounding box point.
            ScalarFloat alpha = angles[i];
            ScalarFloat tan_alpha = dr::tan(alpha);

            std::vector<ScalarPoint3f> intersections;
            
            //  TODO: Fix tangent calculation for alpha = (-) pi / 2
            if (alpha == dr::Pi<ScalarFloat> / 2.0f || alpha == -dr::Pi<ScalarFloat> / 2.0f) {
                //  If alpha is pi/2, we have a vertical ray. The intersection
                //  with the spherical shell is then given by y^2 = r^2
                for (int32_t r_idx = 1; r_idx <= (int32_t) resolution.z(); ++r_idx) {
                    ScalarFloat r = r_idx * voxel_size.z();
                    ScalarFloat z = r;

                    //  Calculate the intersection points with the spherical shell
                    intersections.push_back(ScalarPoint3f{p.x(), p.y(), z});
                    intersections.push_back(ScalarPoint3f{p.x(), p.y(), -z});
                }
            } else {
                //  Since coefficient a is constant for all radii, we can precompute it
                ScalarFloat a = (1 + (tan_alpha * tan_alpha));

                //  Loop over all radii
                for (int32_t r_idx = 1; r_idx <= (int32_t) resolution.z(); ++r_idx) {
                    ScalarFloat r = r_idx * voxel_size.z();

                    //  Calculate the intersection point with the spherical shell 
                    ScalarFloat b = 2.0f * tan_alpha * r_max;
                    ScalarFloat c = r_max * r_max - r * r;
                    
                    //Log(Warn, "alpha = ", alpha * to_deg, ", tan(alpha) = ", tan_alpha);
                    //Log(Warn, "     r_max = ", r_max, ", r = ", r, ", b = ", b, ", c = ", c);

                    ScalarFloat discriminant = b * b - 4 * a * c;
                    //Log(Warn, "     discriminant = ", discriminant);

                    if (discriminant < 0) {
                        // No intersection with the spherical shell
                        continue;   
                    }

                    if (discriminant == 0) {
                        // One intersection point, calculate it
                        ScalarFloat t = -b / (2 * a);
                        
                        // Intersection behind the ray origin?
                        if (t < 0) 
                            continue; 

                        ScalarFloat x = t;
                        ScalarFloat z = r_max + x * tan_alpha;
                        
                        // Add the intersection point to the list
                        intersections.push_back(ScalarPoint3f{x, p.y(), z});
                    } else {
                        // Two intersection points, take the first one
                        ScalarFloat sqrt_discriminant = dr::sqrt(discriminant);
                        ScalarFloat t1 = (-b + sqrt_discriminant) / (2 * a);
                        ScalarFloat t2 = (-b - sqrt_discriminant) / (2 * a);
                        
                        // TODO: Check if this can actually occur?
                        // Both intersections behind the ray origin?
                        if (t1 < 0 && t2 < 0) 
                            continue;

                        ScalarFloat x_1 = std::min(t1, t2);
                        ScalarFloat z_1 = r_max + x_1 * tan_alpha;

                        ScalarFloat x_2 = std::max(t1, t2);
                        ScalarFloat z_2 = r_max + x_2 * tan_alpha;

                        //  Add the intersection points to the list
                        intersections.push_back(ScalarPoint3f{x_1, p.y(), z_1});
                        intersections.push_back(ScalarPoint3f{x_2, p.y(), z_2});    
                    }
                }
            }

            //  Sort the intersection points by z-coordinate
            std::sort(intersections.begin(), intersections.end(),
                      [](const ScalarPoint3f &a, const ScalarPoint3f &b) {
                          return a.z() > b.z();
                      });

            compute_cdf(intersections);
        }
    }

    void compute_cdf(const std::vector<ScalarPoint3f> intersections) const {
        static constexpr size_t SpectralSize = dr::array_size_v<UnpolarizedSpectrum>;
        using ScalarUnpolarized = dr::Array<dr::scalar_t<UnpolarizedSpectrum>, SpectralSize>;

        MediumInteraction3f mei = dr::zeros<MediumInteraction3f>();

        for (size_t i = 0; i < intersections.size() - 1; i++) {
            ScalarPoint3f p1 = intersections[i];
            ScalarPoint3f p2 = intersections[i + 1];

            //  Calculate the length of the segment
            ScalarFloat segment_length = dr::norm(p2 - p1);

            //  Calculate the middle point to use as sample point 
            //  to get the medium interaction
            ScalarPoint3f sample_point = (p1 + p2) * 0.5f;

            //  Set the sample point z-coordinate to a POSITIVE value
            //  to ensure we are sampling the "spherical shell" medium
            //  (medium is not defined below the z-coordinate of the bounding box)
            sample_point.z() = std::abs(sample_point.z());

            //  Set the medium interaction point
            mei.p = sample_point;

            ScalarUnpolarized cumulative = dr::zeros<ScalarUnpolarized>();
            std::vector<ScalarFloat> cum_opt_thickness(intersections.size() * SpectralSize);

            //  Get the scattering coefficients at the sample point
            std::tie(mei.sigma_s, mei.sigma_n, mei.sigma_t) =
                get_scattering_coefficients(mei, true);

            //  Calculate the optical thickness per wavelength channel
            for (size_t k = 0; k < SpectralSize; ++k) {
                if constexpr (dr::is_jit_v<UnpolarizedSpectrum>)
                    cumulative[k] += ScalarFloat(mei.sigma_t[k][0]) * segment_length;
                else
                    cumulative[k] += ScalarFloat(mei.sigma_t[k]) * segment_length;

                Log(Warn, "Cumulative optical thickness at index ", i, " and channel ", k,
                    " = ", cumulative[k]);

                cum_opt_thickness[i * SpectralSize + k] = cumulative[k];
            }
        }
    }

    UnpolarizedSpectrum get_majorant(const MediumInteraction3f &mi,
                                     Mask active) const override {
    }

    std::tuple<UnpolarizedSpectrum, UnpolarizedSpectrum, UnpolarizedSpectrum>
    get_scattering_coefficients(const MediumInteraction3f &mi,
                                Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::MediumEvaluate, active);

        auto sigmat = m_sigmat->eval(mi, active);
        if (has_flag(m_phase_function->flags(), PhaseFunctionFlags::Microflake))
            sigmat *= m_phase_function->projected_area(mi, active);

        auto sigmas = sigmat * m_albedo->eval(mi, active);
        auto sigman = m_max_density - sigmat;

        return { sigmas, sigman, sigmat };
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