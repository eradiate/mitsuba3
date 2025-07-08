#include "drjit/array_router.h"
#include <drjit/texture.h>
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
#include <mitsuba/render/texture.h>
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
    using ShellIntersection = std::tuple<int32_t, int32_t, ScalarPoint3f>;
    using CDFEntry          = std::tuple<int32_t, int32_t, ScalarPoint3f, int32_t, ScalarFloat>;

    static constexpr size_t SpectralSize = dr::array_size_v<UnpolarizedSpectrum>;
    using ScalarUnpolarized = dr::Array<dr::scalar_t<UnpolarizedSpectrum>, SpectralSize>;

    PiecewiseSphericalMedium(const Properties &props) : Base(props) {
        m_is_homogeneous = false;
        m_albedo         = props.volume<Volume>("albedo", 0.75f);
        m_sigmat         = props.volume<Volume>("sigma_t", 1.f);
        m_angle_samples  = props.get<int32_t>("samples", 0);

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
        MI_MASKED_FUNCTION(ProfilerPhase::MediumSample, active);

        // Initial intersection with the medium
        auto [aabb_its, mint, maxt] = intersect_aabb(ray);
        aabb_its &= (dr::isfinite(mint) || dr::isfinite(maxt));
        active &= aabb_its;
        dr::masked(mint, !active) = 0.f;
        dr::masked(maxt, !active) = dr::Infinity<Float>;
    }

    std::tuple<Float, Float, Mask>
    eval_transmittance_pdf_real(const Ray3f &ray,
                                const SurfaceInteraction3f &si, UInt32 channel,
                                Mask active) const override {
    }

    void traverse(TraversalCallback *callback) override {
        callback->put_object("albedo", m_albedo.get(),
                             +ParamFlags::Differentiable);
        callback->put_object("sigma_t", m_sigmat.get(),
                             +ParamFlags::Differentiable);
        callback->put_parameter("optical_thickness_cdf", m_opt_thickness_cdf.tensor(),
                                +ParamFlags::NonDifferentiable);
        callback->put_parameter("angles",
                                m_angles,
                                +ParamFlags::NonDifferentiable);
        Base::traverse(callback);
    }

    void parameters_changed(const std::vector<std::string> & /*keys*/ = {}) override {
    }

    void precompute() const override {
        std::vector<ScalarFloat> angles = precompute_directions();
        precompute_optical_thickness(angles);
    }

    std::vector<ScalarFloat> precompute_directions() const {
        //  Add one to account for NADIR sample (which is always included)
        int32_t sided_samples = (m_angle_samples != 0) ? m_angle_samples + 1 : 0;

        //  Storage for the direction vectors
        std::vector<ScalarFloat> angles;

        ScalarFloat PI_HALF = dr::Pi<ScalarFloat> / 2.0f;
        ScalarFloat NADIR = -PI_HALF; // Nadir direction (downwards)

        //  Always push NADIR
        angles.push_back(NADIR + PI_HALF); // Convert to alpha angle

        if (sided_samples != 0) {
            ScalarFloat step = PI_HALF / sided_samples;

            //  Compute angles theta for sampling
            for (int32_t angle_idx = 1; angle_idx < sided_samples; angle_idx++) {
                //  Based on step, compute angle
                ScalarFloat theta_left = NADIR - (angle_idx * step);
                //ScalarFloat theta_right = NADIR + (angle_idx * step);
                
                //  Ensure that the angles are within valid range
                //if (theta_left < -2.0f * PI_HALF || theta_right > 0.0f)
                if (theta_left < -2.0f * PI_HALF)
                    continue; // Skip angles outside the valid range

                ScalarFloat alpha_left = theta_left + PI_HALF;
                //ScalarFloat alpha_right = theta_right + PI_HALF;

                angles.push_back(alpha_left);
                //angles.push_back(alpha_right);
            }
        }

        //  Sort angles in ascending order
        std::sort(angles.begin(), angles.end());

        ScalarFloat to_rad = dr::Pi<ScalarFloat> / 180.0f;

        std::vector<ScalarFloat> test_angles;
        test_angles.push_back(15.0f * to_rad); // Add 15 degrees for testing
        test_angles.push_back(75.0f * to_rad); // Add 75 degrees for testing

        m_angles = dr::load<FloatStorage>(test_angles.data(), test_angles.size()); //dr::load<FloatStorage>(angles.data(), angles.size());

        return test_angles;
    }

    void precompute_optical_thickness(const std::vector<ScalarFloat> angles) const {
        // Check that the first two dimensions are equal to 1 and that z is one
        // or more.
        const ScalarVector3i resolution = m_sigmat->resolution();
        const ScalarVector3f voxel_size = m_sigmat->voxel_size();
        if (resolution.x() > 1 || resolution.y() > 1)
            Throw("PiecewiseMedium: x or y resolution bigger than one, assumed "
                  "shape is [1,1,n]");

        //  The interaction point from where we calculate directions
        ScalarPoint3f p         = m_sigmat-> bbox().max;
        ScalarFloat r_max       = p.z();
    
        //  The maximum number of intersection is the resolution * the number of angles *
        //  the spectral size (if spectral data is used)
        int32_t max_intersections = resolution.z() * 2;
        std::vector<ScalarFloat> opt_thickness_data(angles.size() * max_intersections * SpectralSize);

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

            std::vector<ShellIntersection> intersections;        

            ScalarFloat angle_deg = ((-dr::Pi<ScalarFloat> / 2.0f) - ((dr::Pi<ScalarFloat> / 2.0f) - alpha)) * to_deg;

            Log(Warn, "Angle = ", angle_deg, ", (alpha = ", alpha, ")");

            if (alpha == 0.0f) {                
                //  If alpha is 0, we have a vertical ray. The intersection
                //  with the spherical shell is then given by y^2 = r^2
                
                //  Loop over all radii from largest to smallest to account for CDF indexing
                for (int32_t r_idx = (int32_t) resolution.z(); r_idx >= 1; --r_idx) {
                    ScalarFloat r = r_idx * voxel_size.z();
                    ScalarFloat z = r;

                    ScalarPoint3f p1{p.x(), p.y(), z};
                    ScalarPoint3f p2{p.x(), p.y(), -z};

                    //  Calculate the intersection points with the spherical shell
                    intersections.emplace_back(r_idx, 0, ScalarPoint3f{p.x(), p.y(), z});
                    intersections.emplace_back(r_idx, 1, ScalarPoint3f{p.x(), p.y(), -z});
                }
            } else {
                //  Since coefficient a is constant for all radii, we can precompute it.
                //  Also, geometric identity: 1 + tan^2(alpha) = sec^2(alpha)
                ScalarFloat a = dr::sqr(dr::sec(alpha));

                //  Loop over all radii from largest to smallest to account for CDF indexing
                for (int32_t r_idx = (int32_t) resolution.z(); r_idx >= 1; --r_idx) { 
                    ScalarFloat r = r_idx * voxel_size.z();

                    //  Calculate the intersection point with the spherical shell 
                    ScalarFloat b = 2.0f * tan_alpha * r_max;
                    ScalarFloat c = (r_max * r_max) - (r * r);
                    ScalarFloat discriminant = (b * b) - (4 * a * c);

                    if (discriminant < 0) {
                        // No intersection with the spherical shell
                        continue;   
                    }

                    //  THIS IS PROBLEMATIC -> Store at 0 index for r_idx for now
                    else if (discriminant == 0) {
                        // One intersection point, calculate it
                        ScalarFloat t = -b / (2 * a);
                        
                        // Intersection behind the ray origin?
                        if (t < 0) 
                            continue; 

                        ScalarFloat x = t;
                        ScalarFloat z = r_max + x * tan_alpha;
                        
                        // Add the intersection point to the list
                        intersections.emplace_back(r_idx, 0, ScalarPoint3f{x, p.y(), z});
                    } 
                    
                    else {
                        // Two intersection points, take the first one
                        ScalarFloat sqrt_discriminant = dr::sqrt(discriminant);
                        ScalarFloat t1 = (-b + sqrt_discriminant) / (2 * a);
                        ScalarFloat t2 = (-b - sqrt_discriminant) / (2 * a);

                        ScalarFloat x_1 = std::min(t1, t2);
                        ScalarFloat z_1 = r_max + x_1 * tan_alpha;

                        ScalarFloat x_2 = std::max(t1, t2);
                        ScalarFloat z_2 = r_max + x_2 * tan_alpha;

                        ScalarPoint3f p1{x_1, p.y(), z_1};
                        ScalarPoint3f p2{x_2, p.y(), z_2};

                        //  Add the intersection points to the list
                        intersections.emplace_back(r_idx, 0, ScalarPoint3f{x_1, p.y(), z_1});
                        intersections.emplace_back(r_idx, 1, ScalarPoint3f{x_2, p.y(), z_2});       
                    }
                }
            }

            //  Sort the intersection points by distance to the ray origin
            std::sort(intersections.begin(), intersections.end(),
                      [](const auto &a, const auto &b) {
                            auto p1 = std::get<2>(a);
                            auto p2 = std::get<2>(b);
                            return p1.z() > p2.z();
                        });

            std::vector<CDFEntry> cdf = compute_cdf(intersections);

            //  Once we obtained the CDF, we can insert it into the array
            //  that will be used for the texture
            for (size_t j = 0; j < cdf.size(); ++j) {
                CDFEntry entry = cdf[j];
                int32_t r_idx               = std::get<0>(entry);
                int32_t intersection_idx    = std::get<1>(entry);
                ScalarFloat cdf_value       = std::get<4>(entry);

                //  Compute the shell index
                size_t cdf_index = (intersection_idx == 0) ? (max_intersections / 2 - r_idx) : (max_intersections / 2 + r_idx - 1); 
                
                //  Compute the final index basedd on the spectral channel and the angle
                //  to obtain a linearized 3D index
                size_t index = cdf_index + (i * max_intersections);

                //  Store the CDF value in the data array
                opt_thickness_data[index] = cdf_value;
            }
        }

        size_t shape[3] = { angles.size(), static_cast<size_t>(max_intersections), 1 };

        //  Publish the optical thickness CDF as a texture
        m_opt_thickness_cdf = Texture2f(TensorXf(opt_thickness_data.data(), 3, shape), true, 
            true, dr::FilterMode::Linear, dr::WrapMode::Clamp);
    }

    std::vector<CDFEntry> compute_cdf(const std::vector<ShellIntersection> intersections) const {
        MediumInteraction3f mei = dr::zeros<MediumInteraction3f>();
        ScalarUnpolarized cumulative = dr::zeros<ScalarUnpolarized>();
        std::vector<CDFEntry> opt_thickness_cdf(intersections.size() * SpectralSize);
        
        if (intersections.size() < 2) {
            Log(Warn, "Not enough intersections to compute CDF");
            return opt_thickness_cdf; // Not enough intersections to compute CDF
        } 

        for (size_t i = 0; i < intersections.size() - 1; i++) {
            ShellIntersection primary = intersections[i];
            ShellIntersection secondary = intersections[i + 1];

            ScalarPoint3f p1 = std::get<2>(primary);
            ScalarPoint3f p2 = std::get<2>(secondary);
            
            //  Get the metadata of the first intersection, for which we have
            //  to store the data in the CDF
            int32_t r_idx = std::get<0>(primary);
            int32_t intersection_idx = std::get<1>(primary);

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

            //  Get the scattering coefficients at the sample point
            std::tie(mei.sigma_s, mei.sigma_n, mei.sigma_t) =
                get_scattering_coefficients(mei, true);
        
            //  Calculate the optical thickness per wavelength channel
            for (size_t k = 0; k < SpectralSize; ++k) {
                if constexpr (dr::is_jit_v<UnpolarizedSpectrum>)
                    cumulative[k] += ScalarFloat(mei.sigma_t[k][0]) * segment_length;
                else 
                    cumulative[k] += ScalarFloat(mei.sigma_t[k]) * segment_length;

                opt_thickness_cdf[i * SpectralSize + k] = std::make_tuple(
                    r_idx, intersection_idx, p1, k, cumulative[k]
                );
            }
        }

        //  At the last intersection, we need to add an element to the CDF with the
        //  optical thickness of the last segment (no more medium after this, but is
        //  necessary for the CDF to be complete and to allow for correct sampling
        //  beyond the last intersection)
        auto last_intersection = intersections.back();
        ScalarPoint3f p1 = std::get<2>(last_intersection);
        int32_t r_idx = std::get<0>(last_intersection);
        int32_t intersection_idx = std::get<1>(last_intersection);

        for (size_t k = 0; k < SpectralSize; ++k) {
            opt_thickness_cdf[(intersections.size() - 1) * SpectralSize + k] = 
                std::make_tuple(r_idx, intersection_idx, p1, k, cumulative[k]);
        }

        return opt_thickness_cdf;
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
    int32_t m_angle_samples;
    mutable FloatStorage m_angles;

    Float m_max_density;
    mutable Texture2f m_opt_thickness_cdf;

    //FloatStorage m_cum_opt_thickness;
    //FloatStorage m_reverse_cum_opt_thickness;
};

MI_IMPLEMENT_CLASS_VARIANT(PiecewiseSphericalMedium, Medium)
MI_EXPORT_PLUGIN(PiecewiseSphericalMedium, "Piecewise Spherical Medium")
NAMESPACE_END(mitsuba)