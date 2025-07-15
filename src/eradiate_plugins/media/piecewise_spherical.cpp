#include <drjit/if_stmt.h>
#include <drjit/texture.h>
#include <drjit/array.h>
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

    static constexpr ScalarFloat PI_HALF    = dr::Pi<Float> / 2.0f;

    using ScalarIndex           = uint32_t;
    using ScalarSize            = uint32_t;
    using FloatStorage          = DynamicBuffer<Float>;
    using ShellIntersection     = std::tuple<int32_t, int32_t, ScalarPoint3f>;
    using CumulativeOTEntry     = std::tuple<int32_t, int32_t, int32_t, ScalarPoint3f, ScalarFloat>;

    PiecewiseSphericalMedium(const Properties &props) : Base(props) {
        m_is_homogeneous = false;
        m_albedo         = props.volume<Volume>("albedo", 0.75f);
        m_sigmat         = props.volume<Volume>("sigma_t", 1.f);
        m_angle_samples  = props.get<int32_t>("angular_samples", 0);
        m_has_spectral_extinction =
            props.get<bool>("has_spectral_extinction", true);
        m_max_density    = dr::opaque<Float>(m_sigmat->max());

        //Log(Warn, "sigma_t bbox: %s",
        //    m_sigmat->bbox());

        //Log(Warn, "sigma_t res: %s", 
        //    m_sigmat->resolution());

        //Log(Warn, "sigma_t voxel size: %s",
        //    m_sigmat->voxel_size());

        //  Precompute OT
        precompute();
    }

    std::tuple<typename Medium<Float, Spectrum>::MediumInteraction3f, Float, Float>
    sample_interaction_real(const Ray3f &ray, const SurfaceInteraction3f &si,
                            Float sample, UInt32 channel,
                            Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::MediumSample, active);

        ScalarFloat two_pi = 2.0f * dr::Pi<Float>;

        // For now, run everything in scalar mode
        for (size_t i = 0; i < dr::width(sample); ++i) {
            ScalarVector3f ray_o = dr::slice(ray.o, i);
            ScalarVector3f ray_dir = dr::normalize(dr::slice(ray.d, i));
            
            // Custom ray/angle intersection, since we can not rely on the 
            // medium bounding box for spherical media.
            ScalarFloat theta = dr::atan2(ray_dir.x(), ray_dir.z()) + two_pi;
            ScalarFloat alpha = (theta > two_pi ? theta - two_pi : theta) - dr::Pi<ScalarFloat>;

            std::vector<ShellIntersection> intersections = compute_ray_intersections(ray_o, alpha);
        
            //Log(Warn, "Origin = %f", ray_o);
            //Log(Warn, "Alpha = %f", alpha * 180.0f / dr::Pi<ScalarFloat>);

            for (uint32_t j = 0; j < intersections.size(); j++) {}
                //Log(Warn, "Intersection %d: r_idx = %d, intersection_idx = %d, p = (%f, %f, %f)",
                //    j, std::get<0>(intersections[j]), std::get<1>(intersections[j]),
                //    std::get<2>(intersections[j]).x(), std::get<2>(intersections[j]).y(),
                //    std::get<2>(intersections[j]).z());
        }

        return { dr::zeros<MediumInteraction3f>(), dr::zeros<Float>(), dr::zeros<Float>() };
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
        callback->put_parameter("cum_opt_thickness", m_cum_opt_thickness.tensor(),
                                +ParamFlags::NonDifferentiable);
        callback->put_parameter("angles", m_angles,
                                +ParamFlags::NonDifferentiable);
        Base::traverse(callback);
    }

    void parameters_changed(const std::vector<std::string> & /*keys*/ = {}) override {
        precompute();
    }

    void precompute() const {
        //  Precompute the different radii of the spherical shells
        precompute_radii();

        //  Precompute the ray angles w.r.t. NADIR
        precompute_directions();

        //  Precompute the (cumulative) optical thickness
        precompute_optical_thickness();
    }

    bool is_almost_zero(ScalarFloat value) const {
        return dr::abs(value) < 1e-04f;
    }

    void precompute_radii() const {
        const ScalarVector3i resolution = m_sigmat->resolution();
        //const ScalarVector3f voxel_size = m_sigmat->voxel_size();

        //  Temporary workaround till Nicolae is back
        const ScalarVector3f voxel_size = m_sigmat->bbox().max / 
            ScalarVector3f(resolution.x(), resolution.y(), resolution.z());

        // Check that the first two dimensions are equal to 1.
        if (resolution.x() > 1 || resolution.y() > 1)
            Throw("PiecewiseMedium: x or y resolution bigger than one, assumed "
                  "shape is [1,1,n]");

        std::vector<ScalarFloat> radii;

        //  Based on the resolution and voxel size, we compute the radii of the 
        //  encompassing spherical shell geometry.
        for (int32_t r_idx = 1; r_idx <= resolution.z(); ++r_idx)
            radii.push_back(r_idx * voxel_size.z());

        //  Save the radius table
        m_radii = dr::load<FloatStorage>(radii.data(), radii.size());
    }

    void precompute_directions() const {
        //  Add one to account for NADIR sample (which is always included)
        int32_t sided_samples = (m_angle_samples != 0) ? m_angle_samples + 1 : 0;

        //  Storage for the direction vectors
        std::vector<ScalarFloat> angles;

        if (sided_samples != 0) {
            ScalarFloat step = PI_HALF / sided_samples;

            //  Compute angles theta for sampling
            for (int32_t angle_idx = 1; angle_idx < sided_samples; ++angle_idx) {
                //  Alpha is simply PI_HALF - (angle_idx * step)
                ScalarFloat alpha_left = PI_HALF - (angle_idx * step);
                ScalarFloat alpha_right = PI_HALF + (angle_idx * step);              

                angles.push_back(alpha_left);
                angles.push_back(alpha_right);
            }
        }

        //  Always push alpha for NADIR (= -PI_HALF)
        angles.push_back(PI_HALF);

        //  Sort angles in ascending order
        std::sort(angles.begin(), angles.end());

        m_angles = dr::load<FloatStorage>(angles.data(), angles.size());
    }

    std::vector<ShellIntersection> compute_ray_sphere_intersection(ScalarPoint3f p, ScalarFloat alpha, int32_t r_idx) const {
        std::vector<ShellIntersection> intersections;
        ScalarFloat r_max   = m_radii[m_radii.size() - 1];  //  The maximum radius is the last element in the radii vector
        ScalarFloat radius  = m_radii[r_idx];               //  The radius of the spherical shell we are currently processing

        //Log(Warn, "Computing intersection for point %f with alpha %f and radius %f",
        //    p, alpha * 180.0f / dr::Pi<ScalarFloat>, radius);

        if (alpha == PI_HALF) {
            //  Calculate the intersection points with the spherical shell
            intersections.emplace_back(r_idx, 0, ScalarPoint3f{p.x(), p.y(), radius});
            intersections.emplace_back(r_idx, 1, ScalarPoint3f{p.x(), p.y(), -radius});
        } else {
            ScalarFloat tan_alpha   = dr::tan(alpha);

            //  Since coefficient a is constant for all radii, we can precompute it.
            //  Also, geometric identity: 1 + tan^2(alpha) = sec^2(alpha)
            ScalarFloat a = 1 + dr::square(tan_alpha);

            //  Calculate the intersection point with the spherical shell 
            ScalarFloat b = 2.0f * tan_alpha * r_max;
            ScalarFloat c = (r_max * r_max) - (radius * radius);
            ScalarFloat discriminant = (b * b) - (4 * a * c);

            if (is_almost_zero(discriminant)) {
                // One intersection point
                ScalarFloat t = -b / (2 * a);

                ScalarFloat x = t;
                ScalarFloat z = r_max + x * tan_alpha;

                intersections.emplace_back(r_idx, 0, ScalarPoint3f{x, p.y(), z});
            } else if (discriminant > 0) {
                // Two intersection points
                ScalarFloat sqrt_discriminant = dr::sqrt(discriminant);
                ScalarFloat t1 = (-b + sqrt_discriminant) / (2 * a);
                ScalarFloat t2 = (-b - sqrt_discriminant) / (2 * a);

                ScalarFloat x_1 = std::min(t1, t2);
                ScalarFloat z_1 = r_max + x_1 * tan_alpha;

                ScalarFloat x_2 = std::max(t1, t2);
                ScalarFloat z_2 = r_max + x_2 * tan_alpha;

                ScalarPoint3f p1{x_1, p.y(), z_1};
                ScalarPoint3f p2{x_2, p.y(), z_2};

                // Compute distances for r_idx assignment
                if (dr::norm(p1 - p) < dr::norm(p2 - p)) {
                    //  p1 is closest
                    intersections.emplace_back(r_idx, 0, p1);
                    intersections.emplace_back(r_idx, 1, p2);
                } else {
                    //  p2 is closest
                    intersections.emplace_back(r_idx, 0, p2);
                    intersections.emplace_back(r_idx, 1, p1);
                }
            }
        }

        return intersections;
    }

    std::vector<ShellIntersection> compute_ray_intersections(ScalarPoint3f p, ScalarFloat alpha) const {
        std::vector<ShellIntersection> intersections;   
        
        for (int32_t r_idx = 0; r_idx < (int32_t) m_radii.size(); ++r_idx) {
            //  Compute the intersection points with the spherical shell
            std::vector<ShellIntersection> shell_intersections = compute_ray_sphere_intersection(p, alpha, r_idx);

            //  Data movement
            intersections.insert(intersections.end(), shell_intersections.begin(), shell_intersections.end());
        }

        return intersections;
    }

    void sort_intersections(std::vector<ShellIntersection> &intersections) const {
        //  Sort the intersection points by distance to the ray origin
        std::sort(intersections.begin(), intersections.end(),
                  [](const auto &a, const auto &b) {
                        auto p1 = std::get<2>(a);
                        auto p2 = std::get<2>(b);
                        return p1.z() > p2.z();
                    });
    }

    std::vector<CumulativeOTEntry> compute_cumulative_ot(const std::vector<ShellIntersection> intersections) const {
        static constexpr size_t SpectralSize = dr::size_v<UnpolarizedSpectrum>;
        using ScalarUnpolarized = dr::Array<dr::scalar_t<UnpolarizedSpectrum>, SpectralSize>;

        ScalarUnpolarized cumulative = dr::zeros<ScalarUnpolarized>();
        std::vector<CumulativeOTEntry> cum_ot(intersections.size() * SpectralSize);
        
        // Not enough intersections to compute cumulative OT -> Either at the 
        // edge of the outer shell (which should not happen given the sample origin)
        // or the ray is completely outside the medium (which should also
        // not happen in the precomputation scenario).
        if (intersections.size() < 2)
            return cum_ot;

        //  The first intersection is always the top shell, 
        //  so we can use it to initialize the CDF
        ShellIntersection first_intersection = intersections[0];
        int32_t r_idx               = std::get<0>(first_intersection);
        int32_t intersection_idx    = std::get<1>(first_intersection);
        ScalarPoint3f p1            = std::get<2>(first_intersection);
        
        for (size_t k = 0; k < SpectralSize; ++k) 
            cum_ot[k] = std::make_tuple(r_idx, intersection_idx, k, p1, 0.0f);

        //  Iterate over the intersections and calculate the optical thickness
        //  for each segment between two intersections.
        MediumInteraction3f mei = dr::zeros<MediumInteraction3f>();
        for (size_t i = 1; i < intersections.size(); ++i) {
            ShellIntersection primary   = intersections[i - 1];
            ShellIntersection secondary = intersections[i];

            ScalarPoint3f p1 = std::get<2>(primary);
            ScalarPoint3f p2 = std::get<2>(secondary);

            //  Get the metadata of the SECOND intersection, for which we have
            //  to store the data in the CDF
            int32_t r_idx               = std::get<0>(secondary);
            int32_t intersection_idx    = std::get<1>(secondary);

            //  Calculate the length of the segment
            ScalarFloat segment_length = dr::norm(p2 - p1);

            //  Calculate the middle point to use as sample point 
            //  to get the medium interaction
            ScalarPoint3f mid_point = (p1 + p2) * 0.5f;

            //  Exploit the fact that the medium is spherical and sample a point
            //  has only to be characterized by it's norm towards the origin.
            //  We add a small epsilon to avoid hitting shell boundaries exactly.
            ScalarPoint3f sample_point{0.0f, 0.0f, dr::abs(dr::norm(mid_point)) + 1e-04f};

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

                cum_ot[i * SpectralSize + k] = std::make_tuple(
                    r_idx, intersection_idx, k, p1, cumulative[k]
                );
            }
        }

        return cum_ot;
    }

    void precompute_optical_thickness() const {
        static constexpr size_t SpectralSize = dr::size_v<UnpolarizedSpectrum>;
        using ScalarUnpolarized = dr::Array<dr::scalar_t<UnpolarizedSpectrum>, SpectralSize>;

        const ScalarVector3i resolution = m_sigmat->resolution();

        // Check that the first two dimensions are equal to 1.
        if (resolution.x() > 1 || resolution.y() > 1)
            Throw("PiecewiseMedium: x or y resolution bigger than one, assumed "
                  "shape is [1,1,n]");

        //  The interaction point from where we calculate directions
        ScalarPoint3f medium_max = m_sigmat->bbox().max;
        ScalarPoint3f p{0.0f, 0.0f, medium_max.z()};
    
        //  The maximum number of intersection is the resolution * the number of angles *
        //  the spectral size (if spectral data is used)
        int32_t max_intersections = resolution.z() * 2;
        std::vector<ScalarFloat> opt_thickness_data(m_angles.size() * max_intersections * SpectralSize);

        //  Here we assume we start at the top shell
        for (int32_t i = 0; i < (int32_t) m_angles.size(); ++i) {
            //  We assume that the voxel size is the radius of a spherical shell!
            //  The intersection point with the shell from P can then be found by 
            //  using the equation of a circle (x^2 + y^2 = r^2) and of the line
            //  defining the directional through P (y = theta * x + r_max) where 
            //  r_max is defined by the z-coordinate of the max bounding box point.
            ScalarFloat alpha = m_angles[i];

            std::vector<ShellIntersection> intersections = compute_ray_intersections(p, alpha);

            //  Sort the intersection points by distance to the ray origin
            sort_intersections(intersections);

            //  Compute cumulative OT from the intersection points
            std::vector<CumulativeOTEntry> cum_ot = compute_cumulative_ot(intersections);

            //  Process the cumulative OT for proper texture indexing
            for (size_t j = 0; j < cum_ot.size(); ++j) {
                CumulativeOTEntry entry     = cum_ot[j];
                int32_t r_idx               = std::get<0>(entry);
                int32_t intersection_idx    = std::get<1>(entry);
                int32_t spectral_idx        = std::get<2>(entry);
                ScalarFloat cum_ot_value    = std::get<4>(entry);

                //  Compute the shell index
                size_t shell_idx = (intersection_idx == 0) 
                                    ? (max_intersections / 2 - r_idx) 
                                    : (max_intersections / 2 + r_idx - 1); 

                //  Compute the flat index based on the shell index, the 
                //  max number of intersections, the angle index and 
                //  the spectral index.
                size_t index = spectral_idx * (m_angles.size() * max_intersections)
                                + i * max_intersections
                                + shell_idx;
                
                //  Store the cumulative OT value in the data array
                opt_thickness_data[index] = cum_ot_value;
            }
        }

        size_t shape[3] = { SpectralSize, m_angles.size(), static_cast<size_t>(max_intersections) };

        //  Store the cumulative OT as a texture
        m_cum_opt_thickness = Texture2f(TensorXf(opt_thickness_data.data(), 3, shape), true, 
            true, dr::FilterMode::Linear, dr::WrapMode::Clamp);
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

    mutable FloatStorage m_radii;
    mutable FloatStorage m_angles;

    Float m_max_density;
    mutable Texture2f m_cum_opt_thickness;
};

MI_IMPLEMENT_CLASS_VARIANT(PiecewiseSphericalMedium, Medium)
MI_EXPORT_PLUGIN(PiecewiseSphericalMedium, "Piecewise Spherical Medium")
NAMESPACE_END(mitsuba)