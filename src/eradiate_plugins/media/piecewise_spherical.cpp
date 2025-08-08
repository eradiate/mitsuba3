#include "drjit/while_loop.h"
#include "mitsuba/core/logger.h"
#include "mitsuba/core/math.h"
#include <algorithm>
#include <cstdint>
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
#include <utility>

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
    using IntStorage            = DynamicBuffer<uint32_t>;
    using ShellIntersection     = std::tuple<int32_t, int32_t, ScalarPoint3f>;
    using CumulativeOTEntry     = std::tuple<int32_t, int32_t, int32_t, ScalarPoint3f, ScalarFloat>;
    using Index                 = dr::uint32_array_t<Float>;

    PiecewiseSphericalMedium(const Properties &props) : Base(props) {
        m_is_homogeneous = false;
        m_albedo         = props.volume<Volume>("albedo", 0.75f);
        m_sigmat         = props.volume<Volume>("sigma_t", 1.f);
        m_angle_samples  = props.get<int32_t>("angular_samples", 0);
        m_angle_step     = m_angle_step = PI_HALF / ((m_angle_samples != 0) ? m_angle_samples + 1 : 0);
        m_has_spectral_extinction =
            props.get<bool>("has_spectral_extinction", true);
        m_max_density    = dr::opaque<Float>(m_sigmat->max());

        m_max_intersections = m_sigmat->resolution().z() * 2;
        m_medium_radius = m_sigmat->bbox().max.z();

        m_texture_cell_angle_step   = 1.0f / (m_angle_samples * 2 + 1);
        m_texture_cell_shell_step   = 1.0f / m_max_intersections;
        m_texture_cell_angle_offset = 0.5f * m_texture_cell_angle_step;
        m_texture_cell_shell_offset = 0.5f * m_texture_cell_shell_step;

        //  Precompute OT
        precompute();
    }

    std::vector<ScalarFloat> build_interpolated_cot_table(uint32_t left_angle_idx, uint32_t right_angle_idx,
                                                          ScalarFloat la_coord, ScalarFloat ra_coord,
                                                          ScalarFloat t) const {
        uint32_t left_noi   = m_intersections_per_angle[left_angle_idx];
        uint32_t right_noi  = m_intersections_per_angle[right_angle_idx];
        uint32_t ilut_size  = std::max(left_noi, right_noi); 
        uint32_t bound      = ilut_size / 2;

        std::vector<ScalarFloat> interpolated_cot(ilut_size, 0.0f);

        for (uint32_t i = 0; i < bound; ++i) {
            uint32_t idx = i;
            uint32_t reverse_idx = m_max_intersections - i - 1;

            //  Interpolate the angle coordinate based on t
            ScalarFloat angle_coord = dr::select(left_angle_idx == right_angle_idx, 
                la_coord, (1.0f - t) * la_coord + t * ra_coord);

            //  Shell coordinates of primary and secondary indices are the same
            ScalarFloat primary_shell_coord = m_texture_cell_shell_offset + idx * 
                    m_texture_cell_shell_step;
            ScalarFloat secondary_shell_coord = m_texture_cell_shell_offset + reverse_idx * 
                    m_texture_cell_shell_step;

            //  Sample position
            Point2f primary_sample_position{ primary_shell_coord, angle_coord };
            Point2f secondary_sample_position{ secondary_shell_coord, angle_coord };

            Float primary_ot_sample, secondary_ot_sample;

            //  Texture Sampling
            m_cot.template eval<Float>(primary_sample_position, &primary_ot_sample);
            m_cot.template eval<Float>(secondary_sample_position, &secondary_ot_sample);

            //  Linear interpolation between primaries and secondaries based on the
            //  interpolation factor t (given by the sample position between two angles)
            interpolated_cot[idx] = dr::slice(primary_ot_sample, 0); 
            interpolated_cot[ilut_size - 1 - i] = dr::slice(secondary_ot_sample, 0);
        }

        return interpolated_cot;
    }

    std::tuple<typename Medium<Float, Spectrum>::MediumInteraction3f, Float, Float>
    sample_interaction_real(const Ray3f &ray, const SurfaceInteraction3f &si,
                            Float sample, UInt32 channel,
                            Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::MediumSample, active);

        static constexpr size_t SpectralSize = dr::size_v<UnpolarizedSpectrum>;
        using ScalarUnpolarized = dr::Array<dr::scalar_t<UnpolarizedSpectrum>, SpectralSize>;

        ScalarFloat pi = dr::Pi<Float>;
        ScalarFloat two_pi = 2.0f * pi;

        //  For now, run everything in scalar mode
        for (size_t i = 0; i < dr::width(sample); ++i) {
            ScalarFloat ksi    = dr::slice(sample, i);

            if (ksi == 0.0f)
                return { dr::zeros<MediumInteraction3f>(), 0.0f, 0.0f };

            ScalarFloat neg_log_ksi = -dr::log(1.0f - ksi);
            ScalarVector3f ray_o    = dr::slice(ray.o, i);
            ScalarVector3f ray_dir  = dr::normalize(dr::slice(ray.d, i));
            
            //  Create the medium interaction, transmittance and pdf variables
            MediumInteraction3f mei = dr::zeros<MediumInteraction3f>();
            ScalarFloat sampled_t = 0.0f;
            ScalarFloat tr = 0.0f;
            ScalarFloat pdf = 1.0f;

            //  Custom ray/angle intersection, since we can not rely on the 
            //  medium bounding box for spherical media.
            ScalarFloat theta = dr::atan2(ray_dir.z(), ray_dir.x());
            ScalarFloat alpha = (theta > two_pi ? theta - two_pi : theta) + dr::Pi<ScalarFloat>;

            //  We can do this, because we know the texture would be symmetrical
            //  and we will, in the future, only keep half of it.
            if (alpha > PI_HALF)
                alpha = pi - alpha;

            uint32_t angle_start_idx    = 0;
            uint32_t angle_end_idx      = m_angles.size() - 1;
            
            //  Find the angle greater than the sample. Based on the binary search,
            //  we will always obtain the index of the angle that is greater than 
            // the sample. 
            uint32_t angle_idx = dr::binary_search<uint32_t>(angle_start_idx, angle_end_idx, 
                [&](int32_t idx) {
                    ScalarFloat angle_cmp = m_angles[idx];
                    return angle_cmp < alpha;
                }
            );

            //  Step 1 - Find the angle coordinate based on the angle index
            uint32_t right_angle_idx = 0;
            uint32_t left_angle_idx  = 0;
            ScalarFloat left_angle_coord = 0.0f;
            ScalarFloat right_angle_coord = 0.0f;

            //  Edge case LEFT -> Snap to left-most angle
            if (angle_idx == 0 && m_angles[angle_idx] > alpha) {
                //  Make sure to use offset!
                left_angle_coord    = m_texture_cell_angle_offset;
                right_angle_coord   = m_texture_cell_angle_offset;
            }

            //  Edge case RIGHT -> Snap to right-most angle
            else if (angle_idx == m_angles.size() - 1 && m_angles[angle_idx] < alpha) {
                right_angle_idx     = angle_idx;
                left_angle_idx      = angle_idx;

                left_angle_coord    = 1.0f - m_texture_cell_angle_offset;
                right_angle_coord   = 1.0f - m_texture_cell_angle_offset;
            }

            //  All other cases
            else {
                right_angle_idx     = angle_idx;
                left_angle_idx      = angle_idx - 1;
            
                right_angle_coord   = m_texture_cell_angle_offset + right_angle_idx *
                    m_texture_cell_angle_step;
                left_angle_coord    = m_texture_cell_angle_offset + left_angle_idx *
                    m_texture_cell_angle_step;
            }

            //  Compute interpolation factor t
            ScalarFloat interpolation_t = (left_angle_idx == right_angle_idx) 
                ? 0.0f
                : (alpha - m_angles[left_angle_idx]) / (m_angles[right_angle_idx] - m_angles[left_angle_idx]);

            //  Step 2 - Compute the interpolated cumulative optical thickness search space
            std::vector<ScalarFloat> interpolated_cot = build_interpolated_cot_table(left_angle_idx, right_angle_idx, 
                left_angle_coord, right_angle_coord, interpolation_t);

            //  Step 3 - Binary search to find between which two intersection points
            //           the sample falls. Currently we use a binary search, but 
            //           future investigation might prove linear time search is
            //           sufficient or more effecient.
            uint32_t start_idx   = 0;
            uint32_t end_idx     = interpolated_cot.size() - 1;
            uint32_t index = dr::binary_search<uint32_t>(start_idx, end_idx,
                [&](uint32_t idx) {
                    ScalarFloat ot_value = interpolated_cot[idx];
                    return ot_value < neg_log_ksi;
                });

            uint32_t interpolated_index = dr::select(index == 0, 0, index - 1);
            ScalarFloat cot = interpolated_cot[interpolated_index];
        
            //  Step 4 - Compute the sphere shell with which to intersect
            //           based on the obtained index.
            uint32_t bound = interpolated_cot.size() / 2;
            bool past_midpoint = (interpolated_index >= bound);
            uint32_t reversed_shell_idx = past_midpoint ? 
                (interpolated_cot.size() - 1 - interpolated_index) : interpolated_index;
            uint32_t shell_idx = m_radii.size() - 1 - reversed_shell_idx;

            //  If the initial index is the max possible index we are past the midpoint,
            //  the ray has escaped the medium.
            if (index >= end_idx && past_midpoint) {    
                sampled_t = dr::Infinity<ScalarFloat>;
                mei.t = sampled_t;

                //  Compute transmittance
                tr = dr::exp(-cot);

                //  Compute pdf
                pdf = tr;

                return { mei, tr, pdf };
            }

            //  Step 5 - Intersect the ray with the given shell index and one shell index up
            std::vector<ShellIntersection> intersections = compute_ray_sphere_intersection(ray_o, alpha, shell_idx);

            if (intersections.size() == 0) {
                //  If we have no intersections, this means we sample from one shell upward!
                //  TODO: Add check to not increase the shell index if we are at the top shell.
                //        This normally should never happen, but we can not guarantee it.
                shell_idx += 1;

                //  If we increase the shell index, we must also update the intersections.
                intersections = compute_ray_sphere_intersection(ray_o, alpha, shell_idx);

                //  Finally, we need to obtain the index in the interpolated COT LUT for this shell.
                //  If we jump up a shell, this means we used interpolated (and thus implicit) points
                //  and therefor, based on even/odd index, decrease the interpolated index by one or two.
                if (interpolated_cot.size() % 4 == 0) {
                    // To even
                    if (interpolated_index % 2 == 0) 
                        interpolated_index -= 2;
                    else
                        interpolated_index -= 1;
                } else if (interpolated_cot.size() % 2 == 0) {
                    // To odd
                    if (interpolated_index % 2 == 0) 
                        interpolated_index -= 1;
                    else
                        interpolated_index -= 2;
                }

                cot = interpolated_cot[interpolated_index];

                //  Finally, we update whether we are past the midpoint or not, based on the new
                //  interpolated index.
                past_midpoint = (interpolated_index >= bound);
            }

            //  Safeguard against no intersections, but this should never happen!
            if (intersections.size() != 2)
                return { dr::zeros<MediumInteraction3f>(), 0.0f, dr::zeros<Float>() };

            ShellIntersection primary_shell_intersection    = intersections[0];
            ShellIntersection secondary_shell_intersection  = intersections[1];

            ScalarPoint3f primary_intersection      = std::get<2>(primary_shell_intersection);
            ScalarPoint3f secondary_intersection    = std::get<2>(secondary_shell_intersection);

            //  Step 6 - Calculate the distance, based on the closest shell intersection
            ScalarPoint3f primary_diff      = primary_intersection - ray_o;
            ScalarPoint3f secondary_diff    = secondary_intersection - ray_o;

            ScalarFloat primary_dist_sq     = dr::dot(primary_diff, primary_diff);
            ScalarFloat secondary_dist_sq   = dr::dot(secondary_diff, secondary_diff);
            
            ScalarFloat min_dist    = dr::select(primary_dist_sq < secondary_dist_sq, primary_dist_sq, secondary_dist_sq);
            ScalarFloat max_dist    = dr::select(primary_dist_sq <= secondary_dist_sq, secondary_dist_sq, primary_dist_sq);
            ScalarFloat min_t       = dr::select(past_midpoint, max_dist, min_dist);

            //  Since we have the shell index from which we sample, we can get the
            //  extinction from the precomputed extinction coefficients.
            ScalarFloat extinction = past_midpoint 
                ? m_extinctions[shell_idx + 1]
                : m_extinctions[shell_idx];

            //  Compute the distance traveled through the shell
            ScalarFloat residual_t = (1.0f / extinction) * (neg_log_ksi - cot);
            sampled_t    = residual_t + min_t;

            //  This clipping is currently necessary in some cases. We could not yet find why
            //  sometimes the sampled distance exceeds its maximum distance.
            ScalarFloat next_cot = interpolated_cot[interpolated_index + 1];
            ScalarFloat max_residual_t = (next_cot - cot) / extinction;

            if (residual_t > max_residual_t) {
                Log(Debug, "Clipping distance from %s to %s", sampled_t, max_residual_t);
                sampled_t = min_t + max_residual_t;
            }

            //  Populate the medium interaction
            mei.t = sampled_t;
            mei.p = ray(mei.t);

            //  Compute transmittance and pdf
            tr = dr::exp(-cot - extinction * residual_t);
            pdf = extinction * tr;

            return { mei, tr, pdf };
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
        callback->put_parameter("cum_opt_thickness", m_cot.tensor(),
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

        //  Precompute extinctions for each shell
        precompute_extinctions();

        //  Precompute the (cumulative) optical thickness
        precompute_optical_thickness();
    }

    bool is_almost_zero(ScalarFloat value) const {
        return dr::abs(value) < 1e-08f;
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
            //  Compute angles theta for sampling
            for (int32_t angle_idx = 1; angle_idx < sided_samples; ++angle_idx) {
                //  Alpha is simply PI_HALF - (angle_idx * step)
                ScalarFloat alpha_left = PI_HALF - (angle_idx * m_angle_step);
                ScalarFloat alpha_right = PI_HALF + (angle_idx * m_angle_step);              

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

    //  TODO: This can probably be simplified thanks to the current medium rep.
    void precompute_extinctions() const {
        std::vector<ScalarFloat> extinctions;

        //  Conditions
        ScalarFloat previous_radius = 0.0f;

        for (uint32_t r_idx = 0; r_idx < m_radii.size(); ++r_idx) {
            //  Get the radius of the current and next shell
            ScalarFloat radius = m_radii[r_idx];

            //  Construct sample point
            ScalarPoint3f lower_sample{0.0f, 0.0f, previous_radius};
            ScalarPoint3f upper_sample{0.0f, 0.0f, radius};
            ScalarPoint3f mid_point = (lower_sample + upper_sample) * 0.5f;
            ScalarPoint3f sample_point{0.0f, 0.0f, 2.0f * (dr::abs(dr::norm(mid_point))) - m_medium_radius + 1e-04f};

            //  Create medium interaction point
            MediumInteraction3f mei = dr::zeros<MediumInteraction3f>();
            mei.p = sample_point;

            //  Get the scattering coefficients at the sample point
            std::tie(mei.sigma_s, mei.sigma_n, mei.sigma_t) =
                get_scattering_coefficients(mei, true);

            //  Store the extinction coefficient for the current shell
            extinctions.push_back(dr::slice(Float(mei.sigma_t[0]), 0));

            //  Update the "previous" radius
            previous_radius = radius;
        }

        //  Save the extinction coefficients
        m_extinctions = dr::load<FloatStorage>(extinctions.data(), extinctions.size());
    }

    std::vector<ShellIntersection> compute_ray_sphere_intersection(ScalarPoint3f p, ScalarFloat alpha, int32_t r_idx) const {
        std::vector<ShellIntersection> intersections;
        ScalarFloat r_max   = m_radii[m_radii.size() - 1];  //  The maximum radius is the last element in the radii vector
        ScalarFloat radius  = m_radii[r_idx];               //  The radius of the spherical shell we are currently processing

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

            //  Important, we skip the shell if the discriminant is negative (no intersections)
            //  or very close to zero. This ensures numeric stability and avoids incompatibility
            //  with the texture format (where 2 entries per sphere are expected). 
            if (discriminant > 0 && !is_almost_zero(discriminant)) {
                //  Mitsuba has a built-in quadratic solver...
                std::tuple<Mask, Float, Float> result = math::solve_quadratic(a, b, c);
                ScalarFloat t_1 = dr::slice(std::get<1>(result), 0);
                ScalarFloat t_2 = dr::slice(std::get<2>(result), 0);

                ScalarFloat x_1 = std::min(t_1, t_2);
                ScalarFloat z_1 = r_max + x_1 * tan_alpha;

                ScalarFloat x_2 = std::max(t_1, t_2);
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

    std::vector<ShellIntersection> compute_ray_intersections(ScalarPoint3f p, ScalarFloat alpha, bool by_distance) const {
        std::vector<ShellIntersection> intersections;   
        
        for (int32_t r_idx = 0; r_idx < (int32_t) m_radii.size(); ++r_idx) {
            //  Compute the intersection points with the spherical shell
            std::vector<ShellIntersection> shell_intersections = compute_ray_sphere_intersection(p, alpha, r_idx);

            //  Data movement
            intersections.insert(intersections.end(), shell_intersections.begin(), shell_intersections.end());
        }

        //  Sort the intersection points by distance to the ray origin
        if (by_distance)
            sort_intersections_by_distance(intersections, p);
        else
            sort_intersections(intersections);

        return intersections;
    }

    void sort_intersections(std::vector<ShellIntersection> &intersections) const {
        //  Sort the intersection points by distance to the ray origin
        std::sort(intersections.begin(), intersections.end(),
                  [](const ShellIntersection &a, const ShellIntersection &b) {
                        ScalarPoint3f p1 = std::get<2>(a);
                        ScalarPoint3f p2 = std::get<2>(b);
                        return p1.z() > p2.z();
                    });
    }

    void sort_intersections_by_distance(std::vector<ShellIntersection> &intersections, ScalarPoint3f p) const {
        //  Sort the intersection points by distance to the ray origin
        std::sort(intersections.begin(), intersections.end(),
                  [p](const ShellIntersection &a, const ShellIntersection &b) {
                        ScalarPoint3f p1 = std::get<2>(a);
                        ScalarPoint3f p2 = std::get<2>(b);

                        return dr::norm(p1 - p) < dr::norm(p2 - p);
                    });
    }

    std::tuple<ScalarFloat, std::vector<CumulativeOTEntry>> compute_cumulative_ot(const uint32_t angle_idx, const std::vector<ShellIntersection> intersections) const {
        static constexpr size_t SpectralSize = dr::size_v<UnpolarizedSpectrum>;
        using ScalarUnpolarized = dr::Array<dr::scalar_t<UnpolarizedSpectrum>, SpectralSize>;

        ScalarUnpolarized cumulative = dr::zeros<ScalarUnpolarized>();
        std::vector<CumulativeOTEntry> cot(intersections.size() * SpectralSize);
        
        //  Midpoint cumulative optical thickness
        ScalarFloat mid_cot = 0.0f;

        // Not enough intersections to compute cumulative OT -> Either at the 
        // edge of the outer shell (which should not happen given the sample origin)
        // or the ray is completely outside the medium (which should also
        // not happen in the precomputation scenario).
        if (intersections.size() < 2)
            return std::make_tuple(0.0f, cot);

        //  The first intersection is always the top shell, 
        //  so we can use it to initialize the CDF
        ShellIntersection first_intersection = intersections[0];
        int32_t r_idx               = std::get<0>(first_intersection);
        int32_t intersection_idx    = std::get<1>(first_intersection);
        ScalarPoint3f p1            = std::get<2>(first_intersection);
        
        for (size_t k = 0; k < SpectralSize; ++k) 
            cot[k] = std::make_tuple(r_idx, intersection_idx, k, p1, 0.0f);

        //  Iterate over the intersections and calculate the optical thickness
        //  for each segment between two intersections.
        MediumInteraction3f mei = dr::zeros<MediumInteraction3f>();
        for (size_t i = 1; i < intersections.size(); ++i) {
            ShellIntersection primary   = intersections[i - 1];
            ShellIntersection secondary = intersections[i];

            ScalarPoint3f p_1 = std::get<2>(primary);
            ScalarPoint3f p_2 = std::get<2>(secondary);

            //  Get the metadata of the SECOND intersection, for which we have
            //  to store the data in the cot LUT
            int32_t r_idx               = std::get<0>(secondary);
            int32_t intersection_idx    = std::get<1>(secondary);

            //  Calculate the length of the segment
            ScalarFloat segment_length = dr::norm(p_2 - p_1);

            //  Calculate the middle point to use as sample point 
            //  to get the medium interaction
            ScalarPoint3f mid_point = (p_1 + p_2) * 0.5f;

            //  Exploit the fact that the medium is spherical and sample a point
            //  has only to be characterized by it's norm towards the origin.
            //  We add a small epsilon to avoid hitting shell boundaries exactly.
            ScalarPoint3f sample_point{0.0f, 0.0f, 2.0f * (dr::abs(dr::norm(mid_point))) - m_medium_radius + 1e-04f};

            //  Set the medium interaction point
            mei.p = sample_point;

            //  Get the scattering coefficients at the sample point
            std::tie(mei.sigma_s, mei.sigma_n, mei.sigma_t) =
                get_scattering_coefficients(mei, true);

            //  If the primary and secondary intersection are on the same shell, 
            //  we use them to compute the midpoint cumulative optical thickness.
            if (r_idx == std::get<0>(primary)) {
                if constexpr (dr::is_jit_v<UnpolarizedSpectrum>)
                    mid_cot = cumulative[0] + ScalarFloat(mei.sigma_t[0][0]) * segment_length * 0.5f;
                else
                    mid_cot = cumulative[0] + ScalarFloat(mei.sigma_t[0]) * segment_length * 0.5f;
            }

            //  Calculate the optical thickness per wavelength channel
            for (size_t k = 0; k < SpectralSize; ++k) {
                if constexpr (dr::is_jit_v<UnpolarizedSpectrum>)
                    cumulative[k] += ScalarFloat(mei.sigma_t[k][0]) * segment_length;
                else
                    cumulative[k] += ScalarFloat(mei.sigma_t[k]) * segment_length;

                cot[i * SpectralSize + k] = std::make_tuple(
                    r_idx, intersection_idx, k, p1, cumulative[k]
                );
            }
        }

        return std::make_tuple(mid_cot, cot);
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
    
        //  Keep track of the maximum number of intersections per angle
        std::vector<uint32_t> number_of_intersections(m_angles.size(), 0.0f);

        //  Keep track of the midpoint cumulative optical thickness for different angles
        std::vector<ScalarFloat> midpoint_cot(m_angles.size(), 0.0f);

        //  The maximum number of intersection is the resolution * the number of angles *
        //  the spectral size (if spectral data is used)
        std::vector<ScalarFloat> opt_thickness_data(m_angles.size() * m_max_intersections * SpectralSize);

        //  Here we assume we start at the top shell
        for (int32_t i = 0; i < (int32_t) m_angles.size(); ++i) {
            //  We assume that the voxel size is the radius of a spherical shell!
            //  The intersection point with the shell from P can then be found by 
            //  using the equation of a circle (x^2 + y^2 = r^2) and of the line
            //  defining the directional through P (y = theta * x + r_max) where 
            //  r_max is defined by the z-coordinate of the max bounding box point.
            ScalarFloat alpha = m_angles[i];

            std::vector<ShellIntersection> intersections = compute_ray_intersections(p, alpha, false);

            //  Add number of intersections per angle
            number_of_intersections[i] = intersections.size();

            //  Compute cumulative OT from the intersection points
            std::tuple<ScalarFloat, std::vector<CumulativeOTEntry>> result = compute_cumulative_ot(i, intersections);
            ScalarFloat mid_cot = std::get<0>(result);
            std::vector<CumulativeOTEntry> cot = std::get<1>(result);

            //  Process the cumulative OT for proper texture indexing
            for (size_t j = 0; j < cot.size(); ++j) {
                CumulativeOTEntry entry     = cot[j];
                int32_t r_idx               = std::get<0>(entry);
                int32_t intersection_idx    = std::get<1>(entry);
                int32_t spectral_idx        = std::get<2>(entry);
                ScalarFloat cot_value       = std::get<4>(entry);

                //  Compute the intersection index (i.e. intersection number of max intersections)
                size_t icts_idx = (intersection_idx == 0) 
                                ? (m_max_intersections / 2 - r_idx - 1) 
                                : (m_max_intersections / 2 + r_idx); 

                //  Compute the flat index based on the intersection index, 
                //  the max number of intersections, the angle index and 
                //  the spectral index.
                size_t index = spectral_idx * (m_angles.size() * m_max_intersections)
                                + i * m_max_intersections
                                + icts_idx;
                
                //  Store the cumulative OT value in the data array
                opt_thickness_data[index] = cot_value;
            }

            //  Save midpoint cumulative optical thickness
            midpoint_cot[i] = mid_cot;

            //Log(Warn, "Angle: %s", m_angles[i]);
            //Log(Warn, "Midpoint COT: %s", mid_cot);

            //  At indices that would be zero in the current angle, we store the 
            //  interpolated midpoint cot value.
            for (size_t j = 1; j < m_max_intersections; ++j) {
                size_t index = i * m_max_intersections + j;

                //  If the current index is zero, we store the midpoint cot value
                if (opt_thickness_data[index] == 0.0f)
                    opt_thickness_data[index] = mid_cot;
            }
        }   

        //  Store the maximum number of intersections per angle
        m_intersections_per_angle = dr::load<IntStorage>(number_of_intersections.data(), number_of_intersections.size());

        //  Store the midpoint cot values
        m_midpoint_cot = dr::load<FloatStorage>(midpoint_cot.data(), midpoint_cot.size());

        size_t shape[] = { m_angles.size(), m_max_intersections, 1 };

        //  Store the cumulative OT as a texture
        m_cot = Texture2f(TensorXf(opt_thickness_data.data(), 3, shape), true,
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
    ref<Volume>     m_sigmat, m_albedo;
    uint32_t        m_angle_samples;
    ScalarFloat     m_angle_step;
    uint32_t        m_max_intersections;
    ScalarFloat     m_medium_radius;

    ScalarFloat     m_texture_cell_angle_step;
    ScalarFloat     m_texture_cell_shell_step;
    ScalarFloat     m_texture_cell_angle_offset;
    ScalarFloat     m_texture_cell_shell_offset;

    mutable FloatStorage m_radii;
    mutable FloatStorage m_angles;
    mutable FloatStorage m_extinctions;
    mutable FloatStorage m_midpoint_cot;
    mutable IntStorage m_intersections_per_angle;

    Float m_max_density;
    mutable Texture2f m_cot;
};

MI_IMPLEMENT_CLASS_VARIANT(PiecewiseSphericalMedium, Medium)
MI_EXPORT_PLUGIN(PiecewiseSphericalMedium, "Piecewise Spherical Medium")
NAMESPACE_END(mitsuba)