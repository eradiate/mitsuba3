#include <cstring>

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

template <typename Type>
struct AttributeCallback : public TraversalCallback {

    AttributeCallback()
        : found_rmin(false), found_rmax(false), rmin(0.f), rmax(0.f) {}

    void put_object(const std::string &/*name*/, Object */*obj*/, uint32_t) override {};

    void put_parameter_impl(const std::string &name, void *val, uint32_t,
                            const std::type_info &type) override {
        if ("rmin" == name) {
            if (strcmp(type.name(), typeid(Type).name()) == 0)
                rmin = *((Type *) val);
            found_rmin = true;
        } else if ("rmax" == name) {
            if (strcmp(type.name(), typeid(Type).name()) == 0)
                rmax = *((Type *) val);
            found_rmax = true;
        }
    };

    bool found() { return found_rmin && found_rmax; };

    bool found_rmin;
    bool found_rmax;
    Type rmin;
    Type rmax;
};


/* Map the shell index to the lower and higher cot map index. */
template <typename Type>
std::tuple<Type, Type> shell_to_cot_idx(Type shell_index, Type shell_resolution) {
    Type lower_index = (shell_resolution-1) - shell_index;
    Type higher_index = shell_resolution + 2 + shell_index;
    return {lower_index,higher_index};
}

/* Map the cot index to the lower and higher shells of the shell index. */
template <typename Type>
std::tuple<Type, Type> cot_to_shell_idx(Type cot_idx, Type cot_size, Type shell_resolution) {
    Type lower_midpoint_idx = cot_size / 2; // lower midpoint index 
    Type inner_idx, outer_idx;

    inner_idx = dr::select(cot_idx < lower_midpoint_idx,
                           shell_resolution - (cot_idx + 1),
                           cot_idx - (shell_resolution + 2));
    outer_idx = dr::select(cot_idx < lower_midpoint_idx, 
                           shell_resolution - cot_idx,
                           cot_idx - (shell_resolution + 1));

    dr::masked(outer_idx, outer_idx == shell_resolution) = -1;
    return {inner_idx, outer_idx};
}

template <typename Float, typename Spectrum>
class PiecewiseSphericalMedium final : public Medium<Float, Spectrum> {
public:
    MI_IMPORT_BASE(Medium, m_is_homogeneous, m_has_spectral_extinction,
                   m_phase_function)
    MI_IMPORT_TYPES(Scene, Sampler, Texture, Volume)

    static constexpr ScalarFloat PI_HALF     = dr::Pi<Float> / 2.0f;
    static constexpr ScalarFloat UPPER_BOUND = PI_HALF - 10*dr::Epsilon<Float>;

    static constexpr size_t SpectralSize = dr::size_v<UnpolarizedSpectrum>;
    
    using ScalarIndex           = uint32_t;
    using ScalarSize            = uint32_t;
    using ScalarUnpolarized = dr::Array<dr::scalar_t<UnpolarizedSpectrum>, SpectralSize>;
    using FloatStorage          = DynamicBuffer<Float>;
    using IntStorage            = DynamicBuffer<uint32_t>;
    using Index                 = dr::uint32_array_t<Float>;

    PiecewiseSphericalMedium(const Properties &props) : Base(props) {
        m_is_homogeneous = false;
        m_albedo         = props.volume<Volume>("albedo", 0.75f);
        m_sigmat         = props.volume<Volume>("sigma_t", 1.f);
        m_has_spectral_extinction =
            props.get<bool>("has_spectral_extinction", true);
        m_max_density    = dr::opaque<Float>(m_sigmat->max());
        m_medium_radius = m_sigmat->bbox().extents().z() * 0.5;

        // angular resolution
        ScalarIndex angle_samples  = props.get<int32_t>("angular_samples", 2);
        m_angle_size  = angle_samples + 1; // account for nadir direction
        m_angle_step  = UPPER_BOUND / ScalarFloat(angle_samples);

        // retrieve rmin and rmax from the volume
        AttributeCallback<ScalarFloat> cb;
        const_cast<Volume*>(m_sigmat.get())->traverse((TraversalCallback *) &cb);
        if (!cb.found()) {
            if constexpr (!dr::is_jit_v<Float>)
                Throw("Invalid attribute requested: rmax or rmax.");
        }
        m_rmin = cb.rmin;
        m_rmax = cb.rmax;
        
        // shell resolution
        if (m_rmax < 1.0f - dr::Epsilon<Float>) {
            m_additional_shells++;
            m_max_shell = true;
        }
        
        if (m_rmin > 0.f + dr::Epsilon<Float>){
            m_additional_shells++;
            m_min_shell = true;
        }
        Log(Debug, "rmin: %f, rmax: %f, m_min_shell: %d, m_max_shell: %d", m_rmin, m_rmax, m_min_shell, m_max_shell);
        
        // Size of COT dimension given by the number of shells in the volume and 
        // the min and max radiux boudaries. Additional two points for the midpoint
        // which will allow for interpolation. 
        const ScalarIndex base_resolution = m_sigmat->resolution().x();
        m_shell_resolution = base_resolution + m_additional_shells;
        m_shell_size = (m_rmax - m_rmin) / ScalarFloat(base_resolution);
        m_cot_size = (m_sigmat->resolution().x() + m_additional_shells + 1) * 2;
        
        Log(Debug, "medium_radius: %f, max_intersections: %f, m_additional_shells: %d", m_medium_radius, m_cot_size, m_additional_shells);

        m_texture_cell_angle_step   = 1.0f / m_angle_size;
        m_texture_cell_shell_step   = 1.0f / m_cot_size;
        m_texture_cell_angle_offset = 0.5f * m_texture_cell_angle_step;
        m_texture_cell_shell_offset = 0.5f * m_texture_cell_shell_step;

        parameters_changed();
    }


    std::tuple<typename Medium<Float, Spectrum>::MediumInteraction3f, Float, Float>
    sample_interaction_real(const Ray3f &ray, const SurfaceInteraction3f &si,
                            Float sample, UInt32 channel,
                            Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::MediumSample, active);

        // Initial intersection with the medium
        Point3f center = m_sigmat->bbox().center();
        Float radius = m_sigmat->bbox().extents().z();
        BoundingSphere3f bsphere(center, radius);
        auto [aabb_its, mint, maxt] = bsphere.ray_intersect(ray);
        aabb_its &= (dr::isfinite(mint) || dr::isfinite(maxt));
        active &= aabb_its;
        dr::masked(mint, !active) = 0.f;
        dr::masked(maxt, !active) = dr::Infinity<Float>;

        // Point where ray enters the sphere 
        Point3f entry_point;
        dr::masked(entry_point, active) = ray.o + ray.d*mint;

        mint = dr::maximum(0.f, mint);
        maxt = dr::minimum(si.t, dr::minimum(ray.maxt, maxt));

        Mask escaped = !active;

        // Initialize basic medium interaction fields
        MediumInteraction3f mei = dr::zeros<MediumInteraction3f>();
        mei.wi                  = -ray.d;
        mei.sh_frame            = Frame3f(mei.wi);
        mei.time                = ray.time;
        mei.wavelengths         = ray.wavelengths;
        mei.mint                = mint;
        mei.t                   = mint;
        mei.medium              = this;

        // get angle for interpolation
        Vector3f nadir = dr::normalize(center - entry_point);
        Float cos_alpha = dr::dot(nadir, ray.d);
        Float sin_alpha = dr::sqrt(1.f - cos_alpha*cos_alpha);
        Float alpha = dr::acos(cos_alpha);
        alpha = dr::clip(alpha, 0.f, UPPER_BOUND);

        // get midpoint in relative coordinates
        Point2f o(sin_alpha, cos_alpha);
        Point2f d(0.f, -1.f);
        Vector2f midpoint = o + d*cos_alpha;
        Float midpoint_dist_to_center = sin_alpha;
        auto [mid_shell_idx, lower_radius, higher_radius] = shell_info<Index, Float>(midpoint_dist_to_center);
        auto [l_mid_idx, h_mid_idx] = shell_to_cot_idx(mid_shell_idx, Index(m_shell_resolution));
        l_mid_idx++;
        h_mid_idx--;
        // intersections[l_mid_idx] = {true, cos_alpha - dr::Epsilon<ScalarFloat>};
        // intersections[h_mid_idx] = {true, cos_alpha + dr::Epsilon<ScalarFloat>};
        
        Log(Debug, "entry_point: %f, center: %f, nadir: %f", entry_point, center, nadir );

        if (dr::any_or<false>(cos_alpha < 0.f)) {
            Log(Warn, "Ray traversing in negative direction; ray.d: %f, cos_alpha: %f", ray.d, cos_alpha);
        }

        // Float cum_opt_thick = dr::zeros<Float>();
        // Float sampled_t     = dr::Infinity<Float>;
        // Float tr            = dr::zeros<Float>();
        // Float pdf           = dr::zeros<Float>();
        
        Float neg_log_ksi = -dr::log(1.0f - sample);

        // convert angle to texture coordinate
        Float angle_coord = m_texture_cell_angle_offset + alpha * dr::rcp(UPPER_BOUND);
        Log(Debug, "log_sample: %f, alpha: %f, angle_coord: %f", neg_log_ksi, alpha, angle_coord );

        Index start_idx   = 0;
        Index end_idx     = m_cot_size -1;
        Index index = dr::binary_search<Index>(
            start_idx, end_idx,
            [&](Index idx) DRJIT_INLINE_LAMBDA {
                // Texture Sampling
                Float shell_coord = m_texture_cell_shell_offset + idx * m_texture_cell_shell_step;
                Point2f p(shell_coord, angle_coord);
                
                UnpolarizedSpectrum cot;
                m_cot.template eval<Float>(p, cot.data());
                Float cot_value = extract_channel(cot, channel);

                Log(Debug, "binary search; idx: %f, shell_coord: %f cot: %f", idx, shell_coord, cot );
                return cot_value < neg_log_ksi;
            });

        Log(Debug, "neg_log_ksi: %f, index: %f", neg_log_ksi, index );
        auto [inner_s_idx, outer_s_idx] = cot_to_shell_idx(index, Index(m_cot_size), Index(m_shell_resolution));

        return {dr::zeros<MediumInteraction3f>(), 0.f, 0.f};
    }

    std::tuple<Float, Float, Mask>
    eval_transmittance_pdf_real(const Ray3f &/*ray*/,
                                const SurfaceInteraction3f &/*si*/, UInt32 /*channel*/,
                                Mask /*active*/) const override {
        return {0.f, 0.f, true};
    }

    void traverse(TraversalCallback *callback) override {
        callback->put_object("albedo", m_albedo.get(),
                             +ParamFlags::Differentiable);
        callback->put_object("sigma_t", m_sigmat.get(),
                             +ParamFlags::Differentiable);
        callback->put_parameter("cum_opt_thickness", m_cot.tensor(),
                                +ParamFlags::NonDifferentiable);
        // callback->put_parameter("angles", m_angles,
        //                         +ParamFlags::NonDifferentiable);
        Base::traverse(callback);   
    }

    void parameters_changed(const std::vector<std::string> & /*keys*/ = {}) override {
        
        
        using ScalarIntersection = typename std::tuple<bool, ScalarFloat>;

        std::vector<ScalarUnpolarized> extinctions;
        precompute_extinctions(extinctions);

        Log(Debug, "=================");
        Log(Debug, "Intersections");
        Log(Debug, "=================");

        std::vector<ScalarUnpolarized> cot_data(m_angle_size * m_cot_size);

        auto write_to_cot_buffer = [&](ScalarUnpolarized data, ScalarIndex angle_idx, ScalarIndex cot_idx) {
            ScalarIndex index = angle_idx * m_cot_size + cot_idx;
            cot_data[index] = data;
        };

        for (ScalarIndex angle_idx = 0; angle_idx < m_angle_size; ++angle_idx) {
            std::vector<ScalarIntersection> intersections(m_cot_size, {false, 0.f});
            // std::vector<ScalarUnpolarized> cots(m_cot_size, -1.f);
            ScalarFloat alpha = m_angle_step * angle_idx;
            
            // 01 - Initialize ray geometry
            auto [sin_alpha, cos_alpha] = dr::sincos(alpha);
            // !Important! Nadir at angle alpha = 0
            ScalarVector2f d( sin_alpha, cos_alpha );
            ScalarVector2f o( 0.f, -1.f );
            Log(Debug, "=========================");
            Log(Debug, "alpha : %f, d: %d", alpha, d);
            Log(Debug, "_________________________");

            // 02 - Calculate midpoint location and distance
            // The midpoint is the point along the ray closest to the center.
            // This will serve as starting point in shell intersections
            ScalarVector2f midpoint = o + d*cos_alpha;
            ScalarFloat midpoint_dist_to_center = sin_alpha;
            Log(Debug, "midpoint: %f, midpoint_dist:%f", midpoint, midpoint_dist_to_center);

            auto [mid_shell_idx, lower_radius, higher_radius] = shell_info<ScalarIndex, ScalarFloat>(midpoint_dist_to_center);
            Log(Debug, "mid_shell_idx: %d, lower_r: %f, higher_r:%f ", mid_shell_idx, lower_radius, higher_radius);

            auto [l_mid_idx, h_mid_idx] = shell_to_cot_idx(mid_shell_idx, m_shell_resolution);
            l_mid_idx++;
            h_mid_idx--;
            intersections[l_mid_idx] = {true, cos_alpha - dr::Epsilon<ScalarFloat>};
            intersections[h_mid_idx] = {true, cos_alpha + dr::Epsilon<ScalarFloat>};
            Log(Debug, "l_mid_idx: %d, h_mid_idx", l_mid_idx, h_mid_idx);
            
            // 03 - Calculate intersections with shell boundaries
            for( ScalarIndex s_idx = mid_shell_idx; s_idx < m_shell_resolution; ++s_idx) {
                std::tie(lower_radius, higher_radius) = shell_bounds<ScalarIndex, ScalarFloat>(s_idx);
                auto [l_idx, h_idx] = shell_to_cot_idx(s_idx, m_shell_resolution);
                auto [valid, x0, x1] = math::solve_quadratic(
                    dr::squared_norm(d),
                    2.f * dr::dot(o, d),
                    dr::squared_norm(o) - dr::square(higher_radius)
                );

                Log(Debug, "s_idx: %d, l_idx: %d, h_idx: %d, x0: %f, x1: %f", s_idx, l_idx, h_idx, x0, x1);

                // should be valid since we are starting from the midpoint
                if(valid) {
                    intersections[l_idx] = {true, x0};
                    intersections[h_idx] = {true, x1};
                } else {
                    Log(Warn, "Invalid intersection at s_idx: %d", s_idx);
                }
            }

            Log(Debug, "______________");
            // 04 - Calculate Cumulative optical thickness from intersections
            ScalarFloat prev_distance = 0.f;
            ScalarUnpolarized cot = 0.f;
            // cots[0] = cot; // first entry is always 0.
            write_to_cot_buffer(ScalarUnpolarized(0.f), angle_idx, 0);

            for (ScalarIndex cot_idx = 1; cot_idx < m_cot_size; ++cot_idx ){
                auto [valid, distance]  = intersections[cot_idx];
                
                if (!valid) {
                    Log(Debug,"No intersections at cot_idx %d, continue", cot_idx);
                    continue;
                }

                auto [inner_s_idx, outer_s_idx] = cot_to_shell_idx(cot_idx, m_cot_size, m_shell_resolution);
                
                // retrieve the shell index from the cot_idx, will depend on 
                // whether cot_idx is the midpoint, before the midpoint or past
                // the midpoint.
                ScalarIndex shell_idx = 0;
                if(cot_idx == l_mid_idx || cot_idx == h_mid_idx) {
                    shell_idx = outer_s_idx;
                } else if (cot_idx <= (m_cot_size / 2) ) {
                    shell_idx = outer_s_idx;
                } else {
                    shell_idx = inner_s_idx;
                }
                Log(Debug, "valid: %f, distance: %f", valid, distance);

                Log(Debug, "cot_idx: %f, inner_s_idx: %f, outer_s_idx: %f, shell_idx: %f", cot_idx, inner_s_idx, outer_s_idx, shell_idx);

                // accumulate the optical thickness 
                cot += (distance - prev_distance)*extinctions[shell_idx];

                //  Compute the flat index based on the intersection index, 
                //  the max number of intersections, the angle index and 
                //  the spectral index.
                
                write_to_cot_buffer(cot, angle_idx, cot_idx);
                // cots[cot_idx] = cot;

                Log(Debug, "cot: %f, dist: %f, sigma_t: %f", cot, (distance - prev_distance), extinctions[shell_idx]);
                prev_distance = distance;
            }

            Log(Debug, "______________");
            // 05 - Fill the invalid entries with the midpoint cot
            ScalarUnpolarized midpoint_cot = cot_data[angle_idx * m_cot_size + h_mid_idx];
            for (ScalarIndex cot_idx = h_mid_idx; cot_idx > l_mid_idx; --cot_idx) {
                write_to_cot_buffer(midpoint_cot, angle_idx, cot_idx);
                // cots[cot_idx] = midpoint_cot;
                Log(Debug, "cot_idx: %d, midpoint_cot: %f", cot_idx, midpoint_cot);
            }

            for (ScalarIndex cot_idx = 0; cot_idx < m_cot_size; ++cot_idx) {
                Log(Debug, "cot_idx: %d, cot: %f", cot_idx, cot_data[angle_idx*m_cot_size+cot_idx]);
            }
        }

        //  Store the cumulative OT as a texture
        size_t shape[] = { m_angle_size, m_cot_size, SpectralSize };

        m_cot = Texture2f(TensorXf(cot_data.data(), 3, shape), true,
            true, dr::FilterMode::Linear, dr::WrapMode::Clamp);
    }

    
    void precompute_extinctions(std::vector<ScalarUnpolarized>& extinctions) const {
        const ScalarVector3f center = m_sigmat->bbox().center();
        const ScalarIndex base_resolution = m_sigmat->resolution().x();
        const ScalarIndex resolution = base_resolution + m_additional_shells;
        
        extinctions.assign(resolution, ScalarUnpolarized(0.f));

        Log(Debug, "base_res %f, res %f, additional_res : %f, shell_size %f", base_resolution, resolution, m_additional_shells, m_shell_size);

        MediumInteraction3f mei = dr::zeros<MediumInteraction3f>();
        ScalarFloat previous_radius = 0.f;

        Log(Debug, "retrieve extinction");
        for (ScalarIndex r_idx = 0 ; r_idx < resolution; ++r_idx) {
            ScalarPoint3f p = center;
            if (m_min_shell && r_idx == 0){
                p[2] += m_medium_radius * m_rmin*0.5;
                mei.p = p;
                previous_radius += m_rmin;
            } else if (m_max_shell && r_idx == (resolution-1)) {
                p[2] += m_medium_radius * (1.f+m_rmax)*0.5;
                mei.p = p;
            } else {
                p[2] += m_medium_radius * (previous_radius+0.5*m_shell_size);
                mei.p = p;
                previous_radius += m_shell_size;
            }

            //  Get the scattering coefficients at the sample point
            std::tie(mei.sigma_s, mei.sigma_n, mei.sigma_t) =
                get_scattering_coefficients(mei, true);
            Log(Debug, "idx : %d, sigma_t: %f, mei.p: %f", r_idx, mei.sigma_t, mei.p);
            //  Store the extinction coefficient for the current shell
            extinctions[r_idx] = dr::slice(mei.sigma_t, 0);
        }

        Log(Debug, "store extinction per shell");
        for (ScalarIndex r_idx = 0 ; r_idx < resolution; ++r_idx) {
            Log(Debug, "idx : %d, sigma_t: %f", r_idx, extinctions[r_idx]);
        }
    }

    /* 
    * Return shell info from the point position. Assume p to be in local
    * volume coordinate.
    * Returns (shell_idx, lower_bound, higher_bound)
    * */
    template<typename _Index, typename _Float>
    std::tuple<_Index, _Float, _Float>
    shell_info(_Float dist) const {
        _Float rmin = _Float(m_rmin);
        _Float rmax = _Float(m_rmax);
        _Float shell_size = _Float(m_shell_size);
        _Index shell_resolution = _Index(m_shell_resolution);
        _Index additional_shells = _Index(m_additional_shells);
        dr::mask_t<_Index> min_shell = dr::mask_t<_Index>(m_min_shell);

        _Index i = dr::floor2int<_Index>((dist - rmin) / shell_size);
        i = dr::clip(i, 0, (shell_resolution-additional_shells)-1);
        _Float lower_bound = rmin + shell_size*i;
        _Float higher_bound = lower_bound + shell_size;
        
        dr::masked(i, min_shell) += 1;
        
        dr::masked(lower_bound, dist < rmin) = 0.f;
        dr::masked(higher_bound, dist < rmin) = rmin;
        dr::masked(i, dist < rmin) = 0;

        dr::masked(lower_bound, dist > rmax) = rmax;
        dr::masked(higher_bound, dist > rmax) = 1.f;
        dr::masked(i, dist > rmax) = shell_resolution -1;

        return {i, lower_bound, higher_bound};
    }

    /* 
    * Return shell bounds for a given shell index. 
    * Returns (shell_idx, lower_bound, higher_bound)
    * */
   template<typename _Index, typename _Float>
    std::tuple<_Float, _Float>
    shell_bounds(_Index shell_index) const {
        dr::mask_t<_Index> min_shell = dr::mask_t<_Index>(m_min_shell);
        dr::mask_t<_Index> max_shell = dr::mask_t<_Index>(m_max_shell);
        _Index shell_resolution = _Index(m_shell_resolution);
        _Float shell_size = _Float(m_shell_size);
        _Float rmin = _Float(m_rmin);
        _Float rmax = _Float(m_rmax);

        // Account for the rmin layer
        _Index i = dr::select(min_shell, shell_index - 1 , shell_index);
        _Float lower_bound  = rmin + i * shell_size;
        _Float higher_bound = rmin + (i + 1) * shell_size;

        dr::masked(lower_bound, shell_index == 0 && min_shell) = 0.f;
        dr::masked(higher_bound, shell_index == 0 && min_shell) = rmin;

        dr::masked(lower_bound, (shell_resolution - 1) == 0 && max_shell) = rmax;
        dr::masked(higher_bound, (shell_resolution - 1) == 0 && max_shell) = 1.f;
        
        return {lower_bound, higher_bound};
    }

    bool is_almost_zero(ScalarFloat value) const {
        return dr::abs(value) < 1e-08f;
    }

 
    UnpolarizedSpectrum get_majorant(const MediumInteraction3f &/*mi*/,
                                     Mask /*active*/) const override {
        return 0.f;
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
            dr::masked(result, channel == 1u) = value[1];
            dr::masked(result, channel == 2u) = value[2];
        } else {
            DRJIT_MARK_USED(channel);
        }

        return result;
    }

    MI_DECLARE_CLASS()

private:
    ref<Volume>     m_sigmat, m_albedo;
    Float m_max_density;
    
    ScalarSize      m_angle_size;
    ScalarSize      m_cot_size;
    ScalarSize      m_shell_resolution;
    ScalarFloat     m_shell_size;

    ScalarFloat     m_angle_step;
    ScalarFloat     m_rmin, m_rmax;
    bool            m_min_shell = false, m_max_shell = false;
    ScalarIndex     m_additional_shells = 0.f;
    ScalarFloat     m_medium_radius;

    ScalarFloat     m_texture_cell_angle_step;
    ScalarFloat     m_texture_cell_shell_step;
    ScalarFloat     m_texture_cell_angle_offset;
    ScalarFloat     m_texture_cell_shell_offset;

    mutable Texture2f m_cot;
};

MI_IMPLEMENT_CLASS_VARIANT(PiecewiseSphericalMedium, Medium)
MI_EXPORT_PLUGIN(PiecewiseSphericalMedium, "Piecewise Spherical Medium")
NAMESPACE_END(mitsuba)