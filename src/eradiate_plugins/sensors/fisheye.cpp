#include <mitsuba/render/sensor.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/transform.h>
#include <mitsuba/core/bbox.h>
#include <mitsuba/core/warp.h>
#include <cmath>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _plugin-sensor-fisheye:

Fisheye camera (:monosp:`fisheye`)
----------------------------------

.. pluginparameters::

 * - to_world
   - |transform|
   - Sensor-to-world transformation matrix.
   - |exposed|

 * - fov
   - |float|
   - The field of view in degrees, i.e. the full angle subtended by the image
     circle. For :monosp:`polynomial` it instead declares the angular domain the
     calibration was fitted over. (Default: :monosp:`180`)
   - —

 * - projection_model
   - |string|
   - The radial mapping from the angle to the optical axis to the normalized
     image radius, see below. One of :monosp:`equidistant`, :monosp:`equisolid`,
     :monosp:`stereographic`, :monosp:`orthographic` or :monosp:`polynomial`.
     (Default: :monosp:`equisolid`)
   - —

 * - lens_coefficients
   - |string|
   - Coefficients of the :monosp:`polynomial` projection, in ascending order, as
     a comma- or whitespace-separated list.
   - —

 * - center_x, center_y
   - |float|
   - Optical centre (principal point) in normalized film coordinates
     :math:`[0, 1]` per axis, measured from the top-left corner. (Default:
     :monosp:`0.5, 0.5`, i.e. the film centre)
   - —

 * - radius
   - |float|
   - Radius of the image circle, expressed as a fraction of the film width.
     (Default: the disk inscribed in the film, i.e.
     :math:`0.5\,\min(W, H)/W`)
   - —

 * - near_clip, far_clip
   - |float|
   - Distance to the near/far clip *spheres*, measured along the ray direction.
     (Default: :monosp:`near_clip=1e-2`, :monosp:`far_clip=1e4`)
   - |exposed|

 * - srf
   - |spectrum|
   - Sensor Response Function that defines the :ref:`spectral sensitivity <explanation_srf_sensor>`
     of the sensor (Default: :monosp:`none`)
   - —

This plugin implements a circular fisheye camera with an infinitely small
aperture (infinite depth of field). The optical centre of the film images the
optical axis and the image circle maps to the field of view, so
:monosp:`fov=180` captures the forward hemisphere and wider settings reach below
the horizon. Samples falling outside the image circle are discarded.
Image orientation follows the :monosp:`perspective` camera convention.

.. note::

   No projection can flatten a hemisphere faithfully, so each model gives up one
   thing to keep another. All map the angle :math:`\theta` from the optical axis
   to the normalized image radius :math:`r`, with
   :math:`\theta_\mathrm{max} = \mathrm{fov}/2`:

   1. :monosp:`equidistant`, :math:`r = \theta / \theta_\mathrm{max}`: radius
      proportional to angle.
   2. :monosp:`equisolid`,
      :math:`r = \sin(\theta/2) / \sin(\theta_\mathrm{max}/2)`: equal solid angle
      per unit image area.
   3. :monosp:`stereographic`,
      :math:`r = \tan(\theta/2) / \tan(\theta_\mathrm{max}/2)`: conformal, i.e.
      shapes stay undistorted locally, at the cost of stretching the periphery.
   4. :monosp:`orthographic`, :math:`r = \sin\theta / \sin\theta_\mathrm{max}`:
      parallel projection of the hemisphere onto the image plane, which
      compresses the periphery, :monosp:`fov` :math:`\le 180`.
   5. :monosp:`polynomial`, :math:`\rho(\theta) = c_0\,\theta + c_1\,\theta^2 +
      \ldots` with up to 16 coefficients usually obtained from callibration.

   The first four are rescaled to reach :math:`r = 1` exactly at
   :math:`\theta_\mathrm{max}`; the polynomial is used as calibrated instead. So
   :math:`\rho` must increase over :math:`[0, \theta_\mathrm{max}]` (the kernel
   inverts it at construction), and :math:`\rho(\theta_\mathrm{max})` decides
   what the rim shows — below 1 the image stops short of the image circle, above
   1 it fills before :monosp:`fov` is reached.

.. tabs::
    .. code-tab:: python

        'type': 'fisheye',
        'fov': 180,
        'projection_model': 'equisolid',
        'to_world': mi.ScalarTransform4f().look_at(
            origin=[0, 0, 0.5],
            target=[0, 0, 1],   # look straight up
            up=[0, 1, 0]
        ),
        'film_id': {
            'type': '<film_type>',
            # ...
        },
        'sampler_id': {
            'type': '<sampler_type>',
            # ...
        }

 */

template <typename Float, typename Spectrum>
class FisheyeCamera final : public Sensor<Float, Spectrum> {
public:
    MI_IMPORT_BASE(Sensor, m_to_world, m_needs_sample_3,
                   m_film, m_sampler, m_resolution, m_shutter_open,
                   m_shutter_open_time, sample_wavelengths)
    MI_IMPORT_TYPES()

    enum class Projection {
        Equidistant, Equisolid, Stereographic, Orthographic,
        EquisolidFull, Polynomial
    };

    /// Generous sanity limit, anything above will probably be due to a bad input
    static constexpr size_t MaxPolynomialOrder = 16;
    /// Entries in the theta(rho) lookup table
    static constexpr size_t TableSize = 512;
    /// Samples used to check monotonicity and to measure the table's error
    static constexpr size_t CheckSamples = 4096;
    /// Angle a pixel subtends on a 4000 px frame; the yardstick for the above
    static constexpr float PixelScale = 8e-4f;

    using FloatStorage = DynamicBuffer<Float>;

    /**
     * Build the theta(rho) lookup table the sampling routine reads.
     *
     * The calibration is a polynomial of any order, but the camera needs the
     * inverse, which has no closed form beyond a quadratic. Rather than solve
     * it per sample, invert it once here onto a uniform rho grid: at render
     * time a lookup is then two gathers and a lerp, whatever the order.
     *
     * Returns the largest angular error the table introduces, measured against
     * the exact calibration -- see the caller, which reports it.
     */
    ScalarFloat build_table(const std::vector<ScalarFloat> &coeffs,
                            ScalarFloat theta_max) {
        // Horner's method
        auto rho = [&](ScalarFloat theta) {
            ScalarFloat r = 0.f;
            for (size_t i = coeffs.size(); i-- > 0;)
                r = (r + coeffs[i]) * theta;
            return r;
        };

        // rho must increase across the domain, or theta(rho) is ambiguous.
        ScalarFloat prev = 0.f;
        for (size_t i = 1; i <= CheckSamples; ++i) {
            ScalarFloat cur = rho(theta_max * ScalarFloat(i) / ScalarFloat(CheckSamples));
            if (cur <= prev)
                Throw("The polynomial projection requires rho(theta) to increase "
                      "over [0, fov/2] = [0, %f] rad, but it stops doing so near "
                      "theta = %f rad.", theta_max,
                      theta_max * ScalarFloat(i) / ScalarFloat(CheckSamples));
            prev = cur;
        }
        m_rho_max = prev;

        // Invert onto a uniform rho grid. Bisection here is a construction-time cost
        auto invert = [&](ScalarFloat target) {
            ScalarFloat lo = 0.f, hi = theta_max;
            for (size_t i = 0; i < 64; ++i) {
                ScalarFloat mid = 0.5f * (lo + hi);
                if (rho(mid) < target) lo = mid; else hi = mid;
            }
            return 0.5f * (lo + hi);
        };

        std::vector<ScalarFloat> table(TableSize);
        for (size_t i = 0; i < TableSize; ++i)
            table[i] = invert(m_rho_max * ScalarFloat(i) / ScalarFloat(TableSize - 1));

        m_theta_table = dr::load<FloatStorage>(table.data(), TableSize);
        m_table_scale = ScalarFloat(TableSize - 1) / m_rho_max;

        // Measure what the interpolation actually costs, so the caller can say
        // so rather than leaving the user to discover it
        ScalarFloat worst = 0.f;
        for (size_t i = 0; i <= CheckSamples; ++i) {
            ScalarFloat exact = theta_max * ScalarFloat(i) / ScalarFloat(CheckSamples);
            ScalarFloat idx   = rho(exact) * m_table_scale;
            size_t k          = dr::minimum((size_t) idx, TableSize - 2);
            ScalarFloat t     = idx - ScalarFloat(k);
            ScalarFloat approx = table[k] * (1.f - t) + table[k + 1] * t;
            worst = dr::maximum(worst, dr::abs(approx - exact));
        }
        return worst;
    }

    /// Read the 'lens_coefficients' calibration: a comma- or whitespace-separated list
    std::vector<ScalarFloat> read_coefficients(const Properties &props) {
        std::vector<ScalarFloat> coeffs;

        if (!props.has_property("lens_coefficients"))
            return coeffs;

        for (const std::string &token : string::tokenize(
                 props.get<std::string_view>("lens_coefficients"), " ,")) {
            try {
                coeffs.push_back((ScalarFloat) std::stod(token));
            } catch (const std::logic_error &) {
                Throw("Could not parse \"%s\" in 'lens_coefficients' as a "
                      "number.", token);
            }
        }

        return coeffs;
    }

    FisheyeCamera(const Properties &props) : Base(props) {
        m_near_clip = props.get<ScalarFloat>("near_clip", 1e-2f);
        m_far_clip  = props.get<ScalarFloat>("far_clip", 1e4f);
        if (m_near_clip <= 0.f)
            Throw("The 'near_clip' parameter must be greater than zero!");
        if (m_near_clip >= m_far_clip)
            Throw("The 'near_clip' parameter must be smaller than 'far_clip'.");

        ScalarFloat fov_deg = props.get<ScalarFloat>("fov", 180.f);
        if (fov_deg <= 0.f || fov_deg >= 360.f)
            Throw("The field of view must be in the range (0, 360) (got %f).",
                  fov_deg);
        m_fov = fov_deg;

        std::string proj =
            string::to_lower(props.get<std::string>("projection_model", "equisolid"));
        if (proj == "equidistant")
            m_projection = Projection::Equidistant;
        else if (proj == "equisolid")
            m_projection = Projection::Equisolid;
        else if (proj == "stereographic")
            m_projection = Projection::Stereographic;
        else if (proj == "orthographic")
            m_projection = Projection::Orthographic;
        // Undocumented developer-only mode, similar to hdistant
        else if (proj == "equisolid_full")
            m_projection = Projection::EquisolidFull;
        else if (proj == "polynomial")
            m_projection = Projection::Polynomial;
        else
            Throw("Unsupported projection \"%s\"! Must be one of: equidistant, "
                  "equisolid, stereographic, orthographic, polynomial.", proj);

        if (m_projection == Projection::Orthographic && fov_deg > 180.f)
            Throw("The orthographic projection requires fov <= 180 (got %f).", fov_deg);

        if (m_projection == Projection::EquisolidFull)
            Log(Warn, "The \"equisolid_full\" projection is a fixed 180 deg "
                      "square-to-hemisphere map that ignores fov, centre, "
                      "radius and crop settings.");

        // Optical centre in film UV and image-circle radius as a fraction of
        // the film width. Defaults: film centre / inscribed disk.
        m_center = ScalarPoint2f(props.get<ScalarFloat>("center_x", 0.5f),
                                 props.get<ScalarFloat>("center_y", 0.5f));
        m_radius_from_film = !props.has_property("radius");
        m_circle_radius = m_radius_from_film
                              ? inscribed_radius()
                              : props.get<ScalarFloat>("radius");
        if (m_circle_radius <= 0.f)
            Throw("The image-circle radius must be > 0 (got %f).", m_circle_radius);
        if (m_center.x() < 0.f || m_center.x() > 1.f ||
            m_center.y() < 0.f || m_center.y() > 1.f)
            Throw("The optical centre must lie within the film, i.e. in "
                  "[0, 1]^2 (got %f, %f).", m_center.x(), m_center.y());

        if (m_projection == Projection::Polynomial) {
            m_lens_coeffs = read_coefficients(props);

            if (m_lens_coeffs.empty())
                Throw("The polynomial projection requires at least one "
                      "coefficient in 'lens_coefficients'.");
            if (m_lens_coeffs.size() > MaxPolynomialOrder)
                Throw("'lens_coefficients' holds %zu coefficients, more than the "
                      "supported maximum of %zu.", m_lens_coeffs.size(),
                      MaxPolynomialOrder);
            for (ScalarFloat value : m_lens_coeffs)
                if (!std::isfinite(value))
                    Throw("'lens_coefficients' contains a non-finite value "
                          "(%f).", value);

            if (m_lens_coeffs[0] <= 0.f)
                Throw("The polynomial projection requires a positive linear "
                      "coefficient (got %f).", m_lens_coeffs[0]);

            ScalarFloat theta_max = 0.5f * dr::deg_to_rad(fov_deg);
            ScalarFloat table_error = build_table(m_lens_coeffs, theta_max);

            if (table_error > PixelScale)
                Throw("This fisheye lens calibration cannot be inverted "
                      "reliably: doing so introduces up to %f rad of error, "
                      "more than %.0fx a pixel on a 4000 px frame. rho(theta) "
                      "is monotonic but nearly stalls, which makes theta(rho) "
                      "ill-conditioned for any method. Check the coefficients "
                      "and 'fov'.", table_error, table_error / PixelScale);

            if (table_error > 0.125f * PixelScale)
                Log(Warn, "This fisheye lens calibration is poorly conditioned: "
                          "inverting it introduces up to %f rad of error "
                          "(%.1f%% of a pixel on a 4000 px frame).",
                    table_error, 100.f * table_error / PixelScale);
        }

        if (m_to_world.scalar().has_scale())
            Throw("Scale factors in the camera-to-world transformation are not allowed!");

        m_needs_sample_3 = false;
        update_mapping();
    }

    void traverse(TraversalCallback *cb) override {
        Base::traverse(cb);
        cb->put("near_clip", m_near_clip, ParamFlags::NonDifferentiable);
        cb->put("far_clip",  m_far_clip,  ParamFlags::NonDifferentiable);
        cb->put("to_world",  m_to_world,  ParamFlags::NonDifferentiable);
    }

    void parameters_changed(const std::vector<std::string> &keys) override {
        Base::parameters_changed(keys);
        if (keys.empty() || string::contains(keys, "to_world")) {
            if (m_to_world.scalar().has_scale())
                Throw("Scale factors in the camera-to-world transformation are not allowed!");
            m_to_world = m_to_world.value().update();
        }

        update_mapping();
    }

    /**
     * Precompute the film-sample -> image-circle coordinate map.
     *
     * `q` is how far a sample lies from the optical centre, in image-circle
     * radii: |q| is the normalized image radius fed to theta(r), q/|q| the
     * azimuth, and |q| > 1 is outside the circle. The y axis carries the film
     * aspect ratio, without which a vertical pixel step would count for more
     * than a horizontal one and the circle would come out elliptical.
     */
    void update_mapping() {
        ScalarVector2f size(m_film->size()),
                       crop_size(m_film->crop_size()),
                       crop_offset(m_film->crop_offset());
        ScalarFloat aspect = size.y() / size.x();

        // The default image circle is the disk inscribed in the film, so it must
        // follow a film resize rather than keeping the size seen at construction
        if (m_radius_from_film)
            m_circle_radius = inscribed_radius();

        ScalarVector2f axis_scale(1.f / m_circle_radius,
                                  aspect / m_circle_radius);
        m_q_scale  = Vector2f(crop_size / size * axis_scale);
        m_q_offset = Vector2f((ScalarVector2f(m_center) - crop_offset / size)
                              * axis_scale);

        // Precompute the per-projection constants read by sample_to_direction.
        m_theta_max = 0.5f * dr::deg_to_rad(m_fov);
        switch (m_projection) {
            case Projection::Equisolid:
                m_trig_const = dr::sin(0.5f * m_theta_max);
                break;
            case Projection::Stereographic:
                m_trig_const = dr::tan(0.5f * m_theta_max);
                break;
            case Projection::Orthographic:
                m_trig_const = dr::sin(m_theta_max);
                break;
            // Equidistant reads theta_max directly, the other two read neither
            case Projection::Equidistant:
            case Projection::Polynomial:
            case Projection::EquisolidFull:
                m_trig_const = 1.f;
                break;
        }

        dr::make_opaque(m_q_scale, m_q_offset);
    }

    Vector3f sample_to_direction(const Point2f &film_sample, Mask &valid) const {
        Vector2f q  = m_q_offset - Vector2f(film_sample) * m_q_scale;
        Float r_raw = dr::norm(q);
        valid       = r_raw <= 1.f;

        // Clamp so out-of-circle samples yield the rim direction, not NaN
        Float r     = dr::minimum(r_raw, 1.f);

        Float theta = 0.f;
        switch (m_projection) {
            case Projection::Equidistant:
                theta = r * m_theta_max;
                break;
            case Projection::Equisolid:
                theta = 2.f * dr::asin(r * m_trig_const);
                break;
            case Projection::Stereographic:
                theta = 2.f * dr::atan(r * m_trig_const);
                break;
            case Projection::Orthographic:
                theta = dr::asin(r * m_trig_const);
                break;
            case Projection::Polynomial: {
                valid &= r <= m_rho_max;

                Float idx  = r * m_table_scale;
                UInt32 i   = dr::minimum(dr::floor2int<UInt32>(idx),
                                         uint32_t(TableSize - 2));
                Float t    = idx - Float(i);
                theta      = dr::lerp(dr::gather<Float>(m_theta_table, i),
                                      dr::gather<Float>(m_theta_table, i + 1u), t);
                break;
            }
            case Projection::EquisolidFull:
                valid = true;
                return warp::square_to_uniform_hemisphere(film_sample);
        }
        auto [sin_t, cos_t] = dr::sincos(theta);

        // Azimuth (cos phi, sin phi) = q / r_raw, guarded against r = 0 (centre).
        // This must use the unclamped radius: q has length r_raw, so dividing by
        // the clamped r would return a direction that is not unit length.
        Vector2f azim = dr::select(r_raw > 0.f, q / r_raw, Vector2f(0.f, 0.f));

        return Vector3f(sin_t * azim.x(), sin_t * azim.y(), cos_t);
    }

    std::pair<Ray3f, Spectrum> sample_ray(Float time, Float wavelength_sample,
                                          const Point2f &film_sample,
                                          const Point2f & /*aperture_sample*/,
                                          Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::EndpointSampleRay, active);

        // Sample spectrum
        auto [wavelengths, wav_weight] =
            sample_wavelengths(dr::zeros<SurfaceInteraction3f>(),
                               wavelength_sample,
                               active);
        Mask valid;
        Ray3f ray;
        ray.time = time;
        ray.wavelengths = wavelengths;

        // Sample ray direction and position the origin on the near clip sphere
        ray.d = m_to_world.value() * sample_to_direction(film_sample, valid);
        ray.o = m_to_world.value().translation() + ray.d * Float(m_near_clip);
        ray.maxt = Float(m_far_clip - m_near_clip);

        return { ray, wav_weight & valid };
    }

    std::pair<RayDifferential3f, Spectrum>
    sample_ray_differential(Float time, Float wavelength_sample, const Point2f &film_sample,
                            const Point2f & /*aperture_sample*/, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::EndpointSampleRay, active);

        RayDifferential3f ray;
        ray.time = time;
        ray.has_differentials = true;

        // Sample spectrum
        auto [wavelengths, wav_weight] =
            sample_wavelengths(dr::zeros<SurfaceInteraction3f>(),
                               wavelength_sample,
                               active);
        ray.wavelengths = wavelengths;

        Mask valid, ignore;
        Point2f dx(1.f / m_resolution.x(), 0.f),
                dy(0.f, 1.f / m_resolution.y());
        ray.d   = m_to_world.value() * sample_to_direction(film_sample, valid);
        ray.d_x = m_to_world.value() * sample_to_direction(film_sample + dx, ignore);
        ray.d_y = m_to_world.value() * sample_to_direction(film_sample + dy, ignore);

        ray.o = m_to_world.value().translation() + ray.d * Float(m_near_clip);
        ray.o_x = ray.o_y = ray.o;
        ray.maxt = Float(m_far_clip - m_near_clip);

        return { ray, wav_weight & valid };
    }

    ScalarBoundingBox3f bbox() const override {
        ScalarPoint3f p = m_to_world.scalar() * ScalarPoint3f(0.f);
        return ScalarBoundingBox3f(p, p);
    }

    std::string to_string() const override {
        using string::indent;

        std::ostringstream oss;
        oss << "FisheyeCamera[" << std::endl
            << "  projection_model = " << projection_name() << "," << std::endl
            << "  fov = " << m_fov << "," << std::endl
            << "  center = " << m_center << "," << std::endl
            << "  radius = " << m_circle_radius << "," << std::endl;

        if (m_projection == Projection::Polynomial) {
            oss << "  lens_coefficients = [";
            for (size_t i = 0; i < m_lens_coeffs.size(); ++i)
                oss << (i ? ", " : "") << m_lens_coeffs[i];
            oss << "]," << std::endl;
        }

        oss << "  near_clip = " << m_near_clip << "," << std::endl
            << "  far_clip = " << m_far_clip << "," << std::endl
            << "  film = " << indent(m_film) << "," << std::endl
            << "  sampler = " << indent(m_sampler) << "," << std::endl
            << "  resolution = " << m_resolution << "," << std::endl
            << "  shutter_open = " << m_shutter_open << "," << std::endl
            << "  shutter_open_time = " << m_shutter_open_time << "," << std::endl
            << "  to_world = " << indent(m_to_world, 13) << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS(FisheyeCamera)

protected:
    ScalarFloat inscribed_radius() const {
        ScalarVector2f size(m_film->size());
        return 0.5f * dr::minimum(size.x(), size.y()) / size.x();
    }

    const char *projection_name() const {
        switch (m_projection) {
            case Projection::Equidistant:   return "equidistant";
            case Projection::Equisolid:     return "equisolid";
            case Projection::Stereographic: return "stereographic";
            case Projection::Orthographic:  return "orthographic";
            case Projection::EquisolidFull: return "equisolid_full";
            case Projection::Polynomial:    return "polynomial";
        }
        return "unknown";
    }

    ScalarFloat m_fov;
    Projection m_projection;
    ScalarFloat m_near_clip, m_far_clip;
    ScalarPoint2f m_center;
    ScalarFloat m_circle_radius;
    bool m_radius_from_film;
    Vector2f m_q_scale, m_q_offset;
    std::vector<ScalarFloat> m_lens_coeffs;
    FloatStorage m_theta_table;
    ScalarFloat m_table_scale = 1.f;
    ScalarFloat m_rho_max = 1.f;
    ScalarFloat m_theta_max;
    ScalarFloat m_trig_const;

    MI_TRAVERSE_CB(Base, m_q_scale, m_q_offset, m_theta_table)
};

MI_EXPORT_PLUGIN(FisheyeCamera)
NAMESPACE_END(mitsuba)
