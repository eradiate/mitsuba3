#include <mitsuba/render/sensor.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/transform.h>
#include <mitsuba/core/bbox.h>
#include <mitsuba/core/warp.h>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _sensor-fisheye:

Fisheye camera (:monosp:`fisheye`)
----------------------------------

.. pluginparameters::

 * - to_world
   - |transform|
   - Specifies an optional camera-to-world transformation.
     (Default: none (i.e. camera space = world space))
   - |exposed|, |differentiable|, |discontinuous|

 * - fov
   - |float|
   - The field of view in degrees, i.e. the full angle subtended by the image
     circle. Maps to the rim of the circle. Ignored when
     :monosp:`projection_model` is :monosp:`polynomial`, where the calibration
     fully defines the mapping. (Default: :monosp:`180`)

 * - projection_model
   - |string|
   - The radial fisheye mapping :math:`r(\theta)` relating a point's angle
     :math:`\theta` from the optical axis to its normalised image radius
     :math:`r \in [0, 1]`. One of:

     1. :monosp:`equidistant` (:math:`r \propto \theta`): radius linear in angle.
     2. :monosp:`equisolid` (:math:`r \propto \sin(\theta/2)`): equal solid angle
        per unit image area (equal-area) -- recommended for radiometry.
     3. :monosp:`stereographic` (:math:`r \propto \tan(\theta/2)`): conformal
        (preserves local shapes).
     4. :monosp:`orthographic` (:math:`r \propto \sin\theta`): strong edge
        compression; requires :monosp:`fov` :math:`\le 180`.
     5. :monosp:`polynomial`: a calibrated lens whose **normalised** image radius
        follows :math:`\rho(\theta) = A\,\theta + B\,\theta^2` for
        :math:`\theta \in [0, \pi/2]` (see :monosp:`lens_a`, :monosp:`lens_b`).

     (Default: :monosp:`equisolid`)

 * - lens_a
   - |float|
   - Polynomial calibration coefficient :math:`A` (linear term) of the
     normalised projection :math:`\rho(\theta) = A\,\theta + B\,\theta^2`. Must
     be :math:`> 0`. Only used when :monosp:`projection_model` is
     :monosp:`polynomial`.

 * - lens_b
   - |float|
   - Polynomial calibration coefficient :math:`B` (quadratic term) of the
     normalised projection :math:`\rho(\theta) = A\,\theta + B\,\theta^2`. Only
     used when :monosp:`projection_model` is :monosp:`polynomial`.

 * - center_x, center_y
   - |float|
   - Optical centre (principal point) in **normalised film coordinates**
     :math:`[0, 1]` per axis, measured from the top-left corner. (Default:
     :monosp:`0.5, 0.5`, i.e. the film centre)

 * - radius
   - |float|
   - Radius of the image circle (where the normalised image radius reaches 1),
     expressed as a **fraction of the film width**. The pixel scale is isotropic:
     the same fraction of the film width covers the same angle in every
     direction, so the image circle stays circular on non-square films.
     (Default: the disk inscribed in the film, i.e.
     :math:`0.5\,\min(W, H)/W`)

 * - near_clip, far_clip
   - |float|
   - Distance to the near/far clip *spheres*: rays start at
     :monosp:`near_clip` and extend to :monosp:`far_clip` measured radially
     along the ray direction. (Default: :monosp:`near_clip=1e-2`,
     :monosp:`far_clip=1e4`)
   - |exposed|

 * - srf
   - |spectrum|
   - Sensor Response Function that defines the :ref:`spectral sensitivity <explanation_srf_sensor>`
     of the sensor (Default: :monosp:`none`)

 * - x_fov
   - |float|
   - The field of view in degrees (alias exposed via ``traverse``).
   - |exposed|

This plugin implements a circular fisheye camera with an infinitely
small aperture (infinite depth of field). Rays are emitted radially, the optical
centre of the film images the optical axis and the image circle maps to the
field of view, so :monosp:`fov=180` captures the full forward hemisphere.

Image orientation follows the same convention as the :monosp:`perspective`
camera: the :monosp:`look_at` :monosp:`up` vector points to the **top** of the
image and the camera :math:`+x` axis to the **left**. For an upward-looking
fisheye (optical axis = zenith), this reproduces a physically faithful
hemispherical photograph viewed from below.

The camera position and orientation are most easily expressed with a
:monosp:`look_at` transform (point :monosp:`target` where you want the centre
of the image to look):

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
class FisheyeCamera final : public ProjectiveCamera<Float, Spectrum> {
public:
    MI_IMPORT_BASE(ProjectiveCamera, m_to_world, m_needs_sample_3,
                   m_film, m_sampler, m_resolution, m_shutter_open,
                   m_shutter_open_time, m_near_clip, m_far_clip,
                   sample_wavelengths)
    MI_IMPORT_TYPES()

    enum class Projection {
        Equidistant, Equisolid, Stereographic, Orthographic,
        EquisolidFull, Polynomial
    };

    FisheyeCamera(const Properties &props) : Base(props) {
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
        ScalarVector2f size(m_film->size());
        ScalarFloat radius_default = 0.5f * dr::minimum(size.x(), size.y()) / size.x();
        m_center = ScalarPoint2f(props.get<ScalarFloat>("center_x", 0.5f),
                                 props.get<ScalarFloat>("center_y", 0.5f));
        m_circle_radius = props.get<ScalarFloat>("radius", radius_default);
        if (m_circle_radius <= 0.f)
            Throw("The image-circle radius must be > 0 (got %f).", m_circle_radius);
        if (m_center.x() < 0.f || m_center.x() > 1.f ||
            m_center.y() < 0.f || m_center.y() > 1.f)
            Throw("The optical centre must lie within the film, i.e. in "
                  "[0, 1]^2 (got %f, %f).", m_center.x(), m_center.y());

        // Read the calibration only in polynomial mode to avoid Mitsuba
        // warning about unqueried properties for the other models.
        if (m_projection == Projection::Polynomial) {
            ScalarFloat lens_a = props.get<ScalarFloat>("lens_a");
            ScalarFloat lens_b = props.get<ScalarFloat>("lens_b");
            if (lens_a <= 0.f)
                Throw("The polynomial projection requires lens_a > 0 (got %f).",
                      lens_a);
            m_lens_a = lens_a;
            m_lens_b = lens_b;
        }

        if (m_to_world.scalar().has_scale())
            Throw("Scale factors in the camera-to-world transformation are not allowed!");

        m_needs_sample_3 = false;
        update_mapping();
    }

    void traverse(TraversalCallback *cb) override {
        Base::traverse(cb);
        cb->put("x_fov",    m_fov,      ParamFlags::NonDifferentiable);
        cb->put("to_world", m_to_world, ParamFlags::NonDifferentiable);
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
     * `q` is a film coordinate measured from the optical centre in units of
     * the image-circle radius, so |q| = 1 on the rim: |q| is the normalised
     * image radius fed to theta(r), q/|q| the azimuth. Its y axis is weighted
     * by the film aspect ratio, so that a pixel step is the same length on
     * both axes and the image circle stays circular on non-square films.
     */
    void update_mapping() {
        ScalarVector2f size(m_film->size()),
                       crop_size(m_film->crop_size()),
                       crop_offset(m_film->crop_offset());
        ScalarFloat aspect = size.y() / size.x();

        ScalarVector2f axis_scale(1.f / m_circle_radius,
                                  aspect / m_circle_radius);
        m_q_scale  = Vector2f(crop_size / size * axis_scale);
        m_q_offset = Vector2f((ScalarVector2f(m_center) - crop_offset / size)
                              * axis_scale);

        // Precompute trig constants to avoid redundant computation in
        // sample_to_direction (called 3x per pixel in sample_ray_differential)
        m_theta_max = 0.5f * m_fov * dr::Pi<ScalarFloat> / 180.f;
        switch (m_projection) {
            case Projection::Equidistant:
                m_trig_const = dr::sin(m_theta_max);
                break;
            case Projection::Equisolid:
                m_trig_const = dr::sin(0.5f * m_theta_max);
                break;
            case Projection::Stereographic:
                m_trig_const = dr::tan(0.5f * m_theta_max);
                break;
            case Projection::Orthographic:
                m_trig_const = dr::sin(m_theta_max);
                break;
            case Projection::Polynomial:
                m_trig_const = 1.f;  // unused, but initialize for consistency
                break;
            case Projection::EquisolidFull:
                m_trig_const = 1.f;  // unused, but initialize for consistency
                break;
        }

        dr::make_opaque(m_q_scale, m_q_offset, m_fov, m_lens_a, m_lens_b,
                        m_theta_max, m_trig_const);
    }

    Vector3f sample_to_direction(const Point2f &position_sample, Mask &valid) const {
        // Centred, isotropic film coordinates: |q| = 1 on the image circle
        // (see update_mapping() for the sign and aspect conventions).
        Vector2f q = m_q_offset - Vector2f(position_sample) * m_q_scale;
        Float r    = dr::norm(q);
        valid      = r <= 1.f;

        Float theta = 0.f;
        // theta_max and m_trig_const are precomputed in update_mapping() to avoid
        // redundant computation (sample_to_direction is called 3x per pixel in
        // sample_ray_differential).
        Float theta_max = m_theta_max;
        switch (m_projection) {
            case Projection::Equidistant:
                theta = r * theta_max;
                break;
            case Projection::Equisolid:
                theta = 2.f * dr::asin(r * m_trig_const);
                break;
            case Projection::Stereographic:
                theta = 2.f * dr::atan(r * m_trig_const);
                break;
            case Projection::Orthographic:
                theta = dr::asin(dr::minimum(r * m_trig_const, Float(1.f)));
                break;
            case Projection::Polynomial: {
                // rho(theta) = A*theta + B*theta^2 is the normalised image
                // radius, so r = rho. Invert B*theta^2 + A*theta - r = 0 for
                // theta in [0, pi/2] in the numerically stable form
                // 2r / (A + sqrt(A^2 + 4Br)), which is exact for B = 0 and
                // whose denominator is >= A > 0.
                Float disc = dr::maximum(m_lens_a * m_lens_a + 4.f * m_lens_b * r,
                                         Float(0.f));
                theta = 2.f * r / (m_lens_a + dr::sqrt(disc));

                // theta is only defined on [0, pi/2]; a calibration whose horizon
                // sits inside the image circle blacks out the outer ring.
                Float half_pi = 0.5f * dr::Pi<Float>;
                valid &= theta <= half_pi;
                theta = dr::minimum(theta, half_pi);
                break;
            }
            case Projection::EquisolidFull:
                valid = true;
                return warp::square_to_uniform_hemisphere(position_sample);
        }
        auto [sin_t, cos_t] = dr::sincos(theta);

        // Azimuth (cos phi, sin phi) = q / r, guarded against r = 0 (centre).
        Vector2f azim = dr::select(r > 0.f, q / r, Vector2f(0.f, 0.f));

        return Vector3f(sin_t * azim.x(), sin_t * azim.y(), cos_t);
    }

    std::pair<Ray3f, Spectrum> sample_ray(Float time, Float wavelength_sample,
                                          const Point2f &position_sample,
                                          const Point2f & /*aperture_sample*/,
                                          Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::EndpointSampleRay, active);

        auto [wavelengths, wav_weight] =
            sample_wavelengths(dr::zeros<SurfaceInteraction3f>(),
                               wavelength_sample,
                               active);
        Mask valid;
        Ray3f ray;
        ray.time = time;
        ray.wavelengths = wavelengths;
        ray.d = m_to_world.value() * sample_to_direction(position_sample, valid);
        ray.o = m_to_world.value().translation() + ray.d * Float(m_near_clip);
        ray.maxt = Float(m_far_clip - m_near_clip);

        return { ray, wav_weight & valid };
    }

    std::pair<RayDifferential3f, Spectrum>
    sample_ray_differential(Float time, Float wavelength_sample, const Point2f &position_sample,
                            const Point2f & /*aperture_sample*/, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::EndpointSampleRay, active);

        RayDifferential3f ray;
        ray.time = time;
        ray.has_differentials = true;

        auto [wavelengths, wav_weight] =
            sample_wavelengths(dr::zeros<SurfaceInteraction3f>(),
                               wavelength_sample,
                               active);
        ray.wavelengths = wavelengths;

        Mask valid;
        Point2f dx(1.f / m_resolution.x(), 0.f),
                dy(0.f, 1.f / m_resolution.y());
        Mask ignore;
        ray.d   = m_to_world.value() * sample_to_direction(position_sample, valid);
        ray.d_x = m_to_world.value() * sample_to_direction(position_sample + dx, ignore);
        ray.d_y = m_to_world.value() * sample_to_direction(position_sample + dy, ignore);

        ray.o = m_to_world.value().translation() + ray.d * Float(m_near_clip);
        ray.o_x = ray.o_y = ray.o;
        ray.maxt = Float(m_far_clip - m_near_clip);

        return { ray, wav_weight & valid };
    }

    ProjectiveTransform4f projection_transform() const override {
        Throw("FisheyeCamera::projection_transform(): the nonlinear fisheye "
              "mapping cannot be represented by a projective transform; "
              "projective-sampling-based integrators are not supported.");
    }

    ScalarBoundingBox3f bbox() const override {
        ScalarPoint3f p = m_to_world.scalar() * ScalarPoint3f(0.f);
        return ScalarBoundingBox3f(p, p);
    }

    std::string to_string() const override {
        using string::indent;

        std::ostringstream oss;
        oss << "FisheyeCamera[" << std::endl
            << "  x_fov = " << m_fov << "," << std::endl
            << "  near_clip = " << m_near_clip << "," << std::endl
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
private:
    Float m_fov;
    Projection m_projection;
    // Optical centre (film UV) and image-circle radius (fraction of film width)
    ScalarPoint2f m_center;
    ScalarFloat m_circle_radius;
    // Affine film-sample -> centred isotropic coordinate map (see update_mapping)
    Vector2f m_q_scale, m_q_offset;
    // Polynomial projection calibration (used only when m_projection == Polynomial)
    Float m_lens_a = 0.f, m_lens_b = 0.f;
    // Precomputed per-projection trig constants (cached to avoid redundant
    // computation in sample_to_direction's 3x per-pixel calls from sample_ray_differential)
    Float m_theta_max;
    Float m_trig_const;  // sin(theta_max/2), tan(theta_max/2), or sin(theta_max)

    MI_TRAVERSE_CB(Base, m_q_scale, m_q_offset, m_fov, m_lens_a, m_lens_b,
                   m_theta_max, m_trig_const)
};

MI_EXPORT_PLUGIN(FisheyeCamera)
NAMESPACE_END(mitsuba)
