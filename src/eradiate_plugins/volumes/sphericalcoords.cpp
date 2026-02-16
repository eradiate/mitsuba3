#include <mitsuba/core/properties.h>
#include <mitsuba/core/transform.h>
#include <mitsuba/render/volume.h>

NAMESPACE_BEGIN(mitsuba)

enum class FilterType { Nearest, Trilinear };
enum class WrapMode { Repeat, Mirror, Clamp };

/**!
.. _plugin-volume-sphericalcoordsvolume:

Mapping to spherical coordinates (:monosp:`sphericalcoordsvolume`)
------------------------------------------------------------------

.. pluginparameters::

 * - volume
   - |volume|
   - Nested volume plugin whose data is to be mapped to spherical coordinates.
   - —

 * - rmin
   - |float|
   - Radius for the inner limit of the spherical shell, relative to the unit
     sphere. Default: 0
   - —

 * - rmax
   - |float|
   - Radius for the outer limit of the spherical shell, relative to the unit
     sphere. Default: 1
   - —

 * - fillmin
   - |float|
   - Constant value to return for points such that :math:`r < r_\mathrm{min}`.
     Default: 0
   - —

 * - fillmax
   - |float|
   - Constant value to return for points such that :math:`r_\mathrm{max} < r`.
     Default: 0
   - —

 * - to_world
   - |transform|
   - Specifies an optional 4x4 transformation matrix that will remap local
     spherical coordinates from the unit sphere (which covers the [-1, 1]³ cube)
     to world coordinates.
   - —

This plugin addresses volume data in spherical coordinates. In practice, it
maps the texture coordinates of a nested volume plugin to the unit sphere using
the following correspondance:

.. math::

    x \in [0, 1] & \leftrightarrow r \in [r_\mathrm{min}, r_\mathrm{max}] \\
    y \in [0, 1] & \leftrightarrow \theta \in [0, \pi] \\
    z \in [0, 1] & \leftrightarrow \phi \in [-\pi, \pi]

where :math:`r` is the radius of the considered point in the unit sphere.
For angles, the default mathematical convention is used:
:math:`\theta` is the zenith angle with respect to the :math:`+Z` unit vector,
and :math:`\phi` is the azimuth angle with respect to the :math:`(+X, +Z)`
plane.

.. note::

    When using this plugin with a nested ``gridvolume``, the data layout
    remains unchanged (*i.e.* zyxc will be interpreted as φθrc).
*/

// Note: This plugin is primarily designed to be used with a spherical stencil,
// but it should work with other shapes.

template <typename Float, typename Spectrum>
class SphericalCoordsVolume final : public Volume<Float, Spectrum> {
public:
    MI_IMPORT_BASE(Volume, m_to_local, m_bbox)
    MI_IMPORT_TYPES(VolumeGrid)

    using VolumeType = Volume<Float, Spectrum>;

    SphericalCoordsVolume(const Properties &props) : Base(props) {
        m_volume = props.get_volume<VolumeType>("volume", 1.f);

        m_rmin = props.get<ScalarFloat>("rmin", 0.f);
        m_rmax = props.get<ScalarFloat>("rmax", 1.f);
        if (m_rmin > m_rmax)
            Throw("rmin must be lower than rmax!");

        m_fillmin = props.get<ScalarFloat>("fillmin", 0.f);
        m_fillmax = props.get<ScalarFloat>("fillmax", 0.f);

        m_to_local = props.get<ScalarAffineTransform4f>("to_world", ScalarAffineTransform4f()).inverse();
        update_bbox_sphere();
    }

    UnpolarizedSpectrum eval(const Interaction3f &it, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::TextureEvaluate, active);

        Point3f p = m_to_local * it.p;
        Float r = dr::norm(p);

        Point3f p_spherical = Point3f(
            (r - m_rmin) / (m_rmax - m_rmin),
            dr::acos(p.z() / r) * dr::InvPi<ScalarFloat>,
            dr::atan2(p.y(), p.x()) * dr::InvTwoPi<ScalarFloat> + .5f
        );
        Interaction3f it_spherical = it;
        it_spherical.p = p_spherical;

        return dr::select(r < m_rmin,
            m_fillmin,
            dr::select(r > m_rmax,
                m_fillmax,
                m_volume->eval(it_spherical, active)
            )
        );
    }

    Float eval_1(const Interaction3f &it, Mask active = true) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::TextureEvaluate, active);

        Point3f p = m_to_local * it.p;
        Float r = dr::norm(p);

        Point3f p_spherical = Point3f(
            (r - m_rmin) / (m_rmax - m_rmin),
            dr::acos(p.z() / r) * dr::InvPi<ScalarFloat>,
            dr::atan2(p.y(), p.x()) * dr::InvTwoPi<ScalarFloat> + .5f
        );
        Interaction3f it_spherical = it;
        it_spherical.p = p_spherical;

        return dr::select(
            r < m_rmin,
            m_fillmin,
            dr::select(
                r > m_rmax,
                m_fillmax,
                m_volume->eval_1(it_spherical, active)
            )
        );
    }

    ScalarFloat max() const override { return dr::maximum(dr::maximum(m_volume->max(), m_fillmin), m_fillmax); }

    ScalarVector3i resolution() const override { return m_volume->resolution(); };

    const ScalarFloat* data() const override {
        return m_volume->data();
    }

    const DynamicBuffer<Float>* array() const override {
        return m_volume->array();
    }

    std::pair<Float, Float>
    extremum(const DynamicBuffer<Float>* array,
             BoundingBox3f local_bounds) const override {
        // local_bounds is in normalized [0,1]^3 space (r_norm, theta_norm, phi_norm).
        // The nested volume uses the same coordinate convention, so we forward
        // directly after clamping to [0,1].

        // Clamp bounds to valid [0,1]^3 range for the nested volume query
        BoundingBox3f clamped_bounds(
            dr::maximum(local_bounds.min, Point3f(0.f)),
            dr::minimum(local_bounds.max, Point3f(1.f))
        );

        auto [maj, min] = m_volume->extremum(array, clamped_bounds);

        // If the r-range extends below 0 (below rmin), include fillmin
        auto below_rmin = local_bounds.min.x() < 0.f;
        dr::masked(maj, below_rmin) = dr::maximum(maj, Float(m_fillmin));
        dr::masked(min, below_rmin) = dr::minimum(min, Float(m_fillmin));

        // If the r-range extends above 1 (above rmax), include fillmax
        auto above_rmax = local_bounds.max.x() > 1.f;
        dr::masked(maj, above_rmax) = dr::maximum(maj, Float(m_fillmax));
        dr::masked(min, above_rmax) = dr::minimum(min, Float(m_fillmax));

        return { maj, min };
    }

    void traverse(TraversalCallback *cb) override {
        cb->put("volume", m_volume.get(), ParamFlags::NonDifferentiable);
        Base::traverse(cb);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "SphericalCoordsVolume[" << std::endl
            << "  to_local = " << string::indent(m_to_local, 13) << "," << std::endl
            << "  bbox = " << string::indent(m_bbox) << "," << std::endl
            << "  volume = " << string::indent(m_volume) << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS(SphericalCoordsVolume)

private:
    ScalarFloat m_rmin, m_rmax, m_fillmin, m_fillmax;
    ref<VolumeType> m_volume;

    void update_bbox_sphere() {
        ScalarAffineTransform4f to_world = m_to_local.inverse();
        ScalarPoint3f a = to_world * ScalarPoint3f(-1.f, -1.f, -1.f);
        ScalarPoint3f b = to_world * ScalarPoint3f(1.f, 1.f, 1.f);
        m_bbox = ScalarBoundingBox3f(a, b);
    }

    MI_TRAVERSE_CB(Base, m_volume)
};

MI_EXPORT_PLUGIN(SphericalCoordsVolume)

NAMESPACE_END(mitsuba)
