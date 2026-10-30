#include <mitsuba/core/bbox.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/string.h>
#include <mitsuba/core/transform.h>
#include <mitsuba/render/fwd.h>
#include <mitsuba/render/interaction.h>
#include <mitsuba/render/shape.h>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _shape-acube:

Analytical cube (:monosp:`acube`)
---------------------------------

.. pluginparameters::

 * - flip_normals
   - |bool|
   - Is the cube inverted, i.e. should the normal vectors be flipped?
     (Default: |false|, i.e. the normals point outside)

 * - to_world
   - |transform|
   - Specifies an optional linear object-to-world transformation.
     (Default: none (i.e. object space = world space))

This shape plugin describes the cube :math:`[-1, 1]^3`, intersected analytically
with a single slab test. It is a drop-in replacement for the
:ref:`cube <shape-cube>` mesh in CPU variants, except as an area emitter
(position sampling is not implemented). The scene acceleration structure holds
a single primitive instead of twelve triangles, which makes intersections
cheaper; the gain is largest when the cube bounds a participating medium
enclosing other geometry.

Only CPU (scalar and LLVM) variants are supported, and the shape is not
differentiable with respect to its parameters.

.. tabs::
    .. code-tab:: python

        'cloud': {
            'type': 'acube',
            'to_world': mi.ScalarTransform4f().scale([2, 2, 1]),
            'bsdf': {'type': 'null'},
            'interior': {'type': 'ref', 'id': 'medium'},
        }

 */

template <typename Float, typename Spectrum>
class AnalyticalCube final : public Shape<Float, Spectrum> {
public:
    MI_IMPORT_BASE(Shape, m_to_world, m_is_instance, initialize,
                   get_children_string)
    MI_IMPORT_TYPES()

    using typename Base::ScalarIndex;

    AnalyticalCube(const Properties &props) : Base(props) {
        if constexpr (dr::is_cuda_v<Float> || dr::is_metal_v<Float>)
            Throw("acube: only CPU variants are supported.");
        m_flip_normals = props.get<bool>("flip_normals", false);
        initialize();
    }

    ScalarBoundingBox3f bbox() const override {
        ScalarBoundingBox3f bbox;
        ScalarAffineTransform4f to_world = m_to_world.scalar();
        for (int i = 0; i < 8; ++i)
            bbox.expand(to_world * ScalarPoint3f(i & 1 ? 1.f : -1.f,
                                                 i & 2 ? 1.f : -1.f,
                                                 i & 4 ? 1.f : -1.f));
        return bbox;
    }

    /* A clip entirely inside the cube does not contain any of the cube
      surface and is excluded. This avoids placing the primitive in every
      kd-tree leaf inside the cube.*/
    ScalarBoundingBox3f bbox(ScalarIndex /*index*/,
                             const ScalarBoundingBox3f &clip) const override {
        ScalarAffineTransform4f to_object = m_to_world.scalar().inverse();
        bool inside = true;
        // Margin keeps nodes whose faces lie on the cube surface
        for (size_t i = 0; i < 8; ++i)
            inside &= dr::all(dr::abs(to_object * clip.corner(i)) < 1.f - 1e-4f);
        if (inside)
            return ScalarBoundingBox3f();

        ScalarBoundingBox3f result = bbox();
        result.clip(clip);
        return result;
    }
    using Base::bbox;

    Float surface_area() const override {
        const AffineTransform4f &to_world = m_to_world.value();
        Vector3f ex = to_world * Vector3f(2.f, 0.f, 0.f),
                 ey = to_world * Vector3f(0.f, 2.f, 0.f),
                 ez = to_world * Vector3f(0.f, 0.f, 2.f);
        return 2.f * (dr::norm(dr::cross(ey, ez)) + dr::norm(dr::cross(ez, ex)) +
                      dr::norm(dr::cross(ex, ey)));
    }

    // =============================================================
    //! @{ \name Ray tracing routines
    // =============================================================

    template <typename FloatP, typename Ray3fP>
    std::tuple<dr::mask_t<FloatP>, FloatP, Point<FloatP, 2>,
               dr::uint32_array_t<FloatP>, dr::uint32_array_t<FloatP>>
    ray_intersect_preliminary_impl(const Ray3fP &ray_,
                                   ScalarIndex /*prim_index*/,
                                   dr::mask_t<FloatP> active) const {
        auto [hit, t] = intersect<FloatP>(ray_, active);
        return { hit, dr::select(hit, t, dr::Infinity<FloatP>),
                 dr::zeros<Point<FloatP, 2>>(), ((uint32_t) -1), 0 };
    }

    template <typename FloatP, typename Ray3fP>
    dr::mask_t<FloatP> ray_test_impl(const Ray3fP &ray_,
                                     ScalarIndex /*prim_index*/,
                                     dr::mask_t<FloatP> active) const {
        return intersect<FloatP>(ray_, active).first;
    }

    MI_SHAPE_DEFINE_RAY_INTERSECT_METHODS()

    SurfaceInteraction3f compute_surface_interaction(const Ray3f &ray,
                                                     const PreliminaryIntersection3f &pi,
                                                     uint32_t /*ray_flags*/,
                                                     uint32_t recursion_depth,
                                                     Mask active) const override {
        MI_MASK_ARGUMENT(active);

        if (!m_is_instance && recursion_depth > 0)
            return dr::zeros<SurfaceInteraction3f>();

        const AffineTransform4f &to_world = m_to_world.value();
        Point3f local = to_world.inverse() * ray(pi.t);

        // The hit face is the axis where the local point is furthest out.
        Vector3f a = dr::abs(local);
        Mask fx = a.x() >= a.y() && a.x() >= a.z(),
             fy = !fx && a.y() >= a.z();
        Vector3f e_n = dr::select(fx, Vector3f(1.f, 0.f, 0.f),
                       dr::select(fy, Vector3f(0.f, 1.f, 0.f), Vector3f(0.f, 0.f, 1.f)));
        Vector3f e_u = dr::select(fx, Vector3f(0.f, 1.f, 0.f),
                       dr::select(fy, Vector3f(0.f, 0.f, 1.f), Vector3f(1.f, 0.f, 0.f)));
        Vector3f e_v = dr::cross(e_n, e_u);
        Float s  = dr::sign(dr::dot(local, e_n)),
              lu = dr::dot(local, e_u),
              lv = dr::dot(local, e_v);

        SurfaceInteraction3f si = dr::zeros<SurfaceInteraction3f>();
        si.t = dr::select(active, pi.t, dr::Infinity<Float>);
        // Re-project onto the face to improve accuracy
        si.p = to_world * Point3f(s * e_n + lu * e_u + lv * e_v);
        si.n = dr::normalize(to_world * Normal3f((m_flip_normals ? -s : s) * e_n));
        si.sh_frame.n = si.n;
        /* v is mirrored on negative faces so that dp_du x dp_dv points outward,
           also with flipped normals, as on meshes */
        si.dp_du = to_world * (2.f * e_u);
        si.dp_dv = to_world * (2.f * s * e_v);
        si.uv    = Point2f(dr::fmadd(lu, .5f, .5f), dr::fmadd(s * lv, .5f, .5f));
        si.dn_du = si.dn_dv = dr::zeros<Vector3f>();
        si.shape    = this;
        si.instance = nullptr;

        return si;
    }

    //! @}
    // =============================================================

    bool has_flipped_normals() const override { return m_flip_normals; }

    std::string to_string() const override {
        using string::indent;
        std::ostringstream oss;
        oss << "AnalyticalCube[" << std::endl
            << "  to_world = " << indent(m_to_world, 13) << "," << std::endl
            << "  flip_normals = " << m_flip_normals << "," << std::endl
            << "  surface_area = " << surface_area() << "," << std::endl
            << "  " << indent(get_children_string()) << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS(AnalyticalCube)

private:
    /// Slab test against [-1, 1]^3 in object space. Returns the entry distance
    /// for rays starting outside and the exit distance for rays starting inside.
    template <typename FloatP, typename Ray3fP>
    std::pair<dr::mask_t<FloatP>, FloatP>
    intersect(const Ray3fP &ray_, dr::mask_t<FloatP> active) const {
        using Point3fP = Point<FloatP, 3>;

        AffineTransform<Point<FloatP, 4>> to_object;
        if constexpr (!dr::is_jit_v<FloatP>)
            to_object = m_to_world.scalar().inverse();
        else
            to_object = m_to_world.value().inverse();

        Ray3fP ray = to_object * ray_;
        BoundingBox<Point3fP> unit(Point3fP(-1.f), Point3fP(1.f));
        auto [hit, mint, maxt] = unit.ray_intersect(ray);

        FloatP t = dr::select(mint >= 0.f, mint, maxt);
        active = active && hit && t >= 0.f && t <= ray.maxt;
        return { active, t };
    }

    bool m_flip_normals;

    MI_TRAVERSE_CB(Base)
};

MI_EXPORT_PLUGIN(AnalyticalCube)
NAMESPACE_END(mitsuba)
