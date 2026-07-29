#include <mitsuba/core/properties.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/eradiate/extremum.h>
#include <mitsuba/render/volume.h>

NAMESPACE_BEGIN(mitsuba)

/**!
.. _extremum-extremum_repeat:

Extremum repeat structure (:monosp:`extremum_repeat`)
-----------------------------------------------------

.. pluginparameters::

 * - (Nested plugin)
   - |extremum|
   - The inner extremum structure to tile.

 * - lattice
   - |vector|
   - Axis-aligned lattice period along x, y, z. Default: the extents of the
     inner structure's domain.

 * - aabb_min, aabb_max
   - |point|
   - Tiled region in world space. Default: infinite.

This plugin periodically tiles an inner extremum structure over a translation
lattice. It is assembled — typically by a :monosp:`repeat` medium — from the
inner medium's *built* structure; \ref ExtremumStructure::build() is never
called on it.

Only next_segment() is implemented: it folds the query point into the inner
structure's tile and delegates a single-segment lookup, clipped to the tile
edge. Traversal itself uses ExtremumStructure's generic default, which steps
through next_segment() one segment at a time — this keeps traversal cheap
whether this structure is used standalone or aggregated as a child of e.g.
extremum_overlap.
*/

template <typename Float, typename Spectrum>
class ExtremumRepeat final : public ExtremumStructure<Float, Spectrum> {
public:
    MI_IMPORT_BASE(ExtremumStructure, m_bbox)
    MI_IMPORT_TYPES(Volume, ExtremumStructure)

    ExtremumRepeat(const Properties &props) : Base(props) {
        for (auto &prop : props.objects()) {
            if (auto *child = prop.try_get<ExtremumStructure>()) {
                if (m_inner)
                    Throw("extremum_repeat accepts a single nested extremum "
                          "structure");
                m_inner = child;
            }
        }
        if (!m_inner)
            Throw("extremum_repeat requires a nested extremum structure");

        ScalarBoundingBox3f tile = m_inner->bbox();
        if (!tile.valid() || !dr::all(dr::isfinite(tile.min) && dr::isfinite(tile.max)))
            Throw("extremum_repeat: the inner structure must have a finite "
                  "domain (one tile), got %s", tile);
        m_origin = tile.min;

        m_lattice = props.get<ScalarVector3f>("lattice", tile.extents());
        m_lattice_rcp = 1.f / m_lattice;

        if (props.has_property("aabb_min") && props.has_property("aabb_max")) {
            m_bbox = ScalarBoundingBox3f(props.get<ScalarPoint3f>("aabb_min"),
                                         props.get<ScalarPoint3f>("aabb_max"));
        } else {
            m_bbox = ScalarBoundingBox3f(
                ScalarPoint3f(-dr::Infinity<ScalarFloat>),
                ScalarPoint3f(dr::Infinity<ScalarFloat>));
        }
    }

    ExtremumSegment next_segment(const Ray3f &ray, Float t,
                                 Mask active) const override {
        auto [hit, d0, d1] = m_bbox.ray_intersect(ray);

        ExtremumSegment segment(t, dr::Infinity<Float>, Vector2f(0.f));

        Float eps = dr::maximum(dr::abs(t), 1.f) * math::RayEpsilon<Float>;
        Float tq  = t + eps;

        active &= hit;
        Mask before = active && (tq < d0);
        Mask inside = active && !before && (tq < d1);

        dr::masked(segment.maxt, before) = d0;

        // early exit
        if (dr::any_or<false>(!inside))
            return segment;

        auto [shift, t_edge] = tile_at(ray, t);
        Ray3f shifted = ray;
        shifted.o -= shift;

        ExtremumSegment s = m_inner->next_segment(shifted, t, inside);
        dr::masked(segment.value, inside) = s.value;
        dr::masked(segment.maxt, inside) = dr::clip(dr::minimum(s.maxt, t_edge), tq, d1);

        return segment;
    }

    std::tuple<Float, Float> eval_1(const Interaction3f &it,
                                    Mask active) const override {
        Vector3f c = (it.p - m_origin) * m_lattice_rcp;
        Interaction3f it_folded = it;
        it_folded.p = it.p - m_lattice * dr::floor(c);
        return m_inner->eval_1(it_folded, active);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "ExtremumRepeat[" << std::endl
            << "  inner = " << string::indent(m_inner) << "," << std::endl
            << "  lattice = " << m_lattice << "," << std::endl
            << "  bbox = " << m_bbox << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS(ExtremumRepeat)

private:
    /// World translation of the lattice cell containing p(t + eps), and the
    /// exit distance of the ray from that cell (exact in t)
    std::pair<Vector3f, Float> tile_at(const Ray3f &ray, const Float &t) const {
        Float eps = dr::maximum(dr::abs(t), 1.f) * math::RayEpsilon<Float>;

        // Ray in lattice coordinates: same t parameterization
        Vector3f o = (ray.o - m_origin) * m_lattice_rcp;
        Vector3f d = ray.d * m_lattice_rcp;
        Vector3f k = dr::floor(dr::fmadd(d, t + eps, o));

        Float t_edge = dr::Infinity<Float>;
        for (size_t i = 0; i < 3; ++i) {
            Float b  = k[i] + dr::select(d[i] > 0.f, 1.f, 0.f);
            Float ti = dr::select(d[i] != 0.f, (b - o[i]) / d[i],
                                  dr::Infinity<Float>);
            dr::masked(ti, ti <= t + eps) = dr::Infinity<Float>;
            t_edge = dr::minimum(t_edge, ti);
        }

        return { k * m_lattice, t_edge };
    }

private:
    ref<ExtremumStructure> m_inner;
    ScalarVector3f m_lattice = ScalarVector3f(1.f);
    ScalarVector3f m_lattice_rcp = ScalarVector3f(1.f);
    ScalarPoint3f m_origin = 0.f;
};

MI_EXPORT_PLUGIN(ExtremumRepeat)
NAMESPACE_END(mitsuba)
