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

 * - lattice_0, lattice_1, lattice_2
   - |vector|
   - Lattice translation vectors. Default: the axis-aligned extents of the
     inner structure's domain.

 * - aabb_min, aabb_max
   - |point|
   - Tiled region in world space. Default: infinite.

This plugin periodically tiles an inner extremum structure over a translation
lattice. It is assembled — typically by a :monosp:`repeat` medium — from the
inner medium's *built* structure; \ref ExtremumStructure::build() is never
called on it.

Traversal loops over the lattice cells intersected by the ray, delegating
each whole cell span to the inner structure's native traversal with a
lattice-translated ray. Translations preserve ray parameterization and
directions, so the tracking callback contract is untouched; the state's ray
remains the world ray.
*/

template <typename Float, typename Spectrum>
class ExtremumRepeat final : public ExtremumStructure<Float, Spectrum> {
public:
    MI_IMPORT_BASE(ExtremumStructure, m_bbox)
    MI_IMPORT_TYPES(Volume, ExtremumStructure)

    using TrackingStateType    = TrackingState<Float, Spectrum>;
    using TrackingFunctionType = TrackingFunction<Float, Spectrum>;

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

        ScalarVector3f l[3];
        for (size_t i = 0; i < 3; ++i) {
            ScalarVector3f def(0.f);
            def[i] = tile.extents()[i];
            l[i] = props.get<ScalarVector3f>("lattice_" + std::to_string(i),
                                             def);
        }
        // Lattice vectors as matrix columns
        m_lattice = ScalarMatrix3f(l[0].x(), l[1].x(), l[2].x(),
                                   l[0].y(), l[1].y(), l[2].y(),
                                   l[0].z(), l[1].z(), l[2].z());
        m_lattice_inv = dr::inverse(m_lattice);

        if (props.has_property("aabb_min") && props.has_property("aabb_max")) {
            m_bbox = ScalarBoundingBox3f(props.get<ScalarPoint3f>("aabb_min"),
                                         props.get<ScalarPoint3f>("aabb_max"));
        } else {
            m_bbox = ScalarBoundingBox3f(
                ScalarPoint3f(-dr::Infinity<ScalarFloat>),
                ScalarPoint3f(dr::Infinity<ScalarFloat>));
        }
    }

    TrackingStateType traverse_extremum(
        const Ray3f &ray,
        Float mint,
        Float maxt,
        UInt32 channel,
        TrackingStateType state,
        TrackingFunctionType* func,
        Mask active
    ) const override {
        struct LoopState {
            TrackingStateType state;
            Mask active;
            Float t;

            DRJIT_STRUCT(LoopState, state, active, t)
        } ls = {
            state,
            active && (maxt > mint),
            mint
        };

        dr::tie(ls) = dr::while_loop(
            dr::make_tuple(ls),
            [](const LoopState &ls) { return ls.active; },
            [this, func, ray, maxt, channel](LoopState &ls) {

            auto [shift, t_edge] = tile_at(ray, ls.t);
            Float t_end = dr::minimum(t_edge, maxt);

            // Delegate the whole tile span with the lattice-translated ray;
            // the state's ray stays the world ray so that interaction
            // positions remain world-space.
            Ray3f shifted = ray;
            shifted.o -= shift;
            ls.state = m_inner->traverse_extremum(
                shifted, ls.t, t_end, channel, ls.state, func, ls.active);

            // A valid interaction is a sampled real scatter — stop
            ls.active &= !ls.state.mei.is_valid();
            ls.t = t_end;
            ls.active &= ls.t < maxt;
        },
        "Repeat Traversal");

        return ls.state;
    }

    ExtremumSegment next_segment(const Ray3f &ray, Float t,
                                 Mask active) const override {
        auto [hit, d0, d1] = m_bbox.ray_intersect(ray);

        Float eps = dr::maximum(dr::abs(t), 1.f) * math::RayEpsilon<Float>;
        Float tq  = t + eps;

        Mask before = hit && (tq < d0);
        Mask inside = hit && !before && (tq < d1);

        auto [shift, t_edge] = tile_at(ray, t);
        Ray3f shifted = ray;
        shifted.o -= shift;

        ExtremumSegment s = m_inner->next_segment(shifted, t, active && inside);

        // Clamp to the tile edge and the domain so segments compose
        Float maxt = dr::select(
            inside,
            dr::maximum(dr::minimum(dr::minimum(s.maxt, t_edge), d1), tq),
            dr::select(before, d0, dr::Infinity<Float>));

        return ExtremumSegment(t, maxt,
                               dr::select(inside, s.value, Vector2f(0.f)));
    }

    std::tuple<Float, Float> eval_1(const Interaction3f &it,
                                    Mask active) const override {
        Vector3f c = Matrix3f(m_lattice_inv) * (it.p - m_origin);
        Interaction3f it_folded = it;
        it_folded.p = it.p - Matrix3f(m_lattice) * dr::floor(c);
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
        Vector3f o = Matrix3f(m_lattice_inv) * (ray.o - m_origin);
        Vector3f d = Matrix3f(m_lattice_inv) * ray.d;
        Vector3f k = dr::floor(dr::fmadd(d, t + eps, o));

        Float t_edge = dr::Infinity<Float>;
        for (size_t i = 0; i < 3; ++i) {
            Float b  = k[i] + dr::select(d[i] > 0.f, 1.f, 0.f);
            Float ti = dr::select(d[i] != 0.f, (b - o[i]) / d[i],
                                  dr::Infinity<Float>);
            dr::masked(ti, ti <= t + eps) = dr::Infinity<Float>;
            t_edge = dr::minimum(t_edge, ti);
        }

        return { Matrix3f(m_lattice) * k, t_edge };
    }

private:
    ref<ExtremumStructure> m_inner;
    ScalarMatrix3f m_lattice = dr::identity<ScalarMatrix3f>();
    ScalarMatrix3f m_lattice_inv = dr::identity<ScalarMatrix3f>();
    ScalarPoint3f m_origin = 0.f;
};

MI_EXPORT_PLUGIN(ExtremumRepeat)
NAMESPACE_END(mitsuba)
