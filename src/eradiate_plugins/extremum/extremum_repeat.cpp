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

The DDA traversal state is the inner structure's, set up on a copy of the ray
shifted into the current lattice cell. The cell is not stored: its exit
distance is recomputed from the state's ray, and the inner traversal is
restarted in the next cell once the cursor reaches it. Where the lattice is
wider than the inner domain, each stretch of a cell outside that domain is
returned as a single zero-valued segment.

Nesting an :monosp:`extremum_repeat` inside another is not supported.
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
        m_tile_extents = tile.extents() * m_lattice_rcp;
        m_to_lattice = ScalarAffineTransform4f::scale(m_lattice_rcp) *
                       ScalarAffineTransform4f::translate(-m_origin) *
                       m_inner->dda_to_local().inverse();

        if (props.has_property("aabb_min") && props.has_property("aabb_max")) {
            m_bbox = ScalarBoundingBox3f(props.get<ScalarPoint3f>("aabb_min"),
                                         props.get<ScalarPoint3f>("aabb_max"));
        } else {
            m_bbox = ScalarBoundingBox3f(
                ScalarPoint3f(-dr::Infinity<ScalarFloat>),
                ScalarPoint3f(dr::Infinity<ScalarFloat>));
        }
    }

    void build(const Volume * /*volume*/) override {}

    DDAState dda_init(const Ray3f &ray, Float mint, Float maxt,
                      Mask active) const override {
        auto [hit, d0, d1] = m_bbox.ray_intersect(ray);

        Float start = dr::select(hit, dr::maximum(mint, d0), mint);
        Float end   = dr::select(hit && active, dr::minimum(maxt, d1), start);

        Vector3f cell = dr::floor(
            (dr::fmadd(ray.d, start, ray.o) - m_origin) * m_lattice_rcp);
        Ray3f shifted = ray;
        shifted.o -= m_lattice * cell;

        DDAState state = m_inner->dda_init(shifted, start, end, active);
        state.mint = start;
        state.maxt = end;
        return state;
    }

    std::pair<ExtremumSegment, DDAState>
    dda_next(const DDAState &state, Mask active) const override {
        DDAState current = state;
        Point3f q   = m_to_lattice * current.o;
        Vector3f dq = m_to_lattice * current.d;
        Vector3f t_exit = cell_exit(q, dq);
        Float t_edge    = dr::min(t_exit);

        Mask stale = active && current.mint >= t_edge;
        if (dr::any_or<true>(stale)) {
            Vector3f tile_step =
                dr::select(t_exit <= t_edge, dr::sign(dq), 0.f);
            Ray3f shifted(m_origin + m_lattice * (q - tile_step),
                          m_lattice * dq);
            DDAState fresh =
                m_inner->dda_init(shifted, current.mint, current.maxt, stale);
            fresh.mint = current.mint;
            fresh.maxt = current.maxt;
            dr::masked(current, stale) = fresh;

            Point3f q_fresh = m_to_lattice * fresh.o;
            dr::masked(q, stale)      = q_fresh;
            dr::masked(t_edge, stale) = dr::min(cell_exit(q_fresh, dq));
        }

        auto [t_enter, t_leave] = tile_span(q, dq);
        t_leave = dr::minimum(t_leave, t_edge);

        Float mint   = current.mint;
        Mask in_tile = mint >= t_enter && mint < t_leave;
        Mask before  = mint < t_enter && t_enter < t_leave;
        Float gap_end = dr::clip(dr::select(before, t_enter, t_edge), mint,
                                 current.maxt);

        ExtremumSegment segment(mint, gap_end, Vector2f(0.f));
        DDAState next = current;
        next.mint = gap_end;

        Mask descend = in_tile && active;
        if (dr::any_or<true>(descend)) {
            DDAState inner = current;
            inner.maxt = dr::minimum(current.maxt, t_leave);
            auto [inner_segment, inner_next] = m_inner->dda_next(inner, descend);
            inner_next.maxt = current.maxt;
            dr::masked(segment, descend) = inner_segment;
            dr::masked(next, descend)    = inner_next;
        }

        return { segment, next };
    }

    ScalarAffineTransform4f dda_to_local() const override {
        NotImplementedError("dda_to_local");
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
    /// Per-axis exit distance of the ray (q, dq) from the unit lattice cell.
    Vector3f cell_exit(const Point3f &q, const Vector3f &dq) const {
        Vector3f face = dr::select(dq > 0.f, 1.f, 0.f);
        return dr::select(dq != 0.f, (face - q) / dq, dr::Infinity<Float>);
    }

    /// Entry and exit distances of the ray (q, dq) through the inner domain.
    std::pair<Float, Float> tile_span(const Point3f &q,
                                      const Vector3f &dq) const {
        Vector3f t_lo = -q / dq,
                 t_hi = (m_tile_extents - q) / dq;
        Vector3f t_near = dr::minimum(t_lo, t_hi),
                 t_far  = dr::maximum(t_lo, t_hi);

        auto flat = dq == 0.f;
        auto flat_inside = q >= 0.f && q < m_tile_extents;
        dr::masked(t_near, flat) =
            dr::select(flat_inside, -dr::Infinity<Float>, dr::Infinity<Float>);
        dr::masked(t_far, flat) =
            dr::select(flat_inside, dr::Infinity<Float>, -dr::Infinity<Float>);

        return { dr::max(t_near), dr::min(t_far) };
    }

private:
    ref<ExtremumStructure> m_inner;
    ScalarVector3f m_lattice = ScalarVector3f(1.f);
    ScalarVector3f m_lattice_rcp = ScalarVector3f(1.f);
    ScalarVector3f m_tile_extents = ScalarVector3f(1.f);
    ScalarPoint3f m_origin = 0.f;
    ScalarAffineTransform4f m_to_lattice;
};

MI_EXPORT_PLUGIN(ExtremumRepeat)
NAMESPACE_END(mitsuba)
