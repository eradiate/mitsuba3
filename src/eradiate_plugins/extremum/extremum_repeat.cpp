#include <mitsuba/core/properties.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/eradiate/extremum.h>

NAMESPACE_BEGIN(mitsuba)

/**!
.. _extremum-extremum_repeat:

Extremum repeat structure (:monosp:`extremum_repeat`)
-----------------------------------------------------

.. pluginparameters::

 * - (Nested plugin)
   - |extremum|
   - The inner extremum structure that is periodically repeated.

 * - lattice
   - |vector|
   - Axis-aligned lattice period along x, y, z. Default: the extents of the
     nested structure's bounding box.

 * - bbox_min, bbox_max
   - |point|
   - Bounds of the tiled domain in world space. Required.

This plugin repeats an extremum structure over the cells of a translation
lattice. It is created by the :ref:`repeat <medium-repeat>` medium for each of
its inner medium's structures and is not meant to be declared in a scene.
The inner medium and integrator are unaware of the tiling.
 */

template <typename Float, typename Spectrum>
class ExtremumRepeat final : public ExtremumStructure<Float, Spectrum> {
public:
    MI_IMPORT_BASE(ExtremumStructure, m_bbox)
    MI_IMPORT_TYPES(Volume)

    ExtremumRepeat(const Properties &props) : Base(props) {
        for (auto &prop : props.objects()) {
            if (auto *inner = prop.try_get<Base>()) {
                if (m_inner)
                    Throw("extremum_repeat accepts a single nested structure");
                m_inner = inner;
            }
        }
        if (!m_inner)
            Throw("extremum_repeat requires a nested extremum structure");

        ScalarBoundingBox3f tile = m_inner->bbox();
        ScalarVector3f lattice = props.get<ScalarVector3f>("lattice", tile.extents());
        m_cell = ScalarBoundingBox3f(tile.min, tile.min + lattice);
        m_bbox = ScalarBoundingBox3f(props.get<ScalarPoint3f>("bbox_min"),
                                     props.get<ScalarPoint3f>("bbox_max"));
    }

    void build(const Volume *) override { NotImplementedError("build"); }

    ExtremumState begin(const Ray3f &ray, Float mint, Float maxt,
                        Mask active) const override {
        auto [hit, d0, d1] = m_bbox.ray_intersect(ray);
        mint = dr::select(hit, dr::maximum(mint, d0), mint);
        maxt = dr::select(hit, dr::minimum(maxt, d1), mint);
        return enter_cell(ray, mint, maxt, active && mint < maxt);
    }

    ExtremumSegment next(ExtremumState &state, const Ray3f &ray,
                         Mask active) const override {
        Mask in_cell = state.mint < state.inner_maxt;
        Mask step    = active && in_cell;
        Mask move    = active && !in_cell;

        ExtremumSegment segment = dr::zeros<ExtremumSegment>();
        if (dr::any_or<true>(step)) {
            Float maxt = state.maxt;
            dr::masked(state.maxt, step) = state.inner_maxt;
            dr::masked(segment, step) = m_inner->next(state, ray, step);
            dr::masked(state.maxt, step) = maxt;
        }

        if (dr::any_or<true>(move)) {
            Float t0 = state.mint;
            dr::masked(state, move) = enter_cell(ray, t0, state.maxt, move);
            dr::masked(segment, move) = ExtremumSegment(t0, state.mint, 0.f, 0.f);
        }

        return segment;
    }

    std::tuple<Float, Float> eval_1(const Interaction3f &,
                                    Mask) const override {
        NotImplementedError("eval_1");
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "ExtremumRepeat[" << std::endl
            << "  inner = " << string::indent(m_inner) << "," << std::endl
            << "  cell = " << m_cell << "," << std::endl
            << "  bbox = " << m_bbox << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS(ExtremumRepeat)

private:
    /**
     * \brief State moved to the lattice cell the ray is in just after
     * \c t0, its inner part clipped to that cell and \c maxt.
     *
     * A point on a cell face belongs to the cell the ray is heading into. A
     * state left with no inner segment in the cell skips to the cell's exit,
     * or to \c maxt on inactive lanes.
     */
    ExtremumState enter_cell(const Ray3f &ray, Float t0, Float maxt,
                                Mask active) const {
        ScalarVector3f lattice = m_cell.extents();

        Vector3f exit_face = dr::select(ray.d > 0.f, 1.f, 0.f);
        Vector3f inv_d     = dr::rcp(ray.d);

        auto cell_exit = [&](const Vector3f &cell) {
            Vector3f t = (m_cell.min + (cell + exit_face) * lattice - ray.o) * inv_d;
            return dr::select(ray.d == 0.f, dr::Infinity<Float>, t);
        };

        // Get next intersection distance, account for precision error.
        Vector3f cell   = dr::floor((ray(t0) - m_cell.min) / lattice);
        Vector3f t_exit = cell_exit(cell);

        auto behind = t_exit <= t0;
        cell  += dr::select(behind, dr::sign(ray.d), 0.f);
        t_exit = dr::select(behind, cell_exit(cell), t_exit);

        Float end    = dr::minimum(dr::min(t_exit), maxt);
        Mask entered = t0 < end && active;

        // Begin inner state
        Ray3f shifted = ray;
        shifted.o     = ray.o - cell * lattice;
        ExtremumState state = m_inner->begin(shifted, t0, end, entered);

        // Inner is not aware of the bbox maxt, set back to current maxt.
        state.inner_maxt = state.maxt;
        state.maxt       = maxt;

        Float skip_to = dr::select(entered, end, maxt);
        Mask empty    = state.mint >= state.inner_maxt;
        dr::masked(state.mint, empty)       = skip_to;
        dr::masked(state.inner_maxt, empty) = skip_to;
        return state;
    }

    ref<const Base> m_inner;
    ScalarBoundingBox3f m_cell;
};

MI_EXPORT_PLUGIN(ExtremumRepeat)
NAMESPACE_END(mitsuba)
