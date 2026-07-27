#include <mitsuba/core/properties.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/eradiate/extremum.h>
#include <mitsuba/render/volume.h>

NAMESPACE_BEGIN(mitsuba)

/**!
.. _extremum-extremum_overlap:

Extremum overlap structure (:monosp:`extremum_overlap`)
-------------------------------------------------------

.. pluginparameters::

 * - (Nested plugin)
   - |extremum|
   - Two or more nested extremum structures to aggregate.

This plugin aggregates the extremum structures of overlapping media. It is
assembled — typically by a :monosp:`multicomponent` medium — from its
children's *built* structures; \ref ExtremumStructure::build() is never
called on it.

Traversal advances through the union of the children's segment boundaries:
at each step, the segment exit is the earliest child exit and the bounds are
the sums of the children's bounds. Summed majorants are valid majorants of
the summed extinction (and likewise for minorants), so the tracking callback
contract is unchanged. Children return empty (zero-valued) segments outside
their domain, which makes entries, exits and gaps between disjoint children
fall out of the same minimum.
*/

template <typename Float, typename Spectrum>
class ExtremumOverlap final : public ExtremumStructure<Float, Spectrum> {
public:
    MI_IMPORT_BASE(ExtremumStructure, m_bbox)
    MI_IMPORT_TYPES(Volume, ExtremumStructure)

    using TrackingStateType    = TrackingState<Float, Spectrum>;
    using TrackingFunctionType = TrackingFunction<Float, Spectrum>;

    struct MediumContext {
        Mask active;
        ExtremumSegment segment;

        DRJIT_STRUCT(MediumContext, active, segment)
    };

    using MediumContextArr = DynamicBuffer<MediumContext>;

    ExtremumOverlap(const Properties &props) : Base(props) {
        for (auto &prop : props.objects()) {
            if (auto *child = prop.try_get<ExtremumStructure>())
                m_children.push_back(child);
        }

        if (m_children.empty())
            Throw("extremum_overlap requires at least one nested extremum "
                  "structure");

        m_bbox.reset();
        for (const auto &child : m_children)
            m_bbox.expand(child->bbox());
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
        active &= maxt > mint;

        MediumContextArr ctx_arr = dr::zeros<MediumContextArr>(m_children.size());

        for (size_t i = 0; i < m_children.size(); ++i) {
            auto [hit, child_mint, child_maxt] = m_children[i]->bbox().ray_intersect(ray);
            ctx_arr[i].active = active && hit;
            ctx_arr[i].segment = ExtremumSegment(mint,
                dr::select(active && hit, mint, dr::Infinity<Float>),
                0.f, 0.f);
        }

        struct LoopState {
            ExtremumSegment segment;
            TrackingStateType state;
            Mask advance;
            Mask active;
            Float t;
            MediumContextArr ctx_arr;

            DRJIT_STRUCT(LoopState, segment, state, advance, active, t, ctx_arr)
        } ls = {
            dr::zeros<ExtremumSegment>(),
            state,
            /*advance=*/active,
            active,
            mint,
            ctx_arr
        };

        dr::tie(ls) = dr::while_loop(
            dr::make_tuple(ls),
            [](const LoopState &ls) { return ls.active; },
            [this, func, ray, maxt, channel](LoopState &ls) {

            // expose combined segment directly in the function
            Float seg_maxt = maxt;
            Vector2f value(0.f);

            // only recompute the segment when requested i.e. advance is true.
            if (dr::any_or<true>(ls.advance)) {

                for (size_t i = 0; i < m_children.size(); ++i) {
                    const auto &child = m_children[i];
                    auto &ctx = ls.ctx_arr[i];

                    Float eps = dr::maximum(dr::abs(ls.t), 1.f) * math::RayEpsilon<Float>;
                    Mask get_next_segment = ctx.active && (ls.t + eps > ctx.segment.maxt);

                    if (dr::any_or<true>(get_next_segment)) {
                        ctx.segment = child->next_segment(ray, ls.t, ls.active);
                    }

                    seg_maxt = dr::minimum(seg_maxt, ctx.segment.maxt);
                    value += ctx.segment.value;
                }

                ls.segment = ExtremumSegment(ls.t, dr::maximum(seg_maxt, ls.t), value);
            }

            // Combined segment at ls.t. When the callback does not advance
            // (null collision within the segment), ls.t is unchanged and the
            // recomputation yields the same segment — the re-feed contract.
            // dr::masked(ls.segment, ls.active) =
            //     combined_segment(ray, ls.t, maxt, ls.active);

            std::tie(ls.advance, ls.active) =
                func(ls.segment, ls.state, channel, ls.active);

            dr::masked(ls.t, ls.advance) = ls.segment.maxt;
            ls.active &= ls.t < maxt;
        },
        "Overlap Traversal");

        return ls.state;
    }

    ExtremumSegment next_segment(const Ray3f &ray, Float t,
                                 Mask active) const override {
        return combined_segment(ray, t, dr::Infinity<Float>, active);
    }

    std::tuple<Float, Float> eval_1(const Interaction3f &it,
                                    Mask active) const override {
        Float minorant = 0.f, majorant = 0.f;
        for (const auto &child : m_children) {
            Mask inside = BoundingBox3f(child->bbox()).contains(it.p);
            auto [mn, mx] = child->eval_1(it, active && inside);
            minorant += dr::select(inside, mn, 0.f);
            majorant += dr::select(inside, mx, 0.f);
        }
        return { minorant, majorant };
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "ExtremumOverlap[" << std::endl
            << "  children = " << m_children.size() << "," << std::endl
            << "  bbox = " << m_bbox << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS(ExtremumOverlap)

private:
    /// Aggregate the children's segments at distance t: earliest exit,
    /// summed bounds
    ExtremumSegment combined_segment(const Ray3f &ray, const Float &t,
                                     const Float &maxt, Mask active) const {
        Float seg_maxt = maxt;
        Vector2f value(0.f);
        for (const auto &child : m_children) {
            ExtremumSegment s = child->next_segment(ray, t, active);
            seg_maxt = dr::minimum(seg_maxt, s.maxt);
            value += s.value;
        }
        return ExtremumSegment(t, dr::maximum(seg_maxt, t), value);
    }

private:
    std::vector<ref<ExtremumStructure>> m_children;
};

MI_EXPORT_PLUGIN(ExtremumOverlap)
NAMESPACE_END(mitsuba)
