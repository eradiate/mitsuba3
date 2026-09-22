#include <mitsuba/core/properties.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/eradiate/extremum.h>
#include <mitsuba/render/volume.h>

NAMESPACE_BEGIN(mitsuba)

/**!
.. _extremum-extremum_global:

Extremum global structure (:monosp:`extremum_global`)
-----------------------------------------------------

This plugin holds the global minorant and majorant values of a volume.
At runtime, traversal is performed via a single segment determined by the
passed ``mint`` and ``maxt`` values.
*/

template <typename Float, typename Spectrum>
class ExtremumGlobal final : public ExtremumStructure<Float, Spectrum> {
public:
    MI_IMPORT_BASE(ExtremumStructure, m_bbox, m_scale)
    MI_IMPORT_TYPES(Volume)

    using TrackingStateType    = TrackingState<Float, Spectrum>;
    using TrackingFunctionType = TrackingFunction<Float, Spectrum>;

    ExtremumGlobal(const Properties &props) : Base(props) {}

    void build(const Volume *volume) override {
        m_minorant = volume->min();
        m_majorant = volume->max();
    }

    ExtremumSegment next_segment(const Ray3f &ray, Float t,
                                 Mask active) const override {
        auto [hit, d0, d1] = m_bbox.ray_intersect(ray);

        ExtremumSegment segment(t, dr::Infinity<Float>, Vector2f(0.f));

        Float eps = dr::maximum(dr::abs(t), 1.f) * math::RayEpsilon<Float>;
        Float tq  = t + eps;

        active &= hit;
        Mask before = hit && (tq < d0);
        Mask inside = hit && !before && (tq < d1);

        dr::masked(segment.maxt, before) = d0;

        // early exit
        if (dr::any_or<false>(!inside)) {
            return segment;
        }

        Float maxt = dr::select(
            inside, dr::maximum(d1, tq),
            dr::select(before, d0, dr::Infinity<Float>));

        Vector2f value = dr::select(
            inside, m_scale*Vector2f(m_minorant, m_majorant), Vector2f(0.f));

        return ExtremumSegment(t, maxt, value);
    }

    TrackingStateType traverse_extremum(
        const Ray3f &/*ray*/,
        Float mint,
        Float maxt,
        UInt32 channel,
        TrackingStateType state,
        TrackingFunctionType* func,
        Mask active
    ) const override {
        ExtremumSegment segment(mint, maxt, m_scale*m_minorant, m_scale*m_majorant);

        struct LoopState {
            ExtremumSegment segment;
            TrackingStateType state;
            Mask advance;
            Mask active;

            DRJIT_STRUCT(LoopState, segment, state, advance, active)
        } ls {
            segment,
            state,
            /*advance =*/true,
            active
        };

        dr::tie(ls) = dr::while_loop(
        dr::make_tuple(ls),
        [](const LoopState &ls){ return ls.active; },
        [func, channel](LoopState &ls){
            std::tie(ls.advance, ls.active) =
                func(ls.segment, ls.state, channel, ls.active);
            ls.active &= !ls.advance;
        });

        return ls.state;
    }

    std::tuple<Float, Float> eval_1(const Interaction3f &/*it*/,
                                    Mask /*active*/) const override {

        return { m_scale*m_minorant, m_scale*m_majorant };
    }

    void traverse(TraversalCallback *cb) override {
        Base::traverse(cb);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "ExtremumGlobal[" << std::endl
            << "  minorant = " << m_minorant << "," << std::endl
            << "  majorant = " << m_majorant << "," << std::endl
            << "  scale = " << m_scale << "," << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS(ExtremumGlobal)

private:
    ScalarFloat m_minorant;
    ScalarFloat m_majorant;
};

MI_EXPORT_PLUGIN(ExtremumGlobal)
NAMESPACE_END(mitsuba)
