#include <mitsuba/core/properties.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/eradiate/extremum.h>
#include <mitsuba/render/volume.h>

NAMESPACE_BEGIN(mitsuba)

/**!
.. _extremum-extremum_grid:

Extremum global structure (:monosp:`extremum_global`)
-------------------------------------------------

.. pluginparameters::

 * - volume
   - |volume|
   - Extinction coefficient volume to build extremum grid from
   - |exposed|

 * - scale
   - |float|
   - Scale factor for the extremum values (Default: 1.0).

This plugin holds the global minorant and majorant values of a volume.
At runtime, traversal is performed via a single segment determined by the 
passed ``mint`` and ``maxt`` values.
*/

template <typename Float, typename Spectrum>
class ExtremumGlobal final : public ExtremumStructure<Float, Spectrum> {
public:
    MI_IMPORT_BASE(ExtremumStructure, m_bbox)
    MI_IMPORT_TYPES(Volume)

    using TrackingState    = TrackingState<Float, Spectrum>;
    using TrackingFunction = TrackingFunction<Float, Spectrum>;

    ExtremumGlobal(const Properties &props) : Base(props) {
        // Volume Parameters
        m_volume = nullptr;
        for (auto &prop : props.objects()) {
            if (auto *vol = prop.try_get<Volume>()) {
                m_volume = vol;
                break;
            }
        }

        if (!m_volume)
            Throw("ExtremumGlobal requires at least one volume");
        
        // Register the extremum structure to the volume to trigger 
        // parameter_changed when the volume is modified.
        m_volume->add_extremum_structure(this);
        m_scale = props.get<ScalarFloat>("scale", 1.0f);
        m_bbox = m_volume->bbox();

        m_majorant = m_volume->max();
        m_minorant = m_volume->min();
    }

    void parameters_changed(const std::vector<std::string> &/*keys*/ = {}) override {
        m_bbox = m_volume->bbox();
        m_majorant = m_volume->max();
        m_minorant = m_volume->min();
    }


    TrackingState traverse_extremum(
        const Ray3f &ray,
        Float mint,
        Float maxt,
        UInt32 channel,
        TrackingState state,
        TrackingFunction* func,
        Mask active
    ) const override {

        // Clip the tracking interval to the volume's bounding box so that
        // get_scattering_coefficients is never queried outside the volume's
        // domain.  This matters when the containing shape is larger than the
        // volume (the medium's AABB != the volume's AABB).
        auto [bbox_hit, bbox_mint, bbox_maxt] = m_bbox.ray_intersect(ray);
        mint  = dr::maximum(mint,  bbox_mint);
        maxt  = dr::minimum(maxt,  bbox_maxt);
        active &= bbox_hit && (mint < maxt);

        ExtremumSegment segment(mint, maxt, m_scale*m_minorant, m_scale*m_majorant);

        struct LoopState {
            ExtremumSegment segment;
            TrackingState state;
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
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS(ExtremumGlobal)

private:
    ref<Volume> m_volume;
    ScalarFloat m_scale;

    ScalarFloat m_minorant;
    ScalarFloat m_majorant;
};

MI_EXPORT_PLUGIN(ExtremumGlobal)
NAMESPACE_END(mitsuba)
