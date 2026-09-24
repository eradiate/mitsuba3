#include <mitsuba/core/properties.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/render/scene.h>
// #ERADIATE_CHANGE_BEGIN: Local extremum structure
#include <mitsuba/render/eradiate/extremum.h>
// #ERADIATE_CHANGE_END
#include <mitsuba/python/python.h>
#include <nanobind/trampoline.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>
#include <drjit/python.h>
// #ERADIATE_CHANGE_BEGIN: DDA support
#include <drjit/while_loop.h>
// #ERADIATE_CHANGE_END

/// Trampoline for derived types implemented in Python
MI_VARIANT class PyMedium : public Medium<Float, Spectrum> {
public:
    MI_IMPORT_TYPES(Medium, Sampler, Scene)
// #ERADIATE_CHANGE_BEGIN: Overlapping Media
    NB_TRAMPOLINE(Medium, 7);
// #ERADIATE_CHANGE_END

    PyMedium(const Properties &props) : Medium(props) {}

    std::tuple<Mask, Float, Float> intersect_aabb(const Ray3f &ray) const override {
        NB_OVERRIDE_PURE(intersect_aabb, ray);
    }

// #ERADIATE_CHANGE_BEGIN: Overlapping media
    Mask in_aabb(const Point3f &pos) const override {
        NB_OVERRIDE_PURE(in_aabb, pos);
    }
// #ERADIATE_CHANGE_END

    UnpolarizedSpectrum get_majorant(const MediumInteraction3f &mi, Mask active = true) const override {
        NB_OVERRIDE_PURE(get_majorant, mi, active);
    }

    std::tuple<UnpolarizedSpectrum, UnpolarizedSpectrum, UnpolarizedSpectrum>
    get_scattering_coefficients(const MediumInteraction3f &mi, Mask active = true) const override {
        NB_OVERRIDE_PURE(get_scattering_coefficients, mi, active);
    }

    std::string to_string() const override {
        NB_OVERRIDE_PURE(to_string);
    }

    void traverse(TraversalCallback *cb) override {
        NB_OVERRIDE(traverse, cb);
    }

    void parameters_changed(const std::vector<std::string> &keys) override {
        NB_OVERRIDE(parameters_changed, keys);
    }

    using Medium::m_sample_emitters;
    using Medium::m_is_homogeneous;
    using Medium::m_has_spectral_extinction;

    DR_TRAMPOLINE_TRAVERSE_CB(Medium)
};

template <typename Ptr, typename Cls> void bind_medium_generic(Cls &cls) {
    MI_PY_IMPORT_TYPES(PhaseFunctionContext)

    using RetPhaseFunction = std::conditional_t<drjit::is_array_v<Ptr>, PhaseFunctionPtr, drjit::scalar_t<PhaseFunctionPtr>>;

    cls.def("phase_function",
            [](Ptr ptr) -> RetPhaseFunction { return ptr->phase_function(); },
            D(Medium, phase_function))
       .def("use_emitter_sampling",
            [](Ptr ptr) { return ptr->use_emitter_sampling(); },
            D(Medium, use_emitter_sampling))
       .def("is_homogeneous",
            [](Ptr ptr) { return ptr->is_homogeneous(); },
            D(Medium, is_homogeneous))
       .def("has_spectral_extinction",
            [](Ptr ptr) { return ptr->has_spectral_extinction(); },
            D(Medium, has_spectral_extinction))
       .def("get_majorant",
            [](Ptr ptr, const MediumInteraction3f &mi, Mask active) {
                return ptr->get_majorant(mi, active); },
            "mi"_a, "active"_a=true,
            D(Medium, get_majorant))
       .def("intersect_aabb",
            [](Ptr ptr, const Ray3f &ray) {
                return ptr->intersect_aabb(ray); },
            "ray"_a,
            D(Medium, intersect_aabb))
// #ERADIATE_CHANGE_BEGIN: Overlapping media
        .def("in_aabb",
            [](Ptr ptr, const Point3f &pos) {
                return ptr->in_aabb(pos); },
            "pos"_a,
            D(Medium, in_aabb))
        .def("phase_function",
            [](Ptr ptr, const UInt32& component, Mask active) {
                return ptr->phase_function(component, active); },
            "component"_a, "active"_a,
            D(Medium, phase_function))
// #ERADIATE_CHANGE_END
       .def("sample_interaction",
            [](Ptr ptr, const Ray3f &ray, Float sample, UInt32 channel, Mask active) {
                return ptr->sample_interaction(ray, sample, channel, active); },
            "ray"_a, "sample"_a, "channel"_a, "active"_a,
            D(Medium, sample_interaction))
       .def("transmittance_eval_pdf",
            [](Ptr ptr, const MediumInteraction3f &mi,
               const SurfaceInteraction3f &si, Mask active) {
                return ptr->transmittance_eval_pdf(mi, si, active); },
            "mi"_a, "si"_a, "active"_a,
            D(Medium, transmittance_eval_pdf))
// #ERADIATE_CHANGE_BEGIN: Add function that calculates the transmittance and pdf
        .def("sample_interaction_analytical",
            [](Ptr ptr, const Ray3f &ray, const Interaction3f &it, Float sample, UInt32 channel, Mask active) {
                return ptr->sample_interaction_analytical(ray, it, sample, channel, active); },
            "ray"_a, "it"_a, "sample"_a, "channel"_a, "active"_a,
            D(Medium, sample_interaction_analytical))
        .def("transmittance_eval_analytical",
            [](Ptr ptr, const Ray3f &ray,
                const Interaction3f &it, Mask active) {
                return ptr->transmittance_eval_analytical(ray, it, active); },
            "ray"_a, "it"_a, "active"_a,
            D(Medium, transmittance_eval_analytical))
// #ERADIATE_CHANGE_END
       .def("get_scattering_coefficients",
            [](Ptr ptr, const MediumInteraction3f &mi, Mask active = true) {
                return ptr->get_scattering_coefficients(mi, active); },
            "mi"_a, "active"_a=true,
            D(Medium, get_scattering_coefficients));
// #ERADIATE_CHANGE_BEGIN: DDA support
    cls.def("track_test",
            [](Ptr ptr, const Ray3f &ray, UInt32 seed, bool ratio,
               bool use_dda, Mask active) {
                using TrackingStateType = TrackingState<Float, Spectrum>;

                auto [mei, mint, maxt] =
                    ptr->prepare_medium_traversal(ray, active);
                active &= dr::isfinite(maxt) && mint < maxt;

                dr::PCG32<UInt32> rng;
                rng.seed(rng.PCG32_DEFAULT_STATE, dr::uint64_array_t<Float>(seed));
                Float target_ot =
                    -dr::log(1.f - rng.template next_float<Float>(active));

                TrackingStateType state{ ray, rng, mei, target_ot,
                                         ptr->use_rrt(),
                                         ptr->has_spectral_extinction(),
                                         0u, UnpolarizedSpectrum(1.f) };
                auto func = ratio ? ratio_track_segment<Float, Spectrum>
                                  : delta_track_segment<Float, Spectrum>;

                if (use_dda)
                    dr::masked(state, active) = ptr->dda_track(
                        ray, mint, maxt, state, UInt32(0), func, active);
                else
                    dr::masked(state, active) =
                        ptr->extremum_structure()->traverse_extremum(
                            ray, mint, maxt, UInt32(0), state, func, active);

                return std::make_tuple(state.mei.t, state.throughput);
            },
            "ray"_a, "seed"_a, "ratio"_a, "use_dda"_a, "active"_a = true,
            "Test utility: delta (or ratio) tracking along `ray`, through "
            "`dda_track` or `traverse_extremum`. Returns (distance, "
            "throughput).");
// #ERADIATE_CHANGE_END
}

MI_PY_EXPORT(Medium) {
    MI_PY_IMPORT_TYPES(Medium, MediumPtr, Scene, Sampler)
    using PyMedium = PyMedium<Float, Spectrum>;
    using Properties = mitsuba::Properties;

    auto medium = MI_PY_TRAMPOLINE_CLASS(PyMedium, Medium, Object)
        .def(nb::init<const Properties &>(), "props"_a)
        .def_field(PyMedium, m_sample_emitters, D(Medium, m_sample_emitters))
        .def_field(PyMedium, m_is_homogeneous, D(Medium, m_is_homogeneous))
        .def_field(PyMedium, m_has_spectral_extinction, D(Medium, m_has_spectral_extinction))
// #ERADIATE_CHANGE_BEGIN: Local extremum structure
        .def("extremum_structure",
             [](Medium *ptr) { return ptr->extremum_structure(); },
             D(Medium, extremum_structure))
// #ERADIATE_CHANGE_END
        .def("__repr__", &Medium::to_string, D(Medium, to_string));

    drjit::bind_traverse(medium);

    bind_medium_generic<Medium *>(medium);
// #ERADIATE_CHANGE_BEGIN: DDA support
    // Host pointer only: `dda_init` / `dda_step` are not in the vcall block.
    medium.def("sample_test_dda",
            [](Medium *ptr, const Ray3f &ray, Float mint, Float maxt,
               Float target_ot, Mask active) {
                struct LoopState {
                    DDAStateList dda;
                    Float t;
                    Float target_ot;
                    Mask active;

                    DRJIT_STRUCT(LoopState, dda, t, target_ot, active)
                };

                DDAStateList dda = ptr->dda_init(ray, mint, maxt, active);
                LoopState ls = { dda, dr::Infinity<Float>, target_ot,
                                 active && (dda.mint < dda.maxt) };

                dr::tie(ls) = dr::while_loop(
                    dr::make_tuple(ls),
                    [](const LoopState &ls) { return ls.active; },
                    [ptr, &ray](LoopState &ls) {
                        ExtremumSegment segment =
                            ptr->dda_step(ls.dda, ray, ls.active);

                        Float segment_ot = (segment.maxt - segment.mint) *
                                           segment.majorant();
                        Mask sampled = (ls.target_ot < segment_ot) && ls.active;

                        dr::masked(ls.t, sampled) =
                            segment.mint +
                            ls.target_ot / dr::maximum(segment.majorant(),
                                                       dr::Epsilon<Float>);
                        dr::masked(ls.target_ot, !sampled && ls.active) -=
                            segment_ot;

                        ls.active &= !sampled && (ls.dda.mint < ls.dda.maxt);
                    },
                    "sample_test_dda");

                return std::make_tuple(ls.t, ls.target_ot);
            },
            "ray"_a, "mint"_a, "maxt"_a, "target_ot"_a, "active"_a = true,
            "Test utility: walk the medium's DDA route, accumulating the "
            "majorant optical thickness until `target_ot` is reached. Returns "
            "(distance, leftover_ot); `distance` is infinite if `target_ot` "
            "is not reached before `maxt`.");
// #ERADIATE_CHANGE_END

    if constexpr (dr::is_array_v<MediumPtr>) {
        dr::ArrayBinding b;
        auto medium_ptr = dr::bind_array_t<MediumPtr>(b, m, "MediumPtr");
        bind_medium_generic<MediumPtr>(medium_ptr);
    }

}
