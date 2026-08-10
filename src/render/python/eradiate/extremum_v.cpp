#include <mitsuba/core/properties.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/eradiate/extremum.h>
#include <mitsuba/render/eradiate/extremum_segment.h>
#include <mitsuba/python/python.h>
#include <nanobind/trampoline.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>
#include <drjit/python.h>


MI_PY_EXPORT(ExtremumSegment) {
    MI_PY_IMPORT_TYPES()

    auto es = nb::class_<ExtremumSegment>(m, "ExtremumSegment", D(ExtremumSegment))
        .def(nb::init<>())
        .def(nb::init<const ExtremumSegment &>(), "other"_a, "Copy constructor")
        .def(nb::init<Float, Float, Float, Float>(),
                 D(ExtremumSegment, ExtremumSegment, 2),
                 "mint"_a, "maxt"_a, "minorant"_a, "majorant"_a)
        .def(nb::init<Float, Float, Vector2f>(),
                 D(ExtremumSegment, ExtremumSegment, 2),
                 "mint"_a, "maxt"_a, "value"_a)
        .def("valid",        &ExtremumSegment::valid,     D(ExtremumSegment, valid))
        .def("reset",        &ExtremumSegment::reset,     D(ExtremumSegment, reset))
        .def("zero_",        &ExtremumSegment::zero_,     "size"_a = 1)
        .def("zero_",        &ExtremumSegment::zero_,     D(ExtremumSegment, zero))
        .def("minorant",     &ExtremumSegment::minorant,  D(ExtremumSegment, minorant))
        .def("majorant",     &ExtremumSegment::majorant,  D(ExtremumSegment, majorant))
        .def_field(ExtremumSegment, mint,   D(ExtremumSegment, mint))
        .def_field(ExtremumSegment, maxt,   D(ExtremumSegment, maxt))
        .def_field(ExtremumSegment, value,  D(ExtremumSegment, value))
        .def_repr(ExtremumSegment);

    MI_PY_DRJIT_STRUCT(es, ExtremumSegment, mint, maxt, value);
}

/// Trampoline for derived types implemented in Python
// Note that `traverse_extremum` does not appear in this list. This is because
// it accepts a concrete function pointer as parameter, the binding of which
// has not been solved yet.
MI_VARIANT class PyExtremumStructure : public ExtremumStructure<Float, Spectrum> {
public:
    MI_IMPORT_TYPES(ExtremumStructure, Volume)
    NB_TRAMPOLINE(ExtremumStructure, 5);

    PyExtremumStructure(const Properties &props) : ExtremumStructure(props) {}

    void build(const Volume * volume) override {
        NB_OVERRIDE_PURE(build, volume);
    }

    std::tuple<Float, Float> eval_1(
        const Interaction3f &it,
        Mask active
    ) const override {
        NB_OVERRIDE_PURE(eval_1, it, active);
    }

    std::string to_string() const override {
        NB_OVERRIDE(to_string);
    }

    void traverse(TraversalCallback *cb) override {
        NB_OVERRIDE(traverse, cb);
    }

    void parameters_changed(const std::vector<std::string> &keys) override {
        NB_OVERRIDE(parameters_changed, keys);
    }

    DR_TRAMPOLINE_TRAVERSE_CB(ExtremumStructure)
};

template <typename Ptr, typename Cls> void bind_extremum_structure_generic(Cls &cls) {
    MI_PY_IMPORT_TYPES(ExtremumStructure, Medium)

    cls.def("eval_1",
            [](Ptr ptr, const Interaction3f &it, Mask active) {
                return ptr->eval_1(it, active);
            },
            "it"_a, "active"_a = true,
            D(ExtremumStructure, eval_1));


    // Test utility: deterministic delta tracking driven by a fixed target
    // optical thickness.
    cls.def("sample_test",
            [](Ptr ptr, const Ray3f &ray, Float mint, Float maxt,
               Float target_ot, UInt32 channel, Mask active) {
                using TrackingStateType = TrackingState<Float, Spectrum>;

                TrackingStateType state = dr::zeros<TrackingStateType>();
                state.ray       = ray;
                state.target_ot = target_ot;
                state.mei       = dr::zeros<MediumInteraction3f>();

                state = ptr->traverse_extremum(
                    ray, mint, maxt, channel, state,
                    [](const ExtremumSegment &segment, TrackingStateType &state,
                       const UInt32 &, Mask active) {
                        Float mint = dr::select(
                            state.mei.is_valid(),
                            dr::maximum(segment.mint, state.mei.t),
                            segment.mint);
                        Float segment_ot =
                            (segment.maxt - mint) * segment.majorant();
                        Mask sampled = (state.target_ot < segment_ot) && active;

                        Float maxt = dr::select(
                            sampled,
                            mint + state.target_ot /
                                       dr::maximum(segment.majorant(),
                                                   dr::Epsilon<Float>),
                            segment.maxt);

                        dr::masked(state.mei.t, sampled)  = maxt;
                        dr::masked(state.mei.t, !sampled) = dr::Infinity<Float>;
                        dr::masked(state.target_ot, !sampled && active) -=
                            segment_ot;

                        return std::pair<Mask, Mask>(/*advance=*/!sampled,
                                                     active && !sampled);
                    },
                    active);

                return std::make_tuple(state.mei.t, state.target_ot);
            },
            "ray"_a, "mint"_a, "maxt"_a, "target_ot"_a, "channel"_a = 0u,
            "active"_a = true,
            "Deterministic delta-tracking test utility. Traverses the extremum "
            "structure's segments, until an interaction is sampled based on "
            "`target_ot`. Returns (distance, leftover_ot); `distance` is "
            "infinite if `target_ot` is not reached before `maxt`.");
}


MI_PY_EXPORT(ExtremumStructure) {
    MI_PY_IMPORT_TYPES(ExtremumStructure, ExtremumStructurePtr)
    using PyExtremumStructure = PyExtremumStructure<Float, Spectrum>;
    using Properties = mitsuba::Properties;

    auto extremum = MI_PY_TRAMPOLINE_CLASS(PyExtremumStructure, ExtremumStructure, Object)
        .def(nb::init<const Properties &>(), "props"_a)
        .def("__repr__", &ExtremumStructure::to_string)
        .def("set_bbox", &ExtremumStructure::set_bbox,
             "bbox"_a, D(ExtremumStructure, set_bbox))
        .def("set_scale", &ExtremumStructure::set_scale,
             "scale"_a, D(ExtremumStructure, set_scale))
        .def("update_extremum", &ExtremumStructure::update_extremum,
             "bbox"_a, "volume"_a, "scale"_a = nb::none(),
             D(ExtremumStructure, update_extremum))
        .def("build", &ExtremumStructure::build,
             "volume"_a, D(ExtremumStructure, build))
        .def("bbox", &ExtremumStructure::bbox, D(ExtremumStructure, bbox));

    drjit::bind_traverse(extremum);

    bind_extremum_structure_generic<ExtremumStructure *>(extremum);

    if constexpr (dr::is_array_v<ExtremumStructurePtr>) {
        dr::ArrayBinding b;
        auto extremum_ptr = dr::bind_array_t<ExtremumStructurePtr>(b, m, "ExtremumStructurePtr");
        bind_extremum_structure_generic<ExtremumStructurePtr>(extremum_ptr);
    }
}
