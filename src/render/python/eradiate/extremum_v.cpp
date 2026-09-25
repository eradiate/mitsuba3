#include <mitsuba/core/properties.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/eradiate/extremum.h>
#include <mitsuba/render/eradiate/extremum_segment.h>
#include <mitsuba/render/eradiate/dda.h>
#include <mitsuba/python/python.h>
#include <nanobind/trampoline.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
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

MI_PY_EXPORT(DDAState) {
    MI_PY_IMPORT_TYPES()

    auto ds = nb::class_<DDAState>(m, "DDAState", D(DDAState))
        .def(nb::init<>())
        .def(nb::init<const DDAState &>(), "other"_a, "Copy constructor")
        .def_field(DDAState, o,    D(DDAState, o))
        .def_field(DDAState, d,    D(DDAState, d))
        .def_field(DDAState, pi,   D(DDAState, pi))
        .def_field(DDAState, mint, D(DDAState, mint))
        .def_field(DDAState, maxt, D(DDAState, maxt));

    MI_PY_DRJIT_STRUCT(ds, DDAState, o, d, pi, mint, maxt);

    nb::class_<DDAStateList>(m, "DDAStateList", D(DDAStateList))
        .def_field(DDAStateList, mint, D(DDAStateList, mint))
        .def_field(DDAStateList, maxt, D(DDAStateList, maxt));
}

/// Trampoline for derived types implemented in Python
MI_VARIANT class PyExtremumStructure : public ExtremumStructure<Float, Spectrum> {
public:
    MI_IMPORT_TYPES(ExtremumStructure, Volume)
    NB_TRAMPOLINE(ExtremumStructure, 7);

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

    DDAState dda_init(const Ray3f &ray, Float mint, Float maxt,
                      Mask active) const override {
        NB_OVERRIDE_PURE(dda_init, ray, mint, maxt, active);
    }

    std::pair<ExtremumSegment, DDAState>
    dda_next(const DDAState &state, Mask active) const override {
        NB_OVERRIDE_PURE(dda_next, state, active);
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

MI_PY_EXPORT(ExtremumStructure) {
    MI_PY_IMPORT_TYPES(ExtremumStructure)
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
        .def("bbox", &ExtremumStructure::bbox, D(ExtremumStructure, bbox))
        .def("dda_init", &ExtremumStructure::dda_init,
             "ray"_a, "mint"_a, "maxt"_a, "active"_a = true,
             D(ExtremumStructure, dda_init))
        .def("dda_next", &ExtremumStructure::dda_next,
             "state"_a, "active"_a = true,
             D(ExtremumStructure, dda_next))
        .def("eval_1", &ExtremumStructure::eval_1,
             "it"_a, "active"_a = true,
             D(ExtremumStructure, eval_1));

    drjit::bind_traverse(extremum);
}
