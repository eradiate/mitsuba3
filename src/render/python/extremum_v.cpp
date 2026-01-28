#include <mitsuba/core/properties.h>
#include <mitsuba/render/eradiate/extremum.h>
#include <mitsuba/python/python.h>
#include <nanobind/trampoline.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <drjit/python.h>

MI_PY_EXPORT(ExtremumSegment) {
    MI_PY_IMPORT_TYPES()

    auto es = nb::class_<ExtremumSegment>(m, "ExtremumSegment", D(ExtremumSegment))
        .def(nb::init<>())
        .def(nb::init<const ExtremumSegment &>(), "other"_a, "Copy constructor")
        .def_rw("tmin",      &ExtremumSegment::tmin,      D(ExtremumSegment, tmin))
        .def_rw("tmax",      &ExtremumSegment::tmax,      D(ExtremumSegment, tmax))
        .def_rw("sigma_maj", &ExtremumSegment::sigma_maj, D(ExtremumSegment, sigma_maj))
        .def_rw("sigma_min", &ExtremumSegment::sigma_min, D(ExtremumSegment, sigma_min))
        .def_repr(ExtremumSegment);

    MI_PY_DRJIT_STRUCT(es, ExtremumSegment, tmin, tmax, sigma_maj, sigma_min);
}

/// Trampoline for derived types implemented in Python
MI_VARIANT class PyExtremumStructure : public ExtremumStructure<Float, Spectrum> {
public:
    MI_IMPORT_TYPES(ExtremumStructure)
    NB_TRAMPOLINE(ExtremumStructure, 4);

    PyExtremumStructure(const Properties &props) : ExtremumStructure(props) {}

    ExtremumSegment sample_segment(
        const Ray3f &ray,
        Float mint,
        Float maxt,
        Float desired_tau,
        Mask active
    ) const override {
        NB_OVERRIDE_PURE(sample_segment, ray, mint, maxt, desired_tau, active);
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
    MI_PY_IMPORT_TYPES(ExtremumStructure)

    cls.def("sample_segment",
            [](Ptr ptr, const Ray3f &ray, Float mint, Float maxt,
               Float desired_tau, Mask active) {
                return ptr->sample_segment(ray, mint, maxt, desired_tau, active);
            },
            "ray"_a, "mint"_a, "maxt"_a, "desired_tau"_a, "active"_a = true,
            D(ExtremumStructure, sample_segment));
}



MI_PY_EXPORT(ExtremumStructure) {
    MI_PY_IMPORT_TYPES(ExtremumStructure, ExtremumStructurePtr)
    using PyExtremumStructure = PyExtremumStructure<Float, Spectrum>;
    using Properties = mitsuba::Properties;

    auto extremum = MI_PY_TRAMPOLINE_CLASS(PyExtremumStructure, ExtremumStructure, Object)
        .def(nb::init<const Properties &>(), "props"_a)
        .def("__repr__", &ExtremumStructure::to_string)
        .def("bbox", &ExtremumStructure::bbox, D(ExtremumStructure, bbox));

    drjit::bind_traverse(extremum);

    bind_extremum_structure_generic<ExtremumStructure *>(extremum);

    if constexpr (dr::is_array_v<ExtremumStructurePtr>) {
        dr::ArrayBinding b;
        auto extremum_ptr = dr::bind_array_t<ExtremumStructurePtr>(b, m, "ExtremumStructurePtr");
        bind_extremum_structure_generic<ExtremumStructurePtr>(extremum_ptr);
    }
}
