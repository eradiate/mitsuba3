#include <mitsuba/core/eradiate/obb.h>
#include <mitsuba/python/python.h>

template <typename OBB>
void bind_obb(nb::module_ &m, const char *name) {
    using BBox3      = typename OBB::BBox3f;
    using Transform4 = typename OBB::Transform4f;

    MI_PY_CHECK_ALIAS(OBB, name) {
        nb::class_<OBB>(m, name)
            .def(nb::init<const BBox3 &, const Transform4 &>(),
                 "bbox"_a, "to_world"_a)
            .def("overlaps",
                 &OBB::overlaps,
                 "other"_a)
            .def_rw("center",       &OBB::center)
            .def_rw("rotation",     &OBB::rotation)
            .def_rw("half_extents", &OBB::half_extents);
    }
}

MI_PY_EXPORT(OrientedBoundingBox) {
    MI_PY_IMPORT_TYPES()

    using ScalarOBB = OrientedBoundingBox<ScalarFloat>;
    bind_obb<ScalarOBB>(m, "ScalarOrientedBoundingBox3f");

    if constexpr (!std::is_same_v<Float, ScalarFloat>) {
        using OBB = OrientedBoundingBox<Float>;
        bind_obb<OBB>(m, "OrientedBoundingBox3f");
    } else {
        m.attr("OrientedBoundingBox3f") = m.attr("ScalarOrientedBoundingBox3f");
    }
}
