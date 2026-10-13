#include <mitsuba/render/interaction.h>
#include <mitsuba/render/eradiate/tracking.h>
#include <mitsuba/python/python.h>
#include <drjit/python.h>

MI_PY_EXPORT(TrackingRecord) {
    MI_PY_IMPORT_TYPES()
    using TrackingRecord = mitsuba::TrackingRecord<Float, Spectrum>;

    auto tr = nb::class_<TrackingRecord>(m, "TrackingRecord", D(TrackingRecord))
        .def(nb::init<>())
        .def(nb::init<const TrackingRecord &>(), "other"_a, "Copy constructor")
        .def_field(TrackingRecord, mei,           D(TrackingRecord, mei))
        .def_field(TrackingRecord, component,     D(TrackingRecord, component))
        .def_field(TrackingRecord, throughput,    D(TrackingRecord, throughput))
        .def_field(TrackingRecord, transmittance,
                   D(TrackingRecord, transmittance))
        .def_field(TrackingRecord, pdf_ratio,     D(TrackingRecord, pdf_ratio));

    MI_PY_DRJIT_STRUCT(tr, TrackingRecord, mei, component, throughput,
                       transmittance, pdf_ratio);
}
