#!/usr/bin/env python
#
# Walks through all eradiate plugin files and extracts documentation
# for inclusion in the plugin reference.
# Adapted from docs/generate_plugin_doc.py.

import os
import re

BSDF_ORDERING = [
    "bilambertian",
    "hapke",
    "maignan",
    "measured_mono",
    "mqdiffuse",
    "ocean_grasp",
    "ocean_legacy",
    "ocean_mishchenko",
    "rpv",
    "rtls",
    "selectbsdf",
]

EMITTER_ORDERING = [
    "astroobject",
]

EXTREMUM_ORDERING = [
    "extremum_grid",
    "extremum_spherical",
]

INTEGRATOR_ORDERING = [
    "piecewise_volpath",
    "eovolpath",
]

MEDIA_ORDERING = [
    "piecewise",
]

PHASE_ORDERING = [
    "multiphase",
    "rayleigh_polarized",
    "tabphase_irregular",
    "tabphase_polarized",
]

SENSOR_ORDERING = [
    "distantflux",
    "hdistant",
    "mdistant",
    "mpdistant",
    "mradiancemeter",
]

SHAPE_ORDERING = [
    "arectangle",
    "instancelist",
]

VOLUME_ORDERING = [
    "sphericalcoords",
]


def find_order_id(filename, ordering):
    f = os.path.split(filename)[-1].split(".")[0]
    if ordering and f in ordering:
        return ordering.index(f)
    elif filename in ordering:
        return ordering.index(filename)
    else:
        return 1000


def extract(target, filename):
    f = open(filename, encoding="utf-8")
    inheader = False
    for line in f.readlines():
        match = re.match(r"^/\*\*! ?(.*)$", line)
        if match is not None:
            print("Processing %s" % filename)
            line = match.group(1).replace("%", "\%")
            target.write(line + "\n")
            inheader = True
            continue
        if not inheader:
            continue
        if re.search(r"^[\s\*]*\*/$", line):
            inheader = False
            continue
        target.write(line)
    f.close()


def extract_python(target, filename):
    f = open(filename, encoding="utf-8")
    inheader = False
    for line in f.readlines():
        # Remove indentation
        if line.startswith("    "):
            line = line[4:]
        match_beg = re.match(r"r\"\"\"", line)
        match_end = re.match(r"\"\"\"", line)
        if not inheader and match_beg is not None:
            print("Processing %s" % filename)
            inheader = True
            continue
        if inheader and match_end is not None:
            inheader = False
            continue
        if not inheader:
            continue
        target.write(line)
    f.close()


def process(path, target, ordering):
    def capture(file_list, dirname, files):
        suffix = os.path.split(dirname)[1]
        if (
            "lib" in suffix
            or suffix == "tests"
            or suffix == "mitsuba"
            or suffix == "utils"
            or suffix == "converter"
        ):
            return
        for filename in files:
            if ".cpp" == os.path.splitext(filename)[1]:
                fname = os.path.join(dirname, filename)
                file_list += [fname]

    file_list = []
    for dirname, _, files in os.walk(path):
        capture(file_list, dirname, files)

    for o in ordering:
        if o.endswith(".py"):
            file_list.append(o)

    ordering = [(find_order_id(fname, ordering), fname) for fname in file_list]
    ordering = sorted(ordering, key=lambda entry: entry[0])

    for entry in ordering:
        if entry[1].endswith(".py"):
            extract_python(target, entry[1])
        else:
            extract(target, entry[1])


def process_src(target, src_subdir, ordering=None):
    section = "section_" + src_subdir

    with open("plugin_reference/" + section + ".rst", "r", encoding="utf-8") as f:
        target.write(f.read())
    process("../src/eradiate_plugins/{0}".format(src_subdir), target, ordering)


def generate(build_dir):
    original_wd = os.getcwd()
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    sections = [
        ("bsdfs", BSDF_ORDERING),
        ("emitters", EMITTER_ORDERING),
        ("extremum", EXTREMUM_ORDERING),
        ("integrators", INTEGRATOR_ORDERING),
        ("media", MEDIA_ORDERING),
        ("phase", PHASE_ORDERING),
        ("sensors", SENSOR_ORDERING),
        ("shapes", SHAPE_ORDERING),
        ("volumes", VOLUME_ORDERING),
    ]

    for section, ordering in sections:
        with open(
            os.path.join(build_dir, f"plugins_{section}.rst"), "w", encoding="utf-8"
        ) as f:
            process_src(f, section, ordering)

    os.chdir(original_wd)
