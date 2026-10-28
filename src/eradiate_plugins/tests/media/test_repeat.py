import drjit as dr
import mitsuba as mi
import pytest

REGION = (-20.0, 20.0)


@pytest.fixture
def component(variant_scalar_rgb):
    """A medium filling [0, 1]^3."""
    return {
        "type": "eoheterogeneous",
        "sigma_t": {"type": "constvolume", "value": 1.0},
        "albedo": 0.8,
        "aabb_min": mi.ScalarPoint3f(0.0),
        "aabb_max": mi.ScalarPoint3f(1.0),
        "extremum": {"type": "globalextremum"},
    }


def _repeat(inner, lattice=None, region=REGION):
    d = {"type": "repeat", "inner": inner}
    if lattice is not None:
        d["lattice"] = mi.ScalarVector3f(lattice)
    if region is not None:
        d["aabb_min"] = mi.ScalarPoint3f(region[0])
        d["aabb_max"] = mi.ScalarPoint3f(region[1])
    return d


def test_repeat_wraps_each_structure(component):
    """``repeat`` tiles every structure of its inner medium over its region."""
    overlap = {"type": "overlap", "a": component, "b": component}
    structures = mi.load_dict(_repeat(overlap, (1.5, 1.5, 1.5))).extrema()
    assert len(structures) == 2
    for structure in structures:
        bbox = structure.bbox()
        assert dr.all(bbox.min == mi.ScalarPoint3f(REGION[0]))
        assert dr.all(bbox.max == mi.ScalarPoint3f(REGION[1]))


def test_repeat_rejects_invalid_tiling(component):
    with pytest.raises(RuntimeError, match="nest repeat media"):
        mi.load_dict(_repeat(_repeat(component)))
    with pytest.raises(RuntimeError, match="do not fit"):
        mi.load_dict(_repeat(component, (1.5, 0.5, 1.0)))
    with pytest.raises(RuntimeError, match="aabb_min"):
        mi.load_dict(_repeat(component, region=None))
