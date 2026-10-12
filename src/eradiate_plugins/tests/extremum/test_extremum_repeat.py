"""Lattice tiling of extremum structures by ``repeat`` media
(``extremum_repeat``)."""

import drjit as dr
import mitsuba as mi
import pytest

from .utils import FAR, walk_iter


def _tiled(lattice, region=FAR):
    """A [0, 1]^3 tile with majorant 1 on x < 0.5 and 2 above, repeated every
    ``lattice`` along each axis inside ``region``."""
    volume = mi.load_dict(
        {
            "type": "gridvolume",
            "grid": mi.VolumeGrid(mi.TensorXf([1.0, 2.0], shape=(1, 1, 2, 1))),
            "filter_type": "nearest",
        }
    )
    inner = mi.load_dict(
        {"type": "extremum_grid", "resolution": mi.ScalarVector3i(2, 1, 1)}
    )
    inner.update_extremum(volume.bbox(), volume)
    return mi.load_dict(
        {
            "type": "extremum_repeat",
            "inner": inner,
            "lattice": mi.ScalarVector3f(lattice),
            "bbox_min": mi.ScalarPoint3f(*region[0]),
            "bbox_max": mi.ScalarPoint3f(*region[1]),
        }
    )


# fmt: off
@pytest.mark.parametrize(
    "lattice, region, o, d, maxt, expected",
    [
        pytest.param(
            1.0, FAR, (0.25, 0.5, 0.5), (1, 0, 0), 2.0,
            [(0.0, 0.25, 1), (0.25, 0.75, 2), (0.75, 1.25, 1), (1.25, 1.75, 2),
             (1.75, 2.0, 1)],
            id="touching",
        ),
        pytest.param(
            1.5, FAR, (0.25, 0.5, 0.5), (1, 0, 0), 3.0,
            [(0.0, 0.25, 1), (0.25, 0.75, 2), (0.75, 1.25, 0), (1.25, 1.75, 1),
             (1.75, 2.25, 2), (2.25, 2.75, 0), (2.75, 3.0, 1)],
            id="gaps",
        ),
        pytest.param(
            1.5, FAR, (1.25, 0.5, 0.5), (1, 0, 0), 3.0,
            [(0.25, 0.75, 1), (0.75, 1.25, 2), (1.25, 1.75, 0), (1.75, 2.25, 1),
             (2.25, 2.75, 2), (2.75, 3.0, 0)],
            id="from_gap",
        ),
        pytest.param(
            1.5, FAR, (2.25, 0.5, 0.5), (-1, 0, 0), 3.0,
            [(0.0, 0.25, 2), (0.25, 0.75, 1), (0.75, 1.25, 0), (1.25, 1.75, 2),
             (1.75, 2.25, 1), (2.25, 2.75, 0), (2.75, 3.0, 2)],
            id="negative",
        ),
        pytest.param(
            1.5, FAR, (1.5, 0.5, 0.5), (1, 0, 0), 1.5,
            [(0.0, 0.5, 1), (0.5, 1.0, 2), (1.0, 1.5, 0)],
            id="on_face_forward",
        ),
        pytest.param(
            1.5, FAR, (1.5, 0.5, 0.5), (-1, 0, 0), 1.5,
            [(0.5, 1.0, 2), (1.0, 1.5, 1)],
            id="on_face_backward",
        ),
        pytest.param(
            1.5, FAR, (0.25, 1.25, 0.5), (1, 0, 0), 3.0,
            [(1.25, 2.75, 0), (2.75, 3.0, 0)],
            id="never_in_tile",
        ),
        pytest.param(
            1.5, ((-0.5, -10, -10), (2.0, 10, 10)), (-1.0, 0.5, 0.5), (1, 0, 0), 4.0,
            [(1.0, 1.5, 1), (1.5, 2.0, 2), (2.0, 2.5, 0), (2.5, 3.0, 1)],
            id="region",
        ),
        pytest.param(
            1.5, FAR, (0.25, 0.25, 0.5), (1, 0.5, 0), 3.0,
            [(0.0, 0.25, 1), (0.25, 0.75, 2), (0.75, 1.25, 0), (1.25, 1.5, 1),
             (1.5, 2.5, 0), (2.5, 2.75, 0), (2.75, 3.0, 1)],
            id="oblique",
        ),
    ],
)
# fmt: on
def test_repeat_walk(variant_scalar_rgb, lattice, region, o, d, maxt, expected):
    """``(mint, maxt, majorant)`` of the walk over ``[0, maxt]``. A gap the
    ray starts in is skipped by ``begin``; float-noise slivers at cell faces
    are dropped."""
    structure = _tiled(lattice, region)
    ray = mi.Ray3f(o=mi.Point3f(o), d=mi.Vector3f(d))
    walk = walk_iter(structure, ray, 0.0, maxt)
    actual = [(lo, hi, majorant) for lo, hi, _, majorant in walk if hi - lo > 1e-5]
    assert dr.allclose(actual, expected, atol=1e-6), actual


def test_repeat_lanes_independent(variants_vec_backends_once_rgb):
    """One ``next`` call where lane 0 steps inside a tile and lane 1 re-enters
    a cell after a gap: neither branch may touch the other lane."""
    structure = _tiled(1.5)
    ray = mi.Ray3f(o=mi.Point3f([0.25, 1.25], 0.5, 0.5), d=mi.Vector3f(1, 0, 0))
    state = structure.begin(ray, mi.Float(0.0), mi.Float(3.0))
    segment, state = structure.next(state, ray)

    # Lane 0 crosses the inner cell at x = 0.5; lane 1 enters the tile at
    # x = 1.5 with an empty segment.
    assert dr.allclose(segment.mint, [0.0, 0.25])
    assert dr.allclose(segment.maxt, [0.25, 0.25])
    assert dr.allclose(segment.majorant(), [1.0, 0.0])
    assert dr.allclose(state.mint, [0.25, 0.25])
    assert dr.allclose(state.maxt, [3.0, 3.0])
