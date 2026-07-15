import drjit as dr
import mitsuba as mi
import numpy as np
import pytest


def test_min_max_wrap(variant_scalar_rgb):
    data = dr.full(mi.TensorXf, 1, [3, 3, 3])
    data[0, 0, 0] = 0.5
    data[1, 1, 1] = 2.0
    grid = mi.VolumeGrid(data)

    vol = mi.load_dict(
        {
            "type": "gridvolume",
            "grid": grid,
            "wrap": True,
        }
    )

    min = vol.min()
    max = vol.max()
    assert dr.allclose(min, 0.5)
    assert dr.allclose(max, 2.0)


def test_min_max_no_wrap(variant_scalar_rgb):
    data = dr.full(mi.TensorXf, 1, [3, 3, 3])
    data[0, 0, 0] = 0.5
    data[1, 1, 1] = 2.0
    grid = mi.VolumeGrid(data)

    vol = mi.load_dict(
        {
            "type": "gridvolume",
            "grid": grid,
            "wrap": False,
        }
    )

    min = vol.min()
    max = vol.max()
    assert dr.allclose(min, 0.0)
    assert dr.allclose(max, 2.0)


def make_gridvolume(data, filter_type="nearest", to_world=None, wrap_mode="clamp"):
    wrap = wrap_mode != "none"
    wrap_mode = wrap_mode if wrap_mode != "none" else "clamp"

    d = {
        "type": "gridvolume",
        "grid": mi.VolumeGrid(data.astype(np.float32)),
        "filter_type": filter_type,
        "accel": False,
        "wrap_mode": wrap_mode,
        "wrap": wrap,
    }
    if to_world is not None:
        d["to_world"] = to_world
    return mi.load_dict(d)


def _sample_points(rng, box_min, box_max, n=64):
    pts = rng.uniform(box_min, box_max, (n, 3))
    corners = np.array(
        [
            [bx, by, bz]
            for bx in (box_min[0], box_max[0])
            for by in (box_min[1], box_max[1])
            for bz in (box_min[2], box_max[2])
        ]
    )
    return np.vstack([pts, corners])


def _eval_1(volume, p):
    it = mi.Interaction3f()
    it.p = mi.Point3f(*p.tolist())
    return volume.eval_1(it)


@pytest.mark.parametrize("wrap_mode", ["clamp", "repeat", "mirror", "none"])
@pytest.mark.parametrize("filter_type", ["nearest", "trilinear"])
def test_extremum_property(variant_scalar_mono, variant_, wrap_mode, filter_type):
    rng = np.random.default_rng(hash((wrap_mode, filter_type)) % 2**32)
    data = rng.uniform(0.1, 3.0, (3, 4, 5))

    transforms = [
        None,
        mi.ScalarTransform4f()
        .translate([0.3, -0.2, 0.1])
        .rotate([1, 2, 3], 37.0)
        .scale([1.5, 0.8, 1.2]),
    ]
    for to_world in transforms:
        volume = make_gridvolume(data, filter_type, to_world, wrap_mode)
        for _ in range(10):
            c = rng.uniform(-1.0, 2.0, 3)
            h = rng.uniform(0.02, 0.8, 3)
            box_min, box_max = c - h, c + h
            mn, mx = volume.extremum(
                mi.BoundingBox3f(box_min.tolist(), box_max.tolist())
            )
            assert mn <= mx
            for p in _sample_points(rng, box_min, box_max):
                v = _eval_1(volume, p)
                assert mn <= v + 1e-5 and v - 1e-5 <= mx, (
                    f"eval {v} outside bounds ({mn}, {mx}) at {p} "
                    f"(wrap={wrap_mode}, filter={filter_type}, "
                    f"to_world={'rotated' if to_world else 'identity'})"
                )


def test_extremum_local_property(variant_scalar_mono):
    rng = np.random.default_rng(5)
    data = rng.uniform(0.1, 3.0, (3, 4, 5))
    to_world = (
        mi.ScalarTransform4f()
        .translate([0.5, 0.0, -0.3])
        .rotate([0, 0, 1], 45.0)
        .scale(2.0)
    )
    M = np.array(to_world.matrix)
    volume = make_gridvolume(data, "trilinear", to_world)

    for _ in range(10):
        c = rng.uniform(0.0, 1.0, 3)
        h = rng.uniform(0.02, 0.4, 3)
        box_min = np.clip(c - h, 0.0, 1.0)
        box_max = np.clip(c + h, 0.0, 1.0)
        mn, mx = volume.extremum_local(
            mi.BoundingBox3f(box_min.tolist(), box_max.tolist())
        )
        # Sample in the native (local) parameterization, evaluate in world
        for p_local in _sample_points(rng, box_min, box_max):
            p_world = M[:3, :3] @ p_local + M[:3, 3]
            v = _eval_1(volume, p_world)
            assert mn <= v + 1e-5 and v - 1e-5 <= mx


def test_extremum_none_wrap_outside(variant_scalar_mono):
    # `none` wrap: eval is zero outside the volume box, so a disjoint query
    # returns (0, 0) and a straddling query has minorant 0
    data = np.full((2, 2, 2), 2.0)
    volume = make_gridvolume(data, wrap_mode="none")

    mn, mx = volume.extremum(mi.BoundingBox3f([2, 2, 2], [3, 3, 3]))
    assert mn == 0.0 and mx == 0.0

    mn, mx = volume.extremum(mi.BoundingBox3f([0.5, 0.5, 0.5], [2, 2, 2]))
    assert mn == 0.0 and np.isclose(mx, 2.0)
