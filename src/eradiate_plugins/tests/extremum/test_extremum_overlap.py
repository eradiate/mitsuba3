import mitsuba as mi
import numpy as np
import pytest


def _make_grid_volume(value, box_min, box_max, n=2):
    # Constant-valued volume occupying [box_min, box_max] in world space.
    extents = np.subtract(box_max, box_min)
    return mi.load_dict(
        {
            "type": "gridvolume",
            "grid": mi.VolumeGrid(np.full((n, n, n), value, dtype=np.float32)),
            "filter_type": "nearest",
            "accel": False,
            "to_world": mi.ScalarTransform4f()
            .translate(list(box_min))
            .scale(extents.tolist()),
        }
    )


def _make_child(value, box_min, box_max):
    volume = _make_grid_volume(value, box_min, box_max)
    structure = mi.load_dict({"type": "extremum_global"})
    structure.update_extremum(volume.bbox(), volume)
    return structure


def test_construct_no_child(variant_scalar_mono):
    with pytest.raises(RuntimeError):
        mi.load_dict({"type": "extremum_overlap"})


def test_construct_one_child(variant_scalar_mono):
    child = _make_child(1.0, [0, 0, 0], [1, 1, 1])
    overlap = mi.load_dict({"type": "extremum_overlap", "s0": child})
    assert np.allclose(overlap.bbox().min, [0, 0, 0])
    assert np.allclose(overlap.bbox().max, [1, 1, 1])


def test_construct_two_children(variant_scalar_mono):
    a = _make_child(1.0, [0, 0, 0], [1, 1, 1])
    b = _make_child(2.0, [2, 0, 0], [3, 1, 1])
    overlap = mi.load_dict({"type": "extremum_overlap", "s0": a, "s1": b})
    assert np.allclose(overlap.bbox().min, [0, 0, 0])
    assert np.allclose(overlap.bbox().max, [3, 1, 1])


def test_update_on_sigma_t_change(variant_scalar_mono):
    # extremum_overlap aggregates its children by reference: rebuilding a
    # component's own extremum structure (as heterogeneous media do on a
    # sigma_t change) must be visible through the overlap without rebuilding
    # the overlap structure itself.
    n = 2
    before = [3.0, 4.0]
    after = [5.0, 6.0]

    vol_a = _make_grid_volume(1.0, [0, 0, 0], [1, 1, 1])
    vol_b = mi.load_dict(
        {
            "type": "gridvolume",
            "grid": mi.VolumeGrid(np.array(before).reshape(n, 1, 1)),
            "filter_type": "nearest",
            "accel": False,
            "to_world": mi.ScalarTransform4f().translate([2, 0, 0]),
        }
    )

    medium = mi.load_dict(
        {
            "type": "overlap",
            "comp_a": {"type": "heterogeneous", "sigma_t": vol_a, "albedo": 0.5},
            "comp_b": {"type": "heterogeneous", "sigma_t": vol_b, "albedo": 0.5},
        }
    )

    # Ground truth: a component built directly from the "after" data,
    # independently of the update mechanism under test.
    ref_volume = mi.load_dict(
        {
            "type": "gridvolume",
            "grid": mi.VolumeGrid(np.array(after).reshape(n, 1, 1)),
            "filter_type": "nearest",
            "accel": False,
            "to_world": mi.ScalarTransform4f().translate([2, 0, 0]),
        }
    )
    ref_medium = mi.load_dict(
        {"type": "heterogeneous", "sigma_t": ref_volume, "albedo": 0.5}
    )
    it = mi.Interaction3f()
    it.p = mi.Point3f(2.5, 0.5, 0.5)
    expected = np.array(ref_medium.extremum_structure().eval_1(it))

    params = mi.eradiate.traverse(medium)
    params["component1.sigma_t.data"] = mi.TensorXf(np.array(after).reshape(n, 1, 1))
    params.update()

    got = np.array(medium.extremum_structure().eval_1(it))
    assert np.allclose(got, expected)


def test_traverse_extremum(variant_scalar_mono):
    # Two disjoint children: majorant 1 over [0, 1], majorant 2 over [2, 3],
    # separated by an empty gap. target_ot is chosen to exhaust the first
    # child's optical thickness (1 * 1 = 1) and part of the gap, then sample
    # inside the second child.
    a = _make_child(1.0, [0, 0, 0], [1, 1, 1])
    b = _make_child(2.0, [2, 0, 0], [3, 1, 1])
    overlap = mi.load_dict({"type": "extremum_overlap", "s0": a, "s1": b})

    ray = mi.Ray3f(o=[-1, 0.5, 0.5], d=[1, 0, 0])
    distance, leftover_ot = overlap.sample_test(ray, 0.0, 4.0, target_ot=1.5)

    assert np.allclose(distance, 3.25)
    assert np.allclose(leftover_ot, 0.5)


def test_next_segment(variant_scalar_mono):
    # Same disjoint setup as above. next_segment tiles: gap, child a, gap,
    # child b, then vacuum to infinity, with summed (here: single-child)
    # bounds in each occupied segment.
    a = _make_child(1.0, [0, 0, 0], [1, 1, 1])
    b = _make_child(2.0, [2, 0, 0], [3, 1, 1])
    overlap = mi.load_dict({"type": "extremum_overlap", "s0": a, "s1": b})

    ray = mi.Ray3f(o=[-1, 0.5, 0.5], d=[1, 0, 0])

    reference = [
        (0.0, 1.0, 0.0, 0.0),
        (1.0, 2.0, 1.0, 1.0),
        (2.0, 3.0, 0.0, 0.0),
        (3.0, 4.0, 2.0, 2.0),
    ]

    t = 0.0
    for mint, maxt, minorant, majorant in reference:
        segment = overlap.next_segment(ray, t)
        assert np.allclose(segment.mint, mint, atol=1e-5)
        assert np.allclose(segment.maxt, maxt, atol=1e-5)
        assert np.allclose(segment.minorant(), minorant, atol=1e-5)
        assert np.allclose(segment.majorant(), majorant, atol=1e-5)
        t = segment.maxt

    last = overlap.next_segment(ray, t)
    assert np.isinf(last.maxt)
    assert np.allclose(last.minorant(), 0.0)
    assert np.allclose(last.majorant(), 0.0)
