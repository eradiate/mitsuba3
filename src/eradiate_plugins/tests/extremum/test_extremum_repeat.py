import mitsuba as mi
import numpy as np


def _make_x_volume(values, n, to_world=None):
    # Radial/X resolution must be the grid's fastest (last) axis.
    data = np.array(values, dtype=float).reshape(1, 1, n)
    d = {
        "type": "gridvolume",
        "grid": mi.VolumeGrid(data),
        "filter_type": "nearest",
        "accel": False,
    }
    if to_world is not None:
        d["to_world"] = to_world
    return mi.load_dict(d)


def _make_inner(values, n, to_world=None):
    volume = _make_x_volume(values, n, to_world=to_world)
    structure = mi.load_dict(
        {"type": "extremum_grid", "resolution": mi.ScalarVector3i(n, 1, 1)}
    )
    structure.update_extremum(volume.bbox(), volume)
    return structure


def test_update_on_sigma_t_change(variant_scalar_mono):
    # extremum_repeat wraps its inner structure by reference: rebuilding the
    # wrapped medium's own extremum structure (as heterogeneous media do on
    # a sigma_t change) must be visible through the repeat without rebuilding
    # the repeat structure itself.
    n = 2
    before = [1.0, 2.0]
    after = [5.0, 6.0]

    volume = _make_x_volume(before, n)
    medium = mi.load_dict(
        {
            "type": "repeat",
            "inner": {"type": "heterogeneous", "sigma_t": volume, "albedo": 0.5},
        }
    )

    # Ground truth: a medium built directly from the "after" data,
    # independently of the update mechanism under test.
    ref_volume = _make_x_volume(after, n)
    ref_medium = mi.load_dict(
        {"type": "heterogeneous", "sigma_t": ref_volume, "albedo": 0.5}
    )
    it = mi.Interaction3f()
    it.p = mi.Point3f(0.25, 0.5, 0.5)
    expected = np.array(ref_medium.extremum_structure().eval_1(it))

    params = mi.eradiate.traverse(medium)
    params["inner.sigma_t.data"] = mi.TensorXf(np.array(after).reshape(n, 1, 1))
    params.update()

    got = np.array(medium.extremum_structure().eval_1(it))
    assert np.allclose(got, expected)


def test_traverse_extremum_through_boundary(variant_scalar_mono):
    # 2 cells of width 0.5 along x, values [1, 2], tiled with the default
    # lattice (the inner structure's own [0, 1] domain extents). target_ot
    # is chosen to exhaust the first tile's optical thickness
    # (0.5*1 + 0.5*2 = 1.5) and land inside the second tile's first cell.
    inner = _make_inner([1.0, 2.0], 2)
    repeat = mi.load_dict({"type": "extremum_repeat", "structure": inner})

    ray = mi.Ray3f(o=[0, 0.5, 0.5], d=[1, 0, 0])
    distance, leftover_ot = repeat.sample_test(ray, 0.0, 4.0, target_ot=1.7)

    assert np.allclose(distance, 1.2)
    assert np.allclose(leftover_ot, 0.2)


def test_next_segment_within_aabb(variant_scalar_mono):
    # Same tiling as above, restricted to a finite [0, 2] aabb. Querying
    # from inside the second tile's first cell returns that cell's segment.
    inner = _make_inner([1.0, 2.0], 2)
    repeat = mi.load_dict(
        {
            "type": "extremum_repeat",
            "structure": inner,
            "aabb_min": [0, 0, 0],
            "aabb_max": [2, 1, 1],
        }
    )

    ray = mi.Ray3f(o=[0, 0.5, 0.5], d=[1, 0, 0])
    segment = repeat.next_segment(ray, 1.2)

    assert np.allclose(segment.mint, 1.2)
    assert np.allclose(segment.maxt, 1.5)
    assert np.allclose(segment.minorant(), 1.0)
    assert np.allclose(segment.majorant(), 1.0)


def test_next_segment_outside_aabb(variant_scalar_mono):
    # Same finite [0, 2] aabb, ray entirely past its far edge and heading
    # further away: the whole ray from t is reported as one empty segment.
    inner = _make_inner([1.0, 2.0], 2)
    repeat = mi.load_dict(
        {
            "type": "extremum_repeat",
            "structure": inner,
            "aabb_min": [0, 0, 0],
            "aabb_max": [2, 1, 1],
        }
    )

    ray = mi.Ray3f(o=[5, 0.5, 0.5], d=[1, 0, 0])
    segment = repeat.next_segment(ray, 0.0)

    assert np.isinf(segment.maxt)
    assert np.allclose(segment.minorant(), 0.0)
    assert np.allclose(segment.majorant(), 0.0)


def test_next_segment_rotated_inner(variant_scalar_mono):
    # Inner volume rotated 90 degrees about z: its local x-variation runs
    # along world y, and so does the tiling. 4 cells of width 0.25, values
    # [1, 2, 3, 4]; the query starts in the last cell of the first tile and
    # crosses the repeat boundary at world y=1 into the folded first cell.
    to_world = mi.ScalarAffineTransform4f.rotate([0, 0, 1], 90)
    inner = _make_inner([1.0, 2.0, 3.0, 4.0], 4, to_world=to_world)
    repeat = mi.load_dict({"type": "extremum_repeat", "structure": inner})

    ray = mi.Ray3f(o=[-0.5, 0, 0.5], d=[0, 1, 0])

    last_of_first_tile = repeat.next_segment(ray, 0.75)
    assert np.allclose(last_of_first_tile.mint, 0.75)
    assert np.allclose(last_of_first_tile.maxt, 1.0)
    assert np.allclose(last_of_first_tile.minorant(), 4.0)
    assert np.allclose(last_of_first_tile.majorant(), 4.0)

    first_of_second_tile = repeat.next_segment(ray, last_of_first_tile.maxt)
    assert np.allclose(first_of_second_tile.mint, 1.0)
    assert np.allclose(first_of_second_tile.maxt, 1.25)
    assert np.allclose(first_of_second_tile.minorant(), 1.0)
    assert np.allclose(first_of_second_tile.majorant(), 1.0)
