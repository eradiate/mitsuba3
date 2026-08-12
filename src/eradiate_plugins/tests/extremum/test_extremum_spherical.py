import mitsuba as mi
import numpy as np
import pytest


def generate_extremum_spherical(
    volume_grid,
    extremum_res,
    filter_type,
    rmin=0.5,
    rmax=1.0,
    fillmin=1.0,
    fillmax=0.0,
    transform=None,
):
    if transform is None:
        transform = mi.ScalarAffineTransform4f()

    volume = mi.load_dict(
        {
            "type": "sphericalcoordsvolume",
            "volume": {
                "type": "gridvolume",
                "grid": volume_grid,
                "filter_type": filter_type,
                "accel": False,
            },
            "rmin": rmin,
            "rmax": rmax,
            "fillmin": fillmin,
            "fillmax": fillmax,
            "to_world": transform,
        }
    )
    extremum_struct = mi.load_dict(
        {"type": "extremum_spherical", "resolution": extremum_res}
    )
    extremum_struct.update_extremum(volume.bbox(), volume)

    extremum_grid = mi.traverse(extremum_struct)["data"].numpy()
    extremum_grid = extremum_grid.reshape(
        extremum_res.z, extremum_res.y, extremum_res.x, 2
    )
    extremum_grid = extremum_grid.transpose(2, 1, 0, 3)
    return extremum_struct, extremum_grid


def test_radial_build_high_res(variant_scalar_rgb):

    n_x = 4
    n_y = 1
    n_z = 1
    n_prod = n_x * n_y * n_z
    data = np.linspace(1, n_prod, n_prod).reshape(n_x, n_y, n_z)
    volume_grid = mi.VolumeGrid(data.transpose(2, 1, 0))

    extremum_resolution = mi.ScalarVector3i(n_x, n_y, n_z)
    _, extremum_grid = generate_extremum_spherical(
        volume_grid, extremum_resolution, "nearest"
    )

    assert np.allclose(data, extremum_grid[:, :, :, 0])
    assert np.allclose(data, extremum_grid[:, :, :, 1])


def test_radial_build_half_res(variant_scalar_rgb):

    n_x = 4
    n_y = 1
    n_z = 1
    n_prod = n_x * n_y * n_z
    data = np.linspace(1, n_prod, n_prod).reshape(n_x, n_y, n_z)
    volume_grid = mi.VolumeGrid(data.transpose(2, 1, 0))

    extremum_resolution = mi.ScalarVector3i(2, 1, 1)
    _, extremum_grid = generate_extremum_spherical(
        volume_grid, extremum_resolution, "nearest"
    )

    assert np.allclose(data[::2, ::2, ::2], extremum_grid[:, :, :, 0])
    assert np.allclose(data[1::2, 1::2, 1::2], extremum_grid[:, :, :, 1])


def _make_spherical_volume(values, n, rmin, rmax, fillmin=1.0, fillmax=0.0):
    # Radial resolution must be the grid's fastest (last) axis.
    data = np.array(values, dtype=float).reshape(1, 1, n)
    grid = mi.load_dict(
        {
            "type": "gridvolume",
            "grid": mi.VolumeGrid(data),
            "filter_type": "nearest",
            "accel": False,
        }
    )
    return mi.load_dict(
        {
            "type": "sphericalcoordsvolume",
            "volume": grid,
            "rmin": rmin,
            "rmax": rmax,
            "fillmin": fillmin,
            "fillmax": fillmax,
        }
    )


def _make_extremum(volume, n):
    extremum = mi.load_dict(
        {"type": "extremum_spherical", "resolution": mi.ScalarVector3i(n, 1, 1)}
    )
    extremum.update_extremum(volume.bbox(), volume)
    return extremum


def _make_medium(medium_type, volume, extremum):
    return mi.load_dict(
        {
            "type": medium_type,
            "sigma_t": volume,
            "albedo": 0.5,
            "extremum": extremum,
        }
    )


@pytest.mark.parametrize("medium_type", ["heterogeneous", "eoheterogeneous"])
def test_update_on_sigma_t_change(variant_scalar_mono, medium_type):
    n = 4
    rmin, rmax = 0.5, 1.0
    resolution = mi.ScalarVector3i(n, 1, 1)
    before = [1.0, 2.0, 3.0, 4.0]
    after = [5.0, 6.0, 7.0, 8.0]

    volume = _make_spherical_volume(before, n, rmin, rmax)
    extremum = mi.load_dict({"type": "extremum_spherical", "resolution": resolution})
    medium = _make_medium(medium_type, volume, extremum)

    # Ground truth: an extremum structure built directly from the "after"
    # data, independently of the update mechanism under test.
    ref_volume = _make_spherical_volume(after, n, rmin, rmax)
    ref_extremum = mi.load_dict(
        {"type": "extremum_spherical", "resolution": resolution}
    )
    ref_extremum.update_extremum(ref_volume.bbox(), ref_volume)
    expected = np.array(mi.traverse(ref_extremum)["data"])

    params = mi.eradiate.traverse(medium)
    params["sigma_t.volume.data"] = mi.TensorXf(np.array(after).reshape(n, 1, 1))
    params.update()

    got = np.array(mi.traverse(medium.extremum_structure())["data"])
    assert np.allclose(got, expected)


def test_sample_multi_layer(variant_scalar_mono):
    # Radial ray entering exactly at rmax, heading to the center. 4 shells of
    # width 0.2, values [1, 2, 3, 4] from rmin=0.2 outward. Sampled inside the
    # 3rd layer (COT per layer: 0.8, 1.4, 1.8) at mid-distance.
    volume = _make_spherical_volume([1.0, 2.0, 3.0, 4.0], 4, rmin=0.2, rmax=1.0)
    extremum = _make_extremum(volume, 4)

    ray = mi.Ray3f(o=[1, 0, 0], d=[-1, 0, 0])
    distance, leftover_ot = extremum.sample_test(ray, 0.0, 2.0, target_ot=1.6)

    assert np.allclose(distance, 0.5)
    assert np.allclose(leftover_ot, 0.2)


def test_sample_escapes(variant_scalar_mono):
    # Same setup as `test_sample_mutli_layer`, but the ray escapes the structure.
    volume = _make_spherical_volume([1.0, 2.0, 3.0, 4.0], 4, rmin=0.2, rmax=1.0)
    extremum = _make_extremum(volume, 4)

    ray = mi.Ray3f(o=[1, 0, 0], d=[-1, 0, 0])
    distance, leftover_ot = extremum.sample_test(ray, 0.0, 3.0, target_ot=6.0)

    assert np.isinf(distance)
    assert np.allclose(leftover_ot, 1.6)


def test_sample_tangent_shell(variant_scalar_mono):
    # Single shell [0.6, 1.0] (value 2.0). Ray tangent to rmin. Samples past
    # tangent point.
    volume = _make_spherical_volume([2.0], 1, rmin=0.6, rmax=1.0)
    extremum = _make_extremum(volume, 1)

    ray = mi.Ray3f(o=[0.6, 0, -3], d=[0, 0, 1])
    distance, leftover_ot = extremum.sample_test(ray, 0.0, 10.0, target_ot=2.4)

    assert np.allclose(distance, 3.4, rtol=0.001)
    assert np.allclose(leftover_ot, 0.8)


def test_sample_through_origin(variant_scalar_mono):
    # Ray straight through the origin is tangent to the degenerate r=0
    # boundary. Sampled past the origin.
    volume = _make_spherical_volume([5.0], 1, rmin=0.0, rmax=1.0)
    extremum = _make_extremum(volume, 1)

    ray = mi.Ray3f(o=[0, 0, -3], d=[0, 0, 1])
    distance, leftover_ot = extremum.sample_test(ray, 0.0, 10.0, target_ot=7.5)

    assert np.allclose(distance, 3.5)
    assert np.allclose(leftover_ot, 7.5)


def test_sample_fillmax_sampled(variant_scalar_mono):
    # fillmax=0.5: the region outside rmax is not vacuum, so an interaction
    # can be sampled there before the ray ever reaches the shell.
    volume = _make_spherical_volume([2.0], 1, rmin=0.5, rmax=1.0, fillmax=0.5)
    extremum = _make_extremum(volume, 1)

    ray = mi.Ray3f(o=[3, 0, 0], d=[-1, 0, 0])
    distance, leftover_ot = extremum.sample_test(ray, 0.0, 10.0, target_ot=0.6)

    assert np.allclose(distance, 1.2)
    assert np.allclose(leftover_ot, 0.6)


def test_next_segment_miss(variant_scalar_mono):
    # Ray entirely outside the [-1, 1]^3 domain bbox, heading away from it:
    # the whole ray from t is reported as one empty segment.
    volume = _make_spherical_volume([2.0], 1, rmin=0.2, rmax=0.8)
    extremum = _make_extremum(volume, 1)

    ray = mi.Ray3f(o=[5, 5, 5], d=[1, 0, 0])
    segment = extremum.next_segment(ray, 0.0)

    assert np.isinf(segment.maxt)
    assert np.allclose(segment.minorant(), 0.0)
    assert np.allclose(segment.majorant(), 0.0)


def test_next_segment_starts_inside_shell(variant_scalar_mono):
    # 4 shells of width 0.2 from rmin=0.2 outward, values [1, 2, 3, 4].
    # Starting at r=0.9 lands in the 4th (outermost) shell, exiting at its
    # inner boundary r=0.8.
    volume = _make_spherical_volume([1.0, 2.0, 3.0, 4.0], 4, rmin=0.2, rmax=1.0)
    extremum = _make_extremum(volume, 4)

    ray = mi.Ray3f(o=[0.9, 0, 0], d=[-1, 0, 0])
    segment = extremum.next_segment(ray, 0.0)

    assert np.allclose(segment.mint, 0.0)
    assert np.allclose(segment.maxt, 0.1)
    assert np.allclose(segment.minorant(), 4.0)
    assert np.allclose(segment.majorant(), 4.0)


def test_next_segment_starts_outside_rmax(variant_scalar_mono):
    # Single shell [0.2, 0.8] (value 2.0), fillmax=0.5. Starting at r=0.95,
    # inside the domain bbox but outside rmax, lands in the fillmax region
    # and exits at the rmax boundary r=0.8.
    volume = _make_spherical_volume([2.0], 1, rmin=0.2, rmax=0.8, fillmax=0.5)
    extremum = _make_extremum(volume, 1)

    ray = mi.Ray3f(o=[0.95, 0, 0], d=[-1, 0, 0])
    segment = extremum.next_segment(ray, 0.0)

    assert np.allclose(segment.mint, 0.0)
    assert np.allclose(segment.maxt, 0.15)
    assert np.allclose(segment.minorant(), 0.5)
    assert np.allclose(segment.majorant(), 0.5)
