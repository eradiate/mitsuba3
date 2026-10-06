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


def _make_spherical_volume(values, n, rmin, rmax):
    data = np.array(values, dtype=float).reshape(n, 1, 1)
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
            "fillmin": 1.0,
            "fillmax": 0.0,
        }
    )


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
def test_update_on_sigma_t_change(variant_scalar_mono_double, medium_type):
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
