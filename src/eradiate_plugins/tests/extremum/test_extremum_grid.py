import mitsuba as mi
import numpy as np
import pytest


def generate_extremum_grid(
    volume_grid,
    extremum_res,
    filter_type,
    transform=None,
):
    if transform is None:
        transform = mi.ScalarAffineTransform4f()

    volume = mi.load_dict(
        {
            "type": "gridvolume",
            "grid": volume_grid,
            "filter_type": filter_type,
            "accel": False,
            "to_world": transform,
        }
    )
    extremum_struct = mi.load_dict(
        {"type": "extremum_grid", "volume": volume, "resolution": extremum_res}
    )

    extremum_grid = mi.traverse(extremum_struct)["data"].numpy()
    extremum_grid = extremum_grid.reshape(
        extremum_res.z, extremum_res.y, extremum_res.x, 2
    )
    extremum_grid = extremum_grid.transpose(2, 1, 0, 3)
    return extremum_struct, extremum_grid


def test_build_high_res(variant_scalar_rgb):

    n_x = 4
    n_y = 8
    n_z = 3
    n_prod = n_x * n_y * n_z
    data = np.linspace(1, n_prod, n_prod).reshape(n_x, n_y, n_z)
    volume_grid = mi.VolumeGrid(data.transpose(2, 1, 0))

    extremum_resolution = mi.ScalarVector3i(n_x, n_y, n_z)
    _, extremum_grid = generate_extremum_grid(
        volume_grid, extremum_resolution, "nearest"
    )

    assert np.allclose(data, extremum_grid[:, :, :, 0])
    assert np.allclose(data, extremum_grid[:, :, :, 1])


def test_build_half_res(variant_scalar_rgb):

    n_x = 4
    n_y = 8
    n_z = 4
    n_prod = n_x * n_y * n_z
    data = np.linspace(1, n_prod, n_prod).reshape(n_x, n_y, n_z)
    volume_grid = mi.VolumeGrid(data.transpose(2, 1, 0))

    extremum_resolution = mi.ScalarVector3i(2, 4, 2)
    _, extremum_grid = generate_extremum_grid(
        volume_grid, extremum_resolution, "nearest"
    )

    assert np.allclose(data[::2, ::2, ::2], extremum_grid[:, :, :, 0])
    assert np.allclose(data[1::2, 1::2, 1::2], extremum_grid[:, :, :, 1])


def test_build_not_multiple(variant_scalar_rgb):
    n_x = 4
    n_y = 9
    n_z = 1
    n_prod = n_x * n_y * n_z
    data = np.linspace(1, n_prod, n_prod).reshape(n_x, n_y, n_z)
    volume_grid = mi.VolumeGrid(data.transpose(2, 1, 0))

    extremum_resolution = mi.ScalarVector3i(3, 4, 1)
    _, extremum_grid = generate_extremum_grid(
        volume_grid, extremum_resolution, "nearest"
    )

    minorant_reference = np.array(
        [1.0, 3.0, 5.0, 7.0, 10.0, 12.0, 14.0, 16.0, 19.0, 21.0, 23.0, 25.0]
    ).reshape(3, 4, 1)
    majorant_reference = np.array(
        [12.0, 14.0, 16.0, 18.0, 21.0, 23.0, 25.0, 27.0, 30.0, 32.0, 34.0, 36.0]
    ).reshape(3, 4, 1)

    assert np.allclose(minorant_reference, extremum_grid[:, :, :, 0])
    assert np.allclose(majorant_reference, extremum_grid[:, :, :, 1])


def test_build_trilinear(variant_scalar_rgb):
    n_x = 3
    n_y = 6
    n_z = 1
    n_prod = n_x * n_y * n_z
    data = np.linspace(1, n_prod, n_prod).reshape(n_x, n_y, n_z)
    volume_grid = mi.VolumeGrid(data.transpose(2, 1, 0))

    extremum_resolution = mi.ScalarVector3i(3, 6, 1)
    _, extremum_grid = generate_extremum_grid(
        volume_grid, extremum_resolution, "trilinear"
    )

    minorant_reference = np.array(
        [
            1.0,
            1.0,
            2.0,
            3.0,
            4.0,
            5.0,
            1.0,
            1.0,
            2.0,
            3.0,
            4.0,
            5.0,
            7.0,
            7.0,
            8.0,
            9.0,
            10.0,
            11.0,
        ]
    ).reshape(3, 6, 1)
    majorant_reference = np.array(
        [
            8.0,
            9.0,
            10.0,
            11.0,
            12.0,
            12.0,
            14.0,
            15.0,
            16.0,
            17.0,
            18.0,
            18.0,
            14.0,
            15.0,
            16.0,
            17.0,
            18.0,
            18.0,
        ]
    ).reshape(3, 6, 1)

    assert np.allclose(minorant_reference, extremum_grid[:, :, :, 0])
    assert np.allclose(majorant_reference, extremum_grid[:, :, :, 1])


def test_build_rotated(variant_scalar_rgb):
    n_x = 3
    n_y = 4
    n_z = 1
    n_prod = n_x * n_y * n_z
    data = np.linspace(1, n_prod, n_prod).reshape(n_x, n_y, n_z)
    volume_grid = mi.VolumeGrid(data.transpose(2, 1, 0))

    extremum_resolution = mi.ScalarVector3i(3, 4, 1)
    _, extremum_grid = generate_extremum_grid(
        volume_grid,
        extremum_resolution,
        "nearest",
        transform=mi.ScalarAffineTransform4f.rotate([0, 0, 1], 180),
    )

    assert np.allclose(data.squeeze(), extremum_grid[:, :, :, 1].squeeze())


def test_build_scaled(variant_scalar_rgb):
    pass


def _make_grid_volume(values, n):
    data = np.array(values, dtype=float).reshape(n, 1, 1)
    return mi.load_dict(
        {
            "type": "gridvolume",
            "grid": mi.VolumeGrid(data),
            "filter_type": "nearest",
            "accel": False,
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
    resolution = mi.ScalarVector3i(1, 1, n)
    before = [1.0, 2.0, 3.0, 4.0]
    after = [5.0, 6.0, 7.0, 8.0]

    volume = _make_grid_volume(before, n)
    extremum = mi.load_dict(
        {"type": "extremum_grid", "volume": volume, "resolution": resolution}
    )
    medium = _make_medium(medium_type, volume, extremum)

    # Ground truth: an extremum structure built directly from the "after"
    # data, independently of the update mechanism under test.
    ref_volume = _make_grid_volume(after, n)
    ref_extremum = mi.load_dict(
        {"type": "extremum_grid", "volume": ref_volume, "resolution": resolution}
    )
    expected = np.array(mi.traverse(ref_extremum)["data"])

    params = mi.eradiate.traverse(medium)
    params["sigma_t.data"] = mi.TensorXf(np.array(after).reshape(n, 1, 1))
    params.update()

    got = np.array(mi.traverse(medium.extremum_structure())["data"])
    assert np.allclose(got, expected)
