import drjit as dr
import mitsuba as mi
import numpy as np
import pytest


def test_construct(variant_scalar_rgb, tmp_path):
    # Construct without parameters
    volume = mi.load_dict({"type": "sphericalcoordsvolume"})
    assert volume is not None

    # Construct with rmin > rmax:
    with pytest.raises(Exception):
        mi.load_dict({"type": "sphericalcoordsvolume", "rmin": 0.8, "rmax": 0.2})

    # Construct with all parameters set to typical values
    volume = mi.load_dict(
        {
            "type": "sphericalcoordsvolume",
            "rmin": 0.0,
            "rmax": 1.0,
            "fillmin": 0.0,
            "fillmax": 0.0,
            "volume": {"type": "gridvolume", "grid": mi.VolumeGrid([[[1.0]]])},
        }
    )
    assert volume is not None


@pytest.mark.parametrize(
    "point,result",
    [
        ([0.9, 0.0, 0.0], 1.0),
        ([0.0, 0.0, 0.1], 3.0),
        ([0.0, 0.5, 0.0], 2.0),
    ],
)
def test_eval_basic(variant_scalar_rgb, tmp_path, point, result):
    volume = mi.load_dict(
        {
            "type": "sphericalcoordsvolume",
            "rmin": 0.2,
            "rmax": 0.8,
            "fillmin": 3.0,
            "fillmax": 1.0,
            "volume": {"type": "gridvolume", "grid": mi.VolumeGrid([[[2.0]]])},
        }
    )

    it = mi.Interaction3f()
    it.p = point

    assert dr.all(volume.eval(it) == result)


def test_eval_advanced(variant_scalar_rgb, tmp_path):
    data = np.broadcast_to(
        np.arange(1, 7, 1, dtype=float).reshape((-1, 1, 1, 1)), (6, 2, 2, 3)
    )
    volume = mi.load_dict(
        {
            "type": "sphericalcoordsvolume",
            "rmin": 0.0,
            "rmax": 1.0,
            "fillmin": 3.0,
            "fillmax": 1.0,
            "volume": {
                "type": "gridvolume",
                "grid": mi.VolumeGrid(data),
                "filter_type": "nearest",
            },
        }
    )

    # Have one point in each sector of the spherical volume.
    # Since the atan2 method maps [-pi, pi] to [0, 1], we choose the following
    # order of inputs.
    phis = np.array([-175, -115, -55, 5, 65, 125])  # degree
    results = np.empty((len(phis), 3))

    for i, phi in enumerate(phis):
        p = [0.5 * np.cos(np.deg2rad(phi)), 0.5 * np.sin(np.deg2rad(phi)), 0.0]
        it = mi.Interaction3f()
        it.p = p
        results[i, :] = volume.eval(it)

    expected = np.broadcast_to(
        np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]).reshape((-1, 1)),
        (6, 3),
    )
    assert np.all(results == expected)


def spherical_volume(data, fillmin=0.0, fillmax=0.0, rmin=0.0, rmax=1.0):
    return mi.load_dict(
        {
            "type": "sphericalcoordsvolume",
            "rmin": rmin,
            "rmax": rmax,
            "fillmin": fillmin,
            "fillmax": fillmax,
            "volume": {
                "type": "gridvolume",
                "grid": mi.VolumeGrid(data),
                "filter_type": "nearest",
            },
        }
    )


def test_min_max(variant_scalar_rgb):
    volume = spherical_volume(
        data=np.arange(1.0, 28.0).reshape(3, 3, 3),
        fillmin=10.0,
        fillmax=10.0,
        rmin=0.5,
        rmax=1.0,
    )

    # fillmin/fillmax minimum
    assert np.allclose(volume.min(), 1.0)
    assert np.allclose(volume.max(), 27.0)


def test_min_max_fill(variant_scalar_rgb):
    volume = spherical_volume(
        data=np.arange(1.0, 28.0).reshape(3, 3, 3),
        fillmin=0.0,
        fillmax=100.0,
        rmin=0.5,
        rmax=1.0,
    )

    # fillmin/fillmax minimum
    assert np.allclose(volume.min(), 0.0)
    assert np.allclose(volume.max(), 100.0)


def test_extremum(variant_scalar_rgb):
    volume = spherical_volume(
        data=np.arange(1.0, 28.0).reshape(3, 3, 3),
        fillmin=0.0,
        fillmax=100.0,
        rmin=0.5,
        rmax=1.0,
    )

    bbox = mi.BoundingBox3f(min=[-0.5, -0.5, -0.5], max=[0.5, 0.5, 0.5])
    minorant, majorant = volume.extremum(bbox)

    # evaluates to global minorant and majorant.
    assert np.allclose(minorant, 0.0)
    assert np.allclose(majorant, 100.0)


def test_extremum_local(variant_scalar_rgb):
    data = np.broadcast_to(np.arange(1.0, 5.0).reshape(-1, 1, 1), (4, 4, 4))
    print(data)
    volume = spherical_volume(
        data=data,
        fillmin=0.0,
        fillmax=1.0,
        # does not affect local extremum calculations
        rmin=0.8,
        rmax=0.9,
    )

    bbox = mi.BoundingBox3f(min=[0.0, 0, 0], max=[0.24, 0.24, 0.24])
    minorant, majorant = volume.extremum_local(bbox)

    # fillmin/fillmax minimum
    assert np.allclose(minorant, 1.0)
    assert np.allclose(majorant, 1.0)

    bbox = mi.BoundingBox3f(min=[0.0, 0, 0], max=[0.49, 0.49, 0.49])
    minorant, majorant = volume.extremum_local(bbox)

    # fillmin/fillmax minimum
    assert np.allclose(minorant, 1.0)
    assert np.allclose(majorant, 2.0)
