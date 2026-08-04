import mitsuba as mi
import numpy as np
import pytest


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


def test_build(variant_scalar_mono_double):
    volume = _make_grid_volume([1.0, 2.0, 3.0, 4.0], 4)
    extremum = mi.load_dict({"type": "extremum_global", "volume": volume})

    it = mi.Interaction3f()
    it.p = mi.Point3f(0.0, 0.0, 0.0)
    assert np.allclose(extremum.eval_1(it), (1.0, 4.0))


@pytest.mark.parametrize(
    "medium_type", ["heterogeneous", "eoheterogeneous", "piecewise", "homogeneous"]
)
def test_update_on_sigma_t_change(variant_scalar_mono_double, medium_type):
    n = 4
    before = [1.0, 2.0, 3.0, 4.0]
    after = [5.0, 6.0, 7.0, 8.0]

    volume = _make_grid_volume(before, n)
    medium = mi.load_dict({"type": medium_type, "sigma_t": volume, "albedo": 0.5})

    # Ground truth: an extremum structure built directly from the "after"
    # data, independently of the update mechanism under test.
    ref_volume = _make_grid_volume(after, n)
    ref_medium = mi.load_dict(
        {"type": medium_type, "sigma_t": ref_volume, "albedo": 0.5}
    )
    it = mi.Interaction3f()
    it.p = mi.Point3f(0.0, 0.0, 0.0)
    expected = np.array(ref_medium.extremum_structure().eval_1(it))

    params = mi.eradiate.traverse(medium)
    params["sigma_t.data"] = mi.TensorXf(np.array(after).reshape(n, 1, 1))
    params.update()

    got = np.array(medium.extremum_structure().eval_1(it))
    assert np.allclose(got, expected)
