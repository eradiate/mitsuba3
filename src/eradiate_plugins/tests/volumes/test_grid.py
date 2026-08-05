import drjit as dr
import mitsuba as mi
import pytest


def test_wrap(variant_all_backends_once):
    grid = dr.full(mi.TensorXf, 1.0, [3, 3, 3])
    grid[0, 0, 0] = 0.0
    grid[1, 1, 1] = 0.5
    volume_grid = mi.VolumeGrid(grid)
    inside_value = (6 * 1.0 + 1 * 0.5 + 1 * 0.0) / 8  # interpolated value

    it = mi.Interaction3f()

    # wrap=False: queries outside [0, 1]^3 evaluate to zero, queries inside
    # evaluate normally
    vol = mi.load_dict({"type": "gridvolume", "grid": volume_grid, "wrap": False})

    it.p = mi.Point3f(1.0) / 3
    assert dr.allclose(vol.eval(it), inside_value)

    for p in [
        (-0.1, 0.5, 0.5),
        (1.1, 0.5, 0.5),
        (0.5, -0.1, 0.5),
        (0.5, 1.1, 0.5),
        (0.5, 0.5, -0.1),
        (0.5, 0.5, 1.1),
    ]:
        it.p = mi.Point3f(*p)
        assert dr.allclose(vol.eval(it), 0.0)

    # The [0, 1]^3 domain is inclusive: querying exactly on the boundary
    # (x=0 or x=1, away from the two special cells so clamp-mode
    # interpolation resolves to the constant 1.0 slice) must not be zeroed.
    for p in [(0.0, 0.5, 0.5), (1.0, 0.5, 0.5)]:
        it.p = mi.Point3f(*p)
        assert dr.allclose(vol.eval(it), 1.0)

    # wrap=True (default): queries outside [0, 1]^3 are handled by
    # wrap_mode instead of evaluating to zero
    vol_wrap = mi.load_dict({"type": "gridvolume", "grid": volume_grid, "wrap": True})

    it.p = mi.Point3f(1.0) / 3
    assert dr.allclose(vol_wrap.eval(it), inside_value)

    it.p = mi.Point3f(-0.1, 0.5, 0.5)
    assert not dr.allclose(vol_wrap.eval(it), 0.0)
