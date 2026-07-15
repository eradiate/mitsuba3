import mitsuba as mi
import numpy as np
import pytest


def test01_constant_min_max(variant_scalar_rgb):
    vol = mi.load_dict(
        {
            "type": "constvolume",
            "value": {"type": "srgb", "color": mi.Color3f([0.5, 1.0, 0.3])},
        }
    )

    assert np.allclose(vol.min(), 0.3)
    assert np.allclose(vol.max(), 1.0)


def test01_constant_extremum(variant_scalar_rgb):
    vol = mi.load_dict(
        {
            "type": "constvolume",
            "value": {"type": "srgb", "color": mi.Color3f([0.5, 1.0, 0.3])},
        }
    )

    bbox = mi.BoundingBox3f()
    min, max = vol.extremum(bbox)
    assert np.allclose(min, 0.3)
    assert np.allclose(max, 1.0)
