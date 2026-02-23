import mitsuba as mi
import numpy as np
import pytest


@pytest.fixture
def shapegroup():
    return mi.load_dict({"type": "shapegroup", "cube": {"type": "cube"}})


def transforms(n_inst):
    patch_scale = 10
    patch_offset = 5
    translate = np.random.rand(n_inst, 3)
    translate *= np.array([[patch_offset, patch_offset, 0]])
    translate -= np.array([patch_scale, patch_scale, 0])
    to_worlds = np.expand_dims(np.identity(4), 0)
    to_worlds = np.tile(to_worlds, (n_inst, 1, 1))
    to_worlds[:, 0:3, 3] = translate
    return to_worlds


def test01_create(variant_scalar_rgb, shapegroup):
    n_inst = 10
    to_worlds = transforms(n_inst)
    s = mi.load_dict(
        {
            "type": "instancelist",
            "shapegroup": shapegroup,
            "transforms": mi.TensorXf(to_worlds),
        }
    )

    assert s is not None
    assert len(s) == n_inst
