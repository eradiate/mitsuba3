import drjit as dr
import mitsuba as mi
import numpy as np

data = np.expand_dims([[0, 1], [2, 3]], axis=-1)
colors = {"red": (1, 0, 0), "green": (0, 1, 0), "blue": (0, 0, 1), "white": (1, 1, 1)}
bsdf = {
    "type": "selectbsdf",
    "id": "mat_select",
    "indices": {
        "type": "bitmap",
        "data": data,
        "raw": True,
        "filter_type": "nearest",
        "wrap_mode": "clamp",
    },
    **{
        f"bsdf_{i}": {"type": "diffuse", "reflectance": {"type": "rgb", "value": color}}
        for i, color in enumerate(colors.values())
    },
}


def test01_create(variant_scalar_rgb):
    b = mi.load_dict(bsdf)
    assert b is not None
    assert b.component_count() == 4
    assert b.flags(0) == mi.BSDFFlags.DiffuseReflection | mi.BSDFFlags.FrontSide
    assert b.flags() == b.flags(0)


def test02_eval_scalar(variant_scalar_rgb):
    b = mi.load_dict(bsdf)
    si = dr.zeros(mi.SurfaceInteraction3f)
    si.t = 0.0
    si.wi = [0, 0, 1]
    si.p = [0, 0, 0]
    wo = mi.ScalarVector3f([0, 0, 1])
    ctx = mi.BSDFContext()

    values = []

    for uv in [[0, 0], [1, 0], [0, 1], [1, 1]]:
        si.uv = uv
        values.append(b.eval(ctx, si, wo, True))

    expected = [
        mi.Color3f(colors[name]) / dr.pi for name in ["red", "green", "blue", "white"]
    ]
    assert dr.allclose(values, expected)


def test03_eval_vec(variants_vec_backends_once_rgb):
    b = mi.load_dict(bsdf)
    si = dr.zeros(mi.SurfaceInteraction3f)
    si.t = [0.0] * 4
    si.wi = mi.Vector3f([[0] * 4, [0] * 4, [1] * 4])
    si.p = mi.Point3f([[0] * 4, [0] * 4, [0] * 4])
    si.uv = mi.Point2f([[0, 1, 0, 1], [0, 0, 1, 1]])
    wo = mi.Vector3f([[0] * 4, [0] * 4, [1] * 4])
    ctx = mi.BSDFContext()

    values = b.eval(ctx, si, wo, True)

    expected = (
        mi.Color3f(
            *[
                [colors[name][i] for name in ["red", "green", "blue", "white"]]
                for i in range(3)
            ]
        )
        / dr.pi
    )
    assert dr.allclose(values, expected)
