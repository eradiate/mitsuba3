import os

import mitsuba as mi

mi.set_variant("llvm_rgb")

import numpy as np
from util import render_cache
from volpathgrid import VolPathGridIntegratorPy

OUTDIR = "."


def make_test_scene(integrator: dict = None):
    n = 10000  # Number of grid points
    zmax = 1.0  # Maximum height
    tau = 1.0  # Total optical thickness
    l = zmax / 10  # Typical cutoff height
    k = 1 / l  # Decay factor
    alpha = tau * k / (1 - np.exp(-k * zmax))  # Pre-exponential factor
    print(f"{alpha = }")

    z = np.linspace(0, zmax, n)
    atmosphere_data = np.reshape(alpha * np.exp(-k * z), (-1, 1, 1))
    atmosphere_data = np.broadcast_to(atmosphere_data, (atmosphere_data.shape[0], 3, 3))

    scene_dict_ref = {
        "type": "scene",
        "atmosphere_medium": {
            "type": "heterogeneous",
            "sigma_t": {
                "type": "gridvolume",
                "grid": mi.VolumeGrid(atmosphere_data),
                "to_world": mi.ScalarTransform4f.translate([-0.5, -0.5, -0.5]),
                "accel": False,
            },
            "albedo": 0.5,
            # "majorant_resolution_factor": [1, 1, 100],
            "phase": {"type": "isotropic"},
        },
        "atmosphere_shape": {
            "type": "cube",
            "interior": {"type": "ref", "id": "atmosphere_medium"},
            "to_world": mi.ScalarTransform4f.translate([0.0, 0.0, 0.0])
            @ mi.ScalarTransform4f.scale(0.5),
            "bsdf": {"type": "null"},
        },
        "surface": {
            "type": "rectangle",
            "to_world": mi.ScalarTransform4f.translate([0.0, 0.0, -0.5])
            @ mi.ScalarTransform4f.scale(0.5),
            "bsdf": {
                "type": "diffuse",
                "reflectance": {
                    "type": "checkerboard",
                    "to_uv": mi.ScalarTransform4f.scale(2),
                },
            },
        },
        "camera": {
            "type": "perspective",
            "to_world": mi.ScalarTransform4f.look_at(
                origin=[0.0, -2.5, 0.0],
                target=[0.0, 0.0, 0.0],
                up=[0, 0, 1],
            ),
            "sampler": {"type": "independent"},
            "film": {
                "type": "hdrfilm",
                "width": 256,
                "height": 256,
                "rfilter": {"type": "box"},
            },
            "near_clip": 1e-8,
            "far_clip": 1e8,
        },
        "illumination": {
            "type": "directional",
            "direction": [0, -1, -1],
            "irradiance": 5.0,
        },
        "integrator": {"type": "path"} if integrator is None else integrator,
    }

    return scene_dict_ref


def test_01_volpathgrid_basic():
    output_dir = os.path.join(OUTDIR, "test_integrators", "test_volpathgrid_basic")
    os.makedirs(output_dir, exist_ok=True)
    mi.set_variant("llvm_rgb")

    scene_dict = make_test_scene()
    scene = mi.load_dict(scene_dict)
    integrators = {
        "volpath": mi.load_dict(
            {
                "type": "volpath",
                "max_depth": 64,
                "rr_depth": 999,
            }
        ),
        "volpathgridpy": mi.load_dict(
            {
                "type": "volpathgridpy",
                "max_depth": 64,
                "rr_depth": 999,
            }
        ),
    }

    results = {}
    for k, integrator in integrators.items():
        fname = os.path.join(output_dir, f"preview_{k}.exr")

        @render_cache(fname, overwrite=True)
        def fn():
            return mi.render(scene, integrator=integrator, spp=32)

        results[k] = fn()

    ref_integrator = "volpath"
    for k, image in results.items():
        if k != ref_integrator:
            a = image
            b = results[ref_integrator]
            assert np.allclose(a, b, atol=5e-2)


if __name__ == "__main__":
    mi.set_variant("llvm_rgb")
    integrator = mi.load_dict({"type": "volpathgridpy"})
