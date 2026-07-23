"""Render-comparison tests for overlapping (multicomponent), clamp-extended,
and repeating media: each scene is rendered with eovolpath and compared
against a reference medium whose extinction field is baked into a single
gridvolume. Images must match within Monte Carlo noise."""

import mitsuba as mi
import numpy as np
import pytest

SPP = 256
TOL = 0.02

_uid = [0]


def render(medium_dict, seed):
    scene = mi.load_dict(
        {
            "type": "scene",
            "integrator": {"type": "eovolpath", "max_depth": 64},
            "sensor": {
                "type": "perspective",
                "fov": 60,
                "to_world": mi.ScalarTransform4f().look_at(
                    origin=[0.75, 0.75, -3.0],
                    target=[0.75, 0.75, 0.75],
                    up=[0, 1, 0],
                ),
                "film": {
                    "type": "hdrfilm",
                    "width": 16,
                    "height": 16,
                    "rfilter": {"type": "box"},
                    "pixel_format": "luminance",
                },
                "sampler": {"type": "independent", "sample_count": SPP},
            },
            "emitter": {"type": "constant", "radiance": 1.0},
            "shape": {
                "type": "cube",
                "to_world": mi.ScalarTransform4f()
                .translate([0.75, 0.75, 0.75])
                .scale(0.75),
                "bsdf": {"type": "null"},
                "interior": medium_dict,
            },
        }
    )
    return np.array(mi.render(scene, seed=seed, spp=SPP))[:, :, 0]


def medium_grid(data, to_world, res=(4, 4, 4), aabb=None):
    # Unique extremum id defeats load_dict deduplication: one structure
    # belongs to one medium
    _uid[0] += 1
    d = {
        "type": "eoheterogeneous",
        "sigma_t": {
            "type": "gridvolume",
            "grid": mi.VolumeGrid(data.astype(np.float32)),
            "filter_type": "nearest",
            "to_world": to_world,
        },
        "albedo": 0.5,
        "extremum": {
            "type": "extremum_grid",
            "resolution": list(res),
            "id": f"extremum_render_{_uid[0]}",
        },
    }
    if aabb is not None:
        d["aabb_min"], d["aabb_max"] = aabb
    return d


def assert_images_match(img, ref, tol=TOL):
    rel = abs(img.mean() - ref.mean()) / ref.mean()
    assert rel < tol, f"mean {img.mean():.5f} vs {ref.mean():.5f} (rel {rel:.4f})"


def test_overlap_full(variant_scalar_mono_double):
    # Two coincident heterogeneous components vs. the voxelwise sum
    rng = np.random.default_rng(42)
    a = rng.uniform(0.05, 2.4, (8, 8, 8))
    b = rng.uniform(0.05, 2.4, (8, 8, 8))
    tw = mi.ScalarTransform4f().scale(1.5)

    overlap = {
        "type": "multicomponent",
        "a": medium_grid(a, tw),
        "b": medium_grid(b, tw),
    }
    assert_images_match(render(overlap, 1), render(medium_grid(a + b, tw), 11))


def test_overlap_partial(variant_scalar_mono_double):
    # A over x in [0, 1], B over x in [0.5, 1.5]: the disjoint parts exercise
    # empty segments in the aggregate
    rng = np.random.default_rng(43)
    a = rng.uniform(0.05, 3.0, (4, 4, 8))  # (z, y, x)
    b = rng.uniform(0.05, 3.0, (4, 4, 8))
    tw_a = mi.ScalarTransform4f().scale([1.0, 1.5, 1.5])
    tw_b = mi.ScalarTransform4f().translate([0.5, 0, 0]).scale([1.0, 1.5, 1.5])

    common = np.zeros((4, 4, 12))
    common[:, :, 0:8] += a
    common[:, :, 4:12] += b
    tw_c = mi.ScalarTransform4f().scale(1.5)

    overlap = {
        "type": "multicomponent",
        "a": medium_grid(a, tw_a),
        "b": medium_grid(b, tw_b),
    }
    assert_images_match(render(overlap, 2), render(medium_grid(common, tw_c), 12))


def test_overlap_homogeneous_child(variant_scalar_mono_double):
    # Homogeneous child with an infinite domain: the aggregate's domain is
    # infinite and traversal is bounded by the enclosing shape
    rng = np.random.default_rng(44)
    a = rng.uniform(0.05, 2.0, (8, 8, 8))
    tw = mi.ScalarTransform4f().scale(1.5)
    c_val = 0.8

    overlap = {
        "type": "multicomponent",
        "a": medium_grid(a, tw),
        "b": {"type": "homogeneous", "sigma_t": c_val, "albedo": 0.5},
    }
    assert_images_match(render(overlap, 3), render(medium_grid(a + c_val, tw), 13))


def test_clamp_extension(variant_scalar_mono_double):
    # Domain much larger than the volume box: clamp wrap extends the edge
    # values; reference bakes the extension into a bigger volume
    rng = np.random.default_rng(45)
    b = rng.uniform(0.1, 2.5, (4, 4, 4))
    tw_small = mi.ScalarTransform4f().translate([0.5, 0.5, 0.5]).scale(0.5)
    tw_full = mi.ScalarTransform4f().scale(1.5)

    clamped = medium_grid(b, tw_small, res=(6, 6, 6), aabb=([0, 0, 0], [1.5, 1.5, 1.5]))
    baked = np.pad(b, 4, mode="edge")  # 12^3 over [0, 1.5]^3, texel size 0.125
    assert_images_match(
        render(clamped, 4), render(medium_grid(baked, tw_full, res=(6, 6, 6)), 14)
    )


def test_repeat_axis_tiling(variant_scalar_mono_double):
    # 3x3x3 axis-aligned tiling of a [0, 0.5]^3 tile vs. the tiled field
    # baked into one volume
    rng = np.random.default_rng(46)
    a = rng.uniform(0.1, 3.0, (4, 4, 4))
    tw_tile = mi.ScalarTransform4f().scale(0.5)
    tw_full = mi.ScalarTransform4f().scale(1.5)

    repeat = {
        "type": "repeat",
        "inner": medium_grid(a, tw_tile),
        "aabb_min": [0, 0, 0],
        "aabb_max": [1.5, 1.5, 1.5],
    }
    reference = medium_grid(np.tile(a, (3, 3, 3)), tw_full, res=(6, 6, 6))
    assert_images_match(render(repeat, 5), render(reference, 15))


def test_repeat_overlap_composition(variant_scalar_mono_double):
    # repeat(multicomponent(a, b)): composition of the wrappers
    rng = np.random.default_rng(47)
    a = rng.uniform(0.1, 3.0, (4, 4, 4))
    b = rng.uniform(0.1, 3.0, (4, 4, 4))
    tw_tile = mi.ScalarTransform4f().scale(0.5)
    tw_full = mi.ScalarTransform4f().scale(1.5)

    repeat = {
        "type": "repeat",
        "inner": {
            "type": "multicomponent",
            "a": medium_grid(a, tw_tile),
            "b": medium_grid(b, tw_tile),
        },
        "aabb_min": [0, 0, 0],
        "aabb_max": [1.5, 1.5, 1.5],
    }
    reference = medium_grid(np.tile(a + b, (3, 3, 3)), tw_full, res=(6, 6, 6))
    assert_images_match(render(repeat, 6), render(reference, 16))


def test_extremum_rebuild_on_update(variant_scalar_mono_double):
    # Updating sigma_t through the parameter traversal must rebuild the
    # medium's extremum structure (design §4: invalidation flows through the
    # medium, there is no observer back-edge anymore)
    rng = np.random.default_rng(49)
    a = rng.uniform(0.05, 2.0, (8, 8, 8))
    tw = mi.ScalarTransform4f().scale(1.5)

    medium = mi.load_dict(medium_grid(a, tw))
    params = mi.traverse(medium)
    params["sigma_t.data"] = np.asarray(params["sigma_t.data"]) * 2.0
    params.update()

    struct_data = np.array(params["extremum_structure.data"])
    ref = mi.load_dict(medium_grid(a * 2.0, tw))
    ref_data = np.array(mi.traverse(ref)["extremum_structure.data"])
    assert np.allclose(struct_data, ref_data)


def test_repeat_tilted_lattice(variant_scalar_mono_double):
    # Field varies along z only; the x/y lattice vectors are tilted so the
    # tiling is non-axis-aligned but the tiled field is unchanged. The inner
    # domain covers the sheared cell via clamp extension.
    rng = np.random.default_rng(48)
    fz = rng.uniform(0.2, 3.0, (4, 1, 1))  # (z, y, x)
    tw_tile = mi.ScalarTransform4f().scale(0.5)
    tw_full = mi.ScalarTransform4f().scale(1.5)

    inner = medium_grid(fz, tw_tile, res=(1, 1, 4), aabb=([0, 0, 0], [0.6, 0.7, 0.5]))
    repeat = {
        "type": "repeat",
        "inner": inner,
        "lattice_0": [0.5, 0.2, 0.0],
        "lattice_1": [0.1, 0.5, 0.0],
        "lattice_2": [0.0, 0.0, 0.5],
        "aabb_min": [0, 0, 0],
        "aabb_max": [1.5, 1.5, 1.5],
    }
    reference = medium_grid(np.tile(fz, (3, 1, 1)), tw_full, res=(1, 1, 12))
    assert_images_match(render(repeat, 7), render(reference, 17))
