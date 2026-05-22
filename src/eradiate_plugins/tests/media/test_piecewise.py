import drjit as dr
import mitsuba as mi
import numpy as np


def gen_sigma_grid(medium_height, num_layers, integral, lbd):
    scale = integral / lbd
    z = np.linspace(0.0, medium_height, num_layers, False)
    ext_coefs = scale * np.exp(-z / lbd)
    ext_coefs = np.reshape(ext_coefs, (num_layers, 1, 1))
    return ext_coefs


# Mitsuba dictionary
def create_medium_dict():
    medium_height = 100000
    num_layers = 10
    integral = 3.0
    lbd_dist = 8300

    ext_coefs = gen_sigma_grid(medium_height, num_layers, integral, lbd_dist)
    exp_volume_grid = mi.VolumeGrid(ext_coefs)
    sigma_scale = 1.0

    medium_transform = (
        mi.ScalarTransform4f()
        .translate([-100000 / 2, -100000 / 2, 0.0])
        .scale([100000.0, 100000.0, 100000.0])
    )

    return mi.load_dict(
        {
            "type": "piecewise",
            "albedo": 0.8,
            "sigma_t": {
                "type": "gridvolume",
                "grid": exp_volume_grid,
                "use_grid_bbox": True,
                "filter_type": "nearest",
                "to_world": medium_transform,
            },
            "scale": sigma_scale,
        }
    )


def test01_sample_distances_down(variant_scalar_mono_double):
    medium = create_medium_dict()

    gt_dists = [0.0, 74578.93910366, 82117.51365975, 88515.28423884, 98460.51859237]
    gt_trs = [1.0, 0.75, 0.5, 0.25, 0.01]
    gt_pdfs = [
        7.06034444e-09,
        2.43563215e-05,
        5.41709938e-05,
        2.70854969e-05,
        3.61445783e-06,
    ]

    samples = [0.0, 0.25, 0.5, 0.75, 0.99]
    si = dr.zeros(mi.SurfaceInteraction3f)

    origins = mi.Point3f([0.0, 0.0, 100000.0])
    directions = mi.Vector3f([0.0, 0.0, -1.0])
    ray = mi.Ray3f(origins, directions)

    dists = []
    trs = []
    pdfs = []
    for sample in samples:
        mei, tr, pdf = medium.sample_interaction_analytical(ray, si, sample, 0, True)
        dists.append(mei.t)
        trs.append(tr)
        pdfs.append(pdf)

    assert dr.allclose(dists, gt_dists)
    assert dr.allclose(trs, gt_trs)
    assert dr.allclose(pdfs, gt_pdfs)


def test02_sample_distances_up(variant_scalar_mono_double):
    medium = create_medium_dict()

    gt_dists = [0.0, 7.95920400e02, 1.91770720e03, 3.83541440e03, 1.91443066e04]
    gt_trs = [1.0, 0.75, 0.5, 0.25, 0.01]
    gt_pdfs = [
        3.61445783e-04,
        2.71084337e-04,
        1.80722892e-04,
        9.03614458e-05,
        1.08341988e-06,
    ]

    samples = [0.0, 0.25, 0.5, 0.75, 0.99]
    si = dr.zeros(mi.SurfaceInteraction3f)

    origins = mi.Point3f([0.0, 0.0, 0.0])
    directions = mi.Vector3f([0.0, 0.0, 1.0])
    ray = mi.Ray3f(origins, directions)

    dists = []
    trs = []
    pdfs = []
    for sample in samples:
        mei, tr, pdf = medium.sample_interaction_analytical(ray, si, sample, 0, True)

        dists.append(mei.t)
        trs.append(tr)
        pdfs.append(pdf)

    assert dr.allclose(dists, gt_dists)
    assert dr.allclose(trs, gt_trs)
    assert dr.allclose(pdfs, gt_pdfs)


def test03_sample_distances_horizontal(variant_scalar_mono_double):
    medium = create_medium_dict()

    gt_dists = [0.0, 2655.31469158, 6397.77055374, 12795.54110748, 42505.86749422]
    gt_trs = [1.0, 0.75, 0.5, 0.25, 0.01]
    gt_pdfs = [
        1.08341988e-04,
        8.12564910e-05,
        5.41709940e-05,
        2.70854970e-05,
        1.08341988e-06,
    ]

    samples = [0.0, 0.25, 0.5, 0.75, 0.99]
    si = dr.zeros(mi.SurfaceInteraction3f)

    origins = mi.Point3f([0.0, 0.0, 15000.0])
    directions = mi.Vector3f([0.0, 1.0, 0.0])
    ray = mi.Ray3f(origins, directions)

    dists = []
    trs = []
    pdfs = []
    for sample in samples:
        mei, tr, pdf = medium.sample_interaction_analytical(ray, si, sample, 0, True)

        dists.append(mei.t)
        trs.append(tr)
        pdfs.append(pdf)

    assert dr.allclose(dists, gt_dists)
    assert dr.allclose(trs, gt_trs)
    assert dr.allclose(pdfs, gt_pdfs)


def test04_sample_distances_diag(variant_scalar_mono_double):
    medium = create_medium_dict()

    gt_dists = [9.23464998e00, 2.65531470e03, 6.39777058e03, 2.24562174e04, dr.inf]

    samples = [0.001, 0.25, 0.5, 0.75, 0.99]
    si = dr.zeros(mi.SurfaceInteraction3f)

    origins = mi.Point3f([0.0, 0.0, 15000.0])
    directions = mi.Vector3f([0.5, 0.5, 0.5])
    ray = mi.Ray3f(origins, directions)

    dists = []
    trs = []
    pdfs = []
    for idx, sample in enumerate(samples):
        mei, tr, pdf = medium.sample_interaction_analytical(ray, si, sample, 0, True)
        dists.append(mei.t)
        trs.append(tr)
        pdfs.append(pdf)

    for idx in range(len(gt_dists)):
        if dr.isinf(gt_dists[idx]) and dr.isinf(gt_dists[idx]) == dr.isinf(dists[idx]):
            gt_dists[idx] = 0.0
            dists[idx] = 0.0

    assert dr.allclose(dists, gt_dists)


def test05_sample_distances_heights(variant_scalar_mono_double):
    medium = create_medium_dict()

    gt_dists = [
        986.80067823,
        2824.88988516,
        3292.12110592,
        dr.inf,
        dr.inf,
        dr.inf,
        dr.inf,
    ]

    heights = [0.0, 9800.0, 10100.0, 27500.0, 40020.0, 75000, 98000.0]
    sample = 0.3
    direction = mi.Vector3f([0.0, 0.0, 1.0])
    rays = []

    si = dr.zeros(mi.SurfaceInteraction3f)

    for h in heights:
        origin = mi.Point3f([0.0, 0.0, h])
        rays.append(mi.Ray3f(origin, direction))

    dists = []
    trs = []
    pdfs = []
    for idx, ray in enumerate(rays):
        mei, tr, pdf = medium.sample_interaction_analytical(ray, si, sample, 0, True)
        dists.append(mei.t)
        trs.append(tr)
        pdfs.append(pdf)

    for idx in range(len(gt_dists)):
        if dr.isinf(gt_dists[idx]) and dr.isinf(gt_dists[idx]) == dr.isinf(dists[idx]):
            gt_dists[idx] = 0.0
            dists[idx] = 0.0

    assert dr.allclose(dists, gt_dists)


def test06_eval_transmittance_up(variant_scalar_mono_double):
    medium = create_medium_dict()

    gt_trs = [
        0.00573247,
        0.03493101,
        0.36588307,
        0.7398143,
        0.97331388,
        0.99930119,
        1.0,
    ]

    heights = [0.0, 5000.0, 15000.0, 25000.0, 45000.0, 75000, 100000.0]
    direction = [0.0, 0.0, 1.0]
    rays = []
    si = dr.zeros(mi.SurfaceInteraction3f, 1)
    si.t = dr.inf

    for h in heights:
        origin = [0.0, 0.0, h]
        rays.append(mi.Ray3f(origin, direction))

    trs = []
    for ray in rays:
        tr = medium.transmittance_eval_analytical(ray, si, True)
        trs.append(tr)

    assert dr.allclose(trs, gt_trs)


def test07_eval_transmittance_down(variant_scalar_mono):
    medium = create_medium_dict()

    gt_trs = [
        1.0,
        0.03865759,
        0.01566748,
        0.00663011,
        0.00588964,
        0.00573872,
        0.00573247,
    ]

    heights = [0.0, 9000.0, 15000.0, 29800.0, 45000.0, 70020, 100000.0]
    direction = [0.0, 0.0, -1.0]
    rays = []
    si = dr.zeros(mi.SurfaceInteraction3f, 1)
    si.t = dr.inf

    for h in heights:
        origin = [0.0, 0.0, h]
        rays.append(mi.Ray3f(origin, direction))

    trs = []
    for ray in rays:
        tr = medium.transmittance_eval_analytical(ray, si, True)
        trs.append(tr)

    assert dr.allclose(trs, gt_trs)


# ── RGB tests ─────────────────────────────────────────────────────────────────

# 2-layer unit medium, z ∈ [0, 1], sigma_t per channel.
_RGB_SIGMA_T = np.array(
    [[1.0, 2.0, 3.0], [3.0, 1.0, 2.0]],
    dtype=np.float64,
)
_LAYER_H = 0.5  # 2 layers over height 1.0


def _create_rgb_medium():
    n = len(_RGB_SIGMA_T)
    grid = _RGB_SIGMA_T.reshape(n, 1, 1, 3).astype(np.float32)
    return mi.load_dict(
        {
            "type": "piecewise",
            "albedo": 0.5,
            "sigma_t": {
                "type": "gridvolume",
                "grid": mi.VolumeGrid(grid),
                "use_grid_bbox": True,
                "filter_type": "nearest",
                "to_world": mi.ScalarTransform4f.translate([-0.5, -0.5, 0.0]),
            },
        }
    )


def _ref_tr_rgb(z_start, z_end):
    """Per-channel transmittance for segment [z_start, z_end] (upward ray)."""
    ot = np.zeros(3)
    for i, sigma in enumerate(_RGB_SIGMA_T):
        seg = min(z_end, (i + 1) * _LAYER_H) - max(z_start, i * _LAYER_H)
        if seg > 0:
            ot += sigma * seg
    return np.exp(-ot)


def _ref_distance_rgb(xi, channel):
    """
    Scatter distance along an upward ray from z=0 for sample xi in channel.
    Returns (t, layer_idx), or (inf, -1) if the ray escapes.
    """
    log_target = -np.log(1.0 - xi)
    cum = 0.0
    for i, sigma in enumerate(_RGB_SIGMA_T):
        contrib = sigma[channel] * _LAYER_H
        if cum + contrib >= log_target:
            return i * _LAYER_H + (log_target - cum) / sigma[channel], i
        cum += contrib
    return np.inf, -1


def test09_sample_interaction_rgb(variant_scalar_rgb):
    medium = _create_rgb_medium()
    si = dr.zeros(mi.SurfaceInteraction3f)
    ray = mi.Ray3f(mi.Point3f(0.0, 0.0, 0.0), mi.Vector3f(0.0, 0.0, 1.0))

    # (xi, channel): chosen so each channel scatters in a different layer
    for xi, channel in ((0.5, 0), (0.5, 1), (0.3, 2)):
        mei, tr, pdf = medium.sample_interaction_analytical(ray, si, xi, channel, True)
        expected_t, layer_idx = _ref_distance_rgb(xi, channel)
        expected_tr = _ref_tr_rgb(0.0, expected_t)
        expected_pdf = expected_tr * _RGB_SIGMA_T[layer_idx]

        assert dr.allclose(float(mei.t), float(expected_t), atol=1e-5)
        assert dr.allclose(tr, expected_tr.tolist(), atol=1e-5)
        assert dr.allclose(pdf, expected_pdf.tolist(), atol=1e-5)


def test08_transmittance_eval_rgb(variant_scalar_rgb):
    medium = _create_rgb_medium()
    si = dr.zeros(mi.SurfaceInteraction3f)

    # Three starting heights: full medium, top layer only, mid-layer start
    for z0 in (0.0, 0.5, 0.25):
        ray = mi.Ray3f(mi.Point3f(0.0, 0.0, z0), mi.Vector3f(0.0, 0.0, 1.0))
        tr = medium.transmittance_eval_analytical(ray, si, True)
        expected = _ref_tr_rgb(z0, 1.0)
        assert dr.allclose(tr, expected.tolist(), atol=1e-5)
