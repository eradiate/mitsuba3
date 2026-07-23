import time

import pytest
from pytest import fixture
import drjit as dr
import mitsuba as mi
import numpy as np


rng = np.random.default_rng(0)


def make_phase_grid(n_r=10, n_v=10, min_pts=50, max_pts=200):
    n_entries = n_r * n_v
    grid_len   = rng.integers(min_pts, max_pts, n_entries).astype(np.uint32)
    grid_start = np.concatenate([[0], np.cumsum(grid_len[:-1])]).astype(np.uint32)
    total_pts  = int(grid_len.sum())

    nodes         = np.empty(total_pts)
    phase_mueller = np.empty((total_pts, 6))

    for i in range(n_entries):
        s, l = int(grid_start[i]), int(grid_len[i])
        cos_vals = np.sort(rng.uniform(-1.0, 1.0, l))
        cos_vals[0], cos_vals[-1] = -1.0, 1.0
        nodes[s : s + l] = cos_vals

        g  = rng.uniform(-0.9, 0.9)
        g2 = g * g
        denom = (1.0 + g2 - 2.0 * g * cos_vals) ** 1.5
        m11 = (1.0 - g2) / denom / (4.0 * np.pi)
        phase_mueller[s : s + l, 0] = m11
        phase_mueller[s : s + l, 1] = m11 * 0.1
        phase_mueller[s : s + l, 2] = m11 * 0.9
        phase_mueller[s : s + l, 3] = m11 * 0.8
        phase_mueller[s : s + l, 4] = m11 * 0.05
        phase_mueller[s : s + l, 5] = m11 * 0.8

    return (
        dr.scalar.ArrayXf64(nodes),
        dr.scalar.ArrayXu(grid_start),
        dr.scalar.ArrayXu(grid_len),
        dr.scalar.ArrayXf64(phase_mueller.flatten()),
    )


def make_rv_volumes(r_eff_bounds, v_eff_bounds, shape=(1, 1, 1)):

    if r_eff_bounds[0] == r_eff_bounds[1]:
        r_eff_val = np.full(shape, r_eff_bounds[0])
    else:
        r_eff_val = r_eff_bounds[0] + rng.random(shape) * (r_eff_bounds[1] - r_eff_bounds[0])
    if v_eff_bounds[0] == v_eff_bounds[1]:
        v_eff_val = np.full(shape, v_eff_bounds[0])
    else:
        v_eff_val = v_eff_bounds[0] + rng.random(shape) * (v_eff_bounds[1] - v_eff_bounds[0])

    return (
        { "type": "gridvolume", "grid": mi.VolumeGrid(r_eff_val), "filter_type": "nearest" },
        { "type": "gridvolume", "grid": mi.VolumeGrid(v_eff_val), "filter_type": "nearest" },
    )


def make_rv_grid(n_r=22, n_v=3, r_bounds=(4.0, 25.0), v_bounds=(0.1, 0.3), sigma_s_weight=None):

    r_eff_grid = np.linspace(r_bounds[0], r_bounds[1], n_r)
    v_eff_grid = np.linspace(v_bounds[0], v_bounds[1], n_v)

    if sigma_s_weight is None:
        sigma_s_weight = np.ones(n_r * n_v)
    sigma_s_weight = np.asarray(sigma_s_weight)

    return (
        mi.VolumeGrid(r_eff_grid.reshape(1, 1, -1, 1)),
        mi.VolumeGrid(v_eff_grid.reshape(1, 1, -1, 1)),
        mi.VolumeGrid(sigma_s_weight.reshape(1, 1, -1, 1)),
    )


def make_interaction(wi=(0, 0, 1), p=(0.5, 0.5, 0.5)):
    mei         = mi.MediumInteraction3f()
    mei.wi      = mi.Vector3f(wi)
    mei.p       = mi.Point3f(p)
    mei.sh_frame = mi.Frame3f(mei.wi)
    return mei


def make_query_interactions(n_points, r_bounds, v_bounds, wi=(0, 0, 1)):
    rv_volumes = make_rv_volumes(r_bounds, v_bounds, shape=(1, 1, n_points))
    meis = [
        make_interaction(wi=wi, p=(float(rng.random()), 0.5, 0.5))
        for _ in range(n_points)
    ]
    return rv_volumes, meis


def make_phase(
        r_eff_volume, v_eff_volume,
        r_eff_grid,   v_eff_grid,
        nodes, phase_mueller,
        grid_start, grid_len,
        blending_method,
        sigma_s_weight,
    ):
    phase = mi.load_dict({
        "type": "particlephase",
        "r_eff_volume": r_eff_volume,
        "v_eff_volume": v_eff_volume,
        "r_eff_grid": r_eff_grid,
        "v_eff_grid": v_eff_grid,
        "n_r": int(r_eff_grid.size()[0]),
        "n_v": int(v_eff_grid.size()[0]),
        "nodes": nodes,
        "phase_mueller": phase_mueller,
        "grid_start": grid_start,
        "grid_len": grid_len,
        "blending_method": blending_method,
        "sigma_s_weight": sigma_s_weight,
    })
    return phase


def build_grid(n_r, n_v, r_bounds=(10.0, 11.0), v_bounds=(0.11, 0.22), sigma_s_weight=None,
               min_pts=10, max_pts=11):
    nodes, grid_start, grid_len, phase_mueller = make_phase_grid(n_r=n_r, n_v=n_v, min_pts=min_pts, max_pts=max_pts)
    r_eff_grid, v_eff_grid, sigma_s_weight = make_rv_grid(
        n_r=n_r, n_v=n_v, r_bounds=r_bounds, v_bounds=v_bounds, sigma_s_weight=sigma_s_weight)
    return nodes, grid_start, grid_len, phase_mueller, r_eff_grid, v_eff_grid, sigma_s_weight


def build_phase_function(grid, rv_volumes, blending_method):
    nodes, grid_start, grid_len, phase_mueller, r_eff_grid, v_eff_grid, sigma_s_weight = grid
    r_eff, v_eff = rv_volumes

    phase = make_phase(
        r_eff, v_eff, r_eff_grid, v_eff_grid,
        nodes, phase_mueller,
        grid_start, grid_len,
        blending_method,
        sigma_s_weight,
    )

    return phase, r_eff, v_eff, r_eff_grid, v_eff_grid, nodes, phase_mueller, grid_start, grid_len, sigma_s_weight


def bilinear_corners(r_eff_volume, v_eff_volume, r_eff_grid, v_eff_grid, sigma_s_weight, mei):
    r_eff_grid_np     = np.array(r_eff_grid).ravel()
    v_eff_grid_np     = np.array(v_eff_grid).ravel()
    sigma_s_weight_np = np.array(sigma_s_weight).ravel()

    n_r = len(r_eff_grid_np)
    n_v = len(v_eff_grid_np)

    r_eff = float(mi.load_dict(r_eff_volume).eval_1(mei))
    v_eff = float(mi.load_dict(v_eff_volume).eval_1(mei))

    if n_r > 1:
        ir = int(np.clip(np.searchsorted(r_eff_grid_np, r_eff, side="right") - 1, 0, n_r - 2))
        tr = (r_eff - r_eff_grid_np[ir]) / (r_eff_grid_np[ir + 1] - r_eff_grid_np[ir])
    else:
        ir, tr = 0, 0.0

    if n_v > 1:
        iv = int(np.clip(np.searchsorted(v_eff_grid_np, v_eff, side="right") - 1, 0, n_v - 2))
        tv = (v_eff - v_eff_grid_np[iv]) / (v_eff_grid_np[iv + 1] - v_eff_grid_np[iv])
    else:
        iv, tv = 0, 0.0

    ir_hi = ir + 1 if n_r > 1 else ir
    iv_hi = iv + 1 if n_v > 1 else iv

    idx = np.array([
        ir    * n_v + iv,
        ir    * n_v + iv_hi,
        ir_hi * n_v + iv,
        ir_hi * n_v + iv_hi,
    ])
    w = np.array([
        (1.0 - tr) * (1.0 - tv),
        (1.0 - tr) * tv,
        tr * (1.0 - tv),
        tr * tv,
    ]) * sigma_s_weight_np[idx]
    w_sum = w.sum()
    if w_sum > 0:
        w = w / w_sum

    return idx, w


def entry_cdf(grid_start_np, grid_len_np, nodes_np, phase_mueller_np, e):
    s, l = int(grid_start_np[e]), int(grid_len_np[e])
    entry_nodes = nodes_np[s : s + l]
    m11 = phase_mueller_np[s : s + l, 0]
    increments = 0.5 * (m11[:-1] + m11[1:]) * np.diff(entry_nodes)
    cdf = np.concatenate(([0.0], np.cumsum(increments)))
    total = cdf[-1]
    if total > 0:
        cdf = cdf / total
    return cdf, total


def union_grid(idx, w, grid_start_np, grid_len_np, nodes_np, lo_bound, hi_bound):
    merged = [lo_bound, hi_bound]
    for e, wgt in zip(idx, w):
        if wgt == 0.0:
            continue
        e = int(e)
        s, l = int(grid_start_np[e]), int(grid_len_np[e])
        entry_nodes = nodes_np[s : s + l]
        mask = (entry_nodes >= lo_bound) & (entry_nodes <= hi_bound)
        merged.extend(entry_nodes[mask].tolist())
    return np.unique(np.array(merged))


def mock_phase(
        r_eff_volume, v_eff_volume,
        r_eff_grid,   v_eff_grid,
        nodes, phase_mueller,
        grid_start, grid_len,
        sigma_s_weight,
        mei,
    ):
    nodes_np         = np.array(nodes)
    grid_start_np    = np.array(grid_start)
    grid_len_np      = np.array(grid_len)
    phase_mueller_np = np.array(phase_mueller).reshape(-1, 6)

    idx, w = bilinear_corners(
        r_eff_volume, v_eff_volume, r_eff_grid, v_eff_grid, sigma_s_weight, mei)

    lo_bound = max(nodes_np[int(grid_start_np[int(e)])] for e in idx)
    hi_bound = min(nodes_np[int(grid_start_np[int(e)]) + int(grid_len_np[int(e)]) - 1] for e in idx)

    grid = union_grid(idx, w, grid_start_np, grid_len_np, nodes_np, lo_bound, hi_bound)

    columns = np.zeros((6, len(grid)))
    for e, wgt in zip(idx, w):
        if wgt == 0.0:
            continue
        e = int(e)
        s, l = int(grid_start_np[e]), int(grid_len_np[e])
        entry_nodes = nodes_np[s : s + l]
        _, total = entry_cdf(grid_start_np, grid_len_np, nodes_np, phase_mueller_np, e)
        norm = 1.0 / total if total > 0 else 0.0
        for c in range(6):
            columns[c] += wgt * norm * np.interp(grid, entry_nodes, phase_mueller_np[s : s + l, c])

    fmt = lambda v: ",".join(repr(float(x)) for x in v)

    return mi.load_dict({
        "type": "tabphase_polarized",
        "nodes": fmt(grid),
        "m11": fmt(columns[0]),
        "m12": fmt(columns[1]),
        "m22": fmt(columns[2]),
        "m33": fmt(columns[3]),
        "m34": fmt(columns[4]),
        "m44": fmt(columns[5]),
    })


def eval_and_sample(phase, mei):
    ctx = mi.PhaseFunctionContext(None)

    sample1 = float(rng.random())
    sample2 = mi.Point2f(float(rng.random()), float(rng.random()))

    wo, weight, pdf_sample = phase.sample(ctx, mei, sample1, sample2)
    val, pdf_eval = phase.eval_pdf(ctx, mei, wo)

    return wo, weight, pdf_sample, val, pdf_eval


def check_phase_consistency(
        phase,
        r_eff, v_eff, r_eff_grid, v_eff_grid,
        nodes, phase_mueller, grid_start, grid_len,
        sigma_s_weight,
        mei,
    ):
    wo, weight, pdf_sample, val, pdf_eval = eval_and_sample(phase, mei)

    assert np.isclose(float(pdf_sample), float(pdf_eval))

    mock = mock_phase(
        r_eff, v_eff, r_eff_grid, v_eff_grid,
        nodes, phase_mueller, grid_start, grid_len, sigma_s_weight, mei)
    ctx = mi.PhaseFunctionContext(None)
    val_mock, pdf_mock = mock.eval_pdf(ctx, mei, wo)

    assert np.isclose(float(pdf_eval), float(pdf_mock), rtol=1e-4)
    assert np.allclose(np.array(val), np.array(val_mock), rtol=1e-3, atol=1e-6)


def cell_query_points(r_nodes, v_nodes):
    points = []
    for i in range(len(r_nodes) - 1):
        r_lo, r_hi = r_nodes[i], r_nodes[i + 1]
        r_mid = 0.5 * (r_lo + r_hi)
        for j in range(len(v_nodes) - 1):
            v_lo, v_hi = v_nodes[j], v_nodes[j + 1]
            v_mid = 0.5 * (v_lo + v_hi)
            points.append((r_mid, v_mid))
            points.append((r_lo, v_mid))
            points.append((r_mid, v_lo))
    return points


def reference_eval_max(query_nodes, grid_start_np, grid_len_np, phase_mueller_np, nodes_np, n_entries):
    values = np.zeros(len(query_nodes))
    for e in range(n_entries):
        s, l = int(grid_start_np[e]), int(grid_len_np[e])
        entry_nodes = nodes_np[s : s + l]
        _, total = entry_cdf(grid_start_np, grid_len_np, nodes_np, phase_mueller_np, e)
        norm = 1.0 / total if total > 0 else 0.0
        m11 = np.interp(query_nodes, entry_nodes, phase_mueller_np[s : s + l, 0])
        values = np.maximum(values, m11 * norm / (2.0 * np.pi))
    return values


@fixture
def phase_grid_single(variant_scalar_mono_polarized_double):
    return make_phase_grid(n_r=1, n_v=1, min_pts=10, max_pts=11)

@fixture
def rv_volumes_single():
    return make_rv_volumes((10.0, 10.0), (0.11, 0.11))

@fixture
def interaction_centrenadir():
    return make_interaction(wi=(0, 0, 1), p=(0.5, 0.5, 0.5))

@fixture
def rv_grid_single():
    return make_rv_grid(
        n_r=1, n_v=1,
        r_bounds=(10.0, 10.0), v_bounds=(0.11, 0.11),
        sigma_s_weight=0.11,
    )

@fixture(params=[
        (0, 0, 1),
        (0, 0, -1),
        (1, 0, 0),
        (-1, 0, 0),
        (0, 1, 0),
        (0, -1, 0),
        (1, 1, 1),
        (0.27, -0.64, 0.53),
    ], ids=["nadir", "zenith", "+x", "-x", "+y", "-y", "oblique", "generic"])
def interaction_directions(request):
    wi = np.array(request.param, dtype=float)
    wi = wi / np.linalg.norm(wi)
    return make_interaction(wi=wi, p=(0.5, 0.5, 0.5))

@fixture(params=["bisect", "search", "stochastic", "tabulate", "reconstruct"])
def phase_function_single(request, phase_grid_single, rv_volumes_single, rv_grid_single):
    grid = (*phase_grid_single, *rv_grid_single)
    return build_phase_function(grid, rv_volumes_single, request.param)

@fixture(params=[(2, 1), (1, 2), (2, 2)], ids=["2x1", "1x2", "2x2"])
def grid_dual(request):
    n_r, n_v = request.param
    return build_grid(n_r, n_v)

@fixture(params=[(10.0, 0.11), (10.5, 0.165), (11.0, 0.22)], ids=["low", "mid", "high"])
def rv_volumes_dual(request):
    r_eff_val, v_eff_val = request.param
    return make_rv_volumes((r_eff_val, r_eff_val), (v_eff_val, v_eff_val))

@fixture(params=["bisect", "search", "stochastic", "tabulate", "reconstruct"])
def phase_function_dual(request, grid_dual, rv_volumes_dual):
    return build_phase_function(grid_dual, rv_volumes_dual, request.param)

@fixture(params=[
        [1.0,    1.0, 1.0, 1.0],
        [1000.0, 1.0, 1.0, 1.0],
        [1.0,    1.0, 1.0, 1000.0],
        [1e-6,   1.0, 1.0, 1.0],
    ], ids=["uniform", "corner0_dominant", "corner3_dominant", "corner0_suppressed"])
def grid_sigma_s(request):
    return build_grid(2, 2, sigma_s_weight=request.param)

@fixture(params=["bisect", "search", "stochastic", "tabulate", "reconstruct"])
def phase_function_sigma_s(request, grid_sigma_s):
    rv_volumes = make_rv_volumes((10.5, 10.5), (0.165, 0.165))
    return build_phase_function(grid_sigma_s, rv_volumes, request.param)

@fixture(params=[
        (100.0, 0.165, "out of dataset range"),
        (10.5, 100.0,  "out of dataset range"),
        (float("nan"), 0.165, "NaN"),
    ], ids=["r_eff_out_of_range", "v_eff_out_of_range", "nan_r_eff"])
def phase_invalid_query(request):
    r_eff_val, v_eff_val, match = request.param
    grid = build_grid(2, 2)
    rv_volumes = make_rv_volumes((r_eff_val, r_eff_val), (v_eff_val, v_eff_val))
    phase, *_ = build_phase_function(grid, rv_volumes, "bisect")
    return phase, match


def test_instantiate(phase_function_single):
    phase, *_ = phase_function_single
    assert phase is not None


def test_single_directions(phase_function_single, interaction_directions):
    check_phase_consistency(*phase_function_single, interaction_directions)


def test_dual_directions(phase_function_dual, interaction_directions):
    check_phase_consistency(*phase_function_dual, interaction_directions)


def test_sigma_s_weight(phase_function_sigma_s, interaction_centrenadir):
    check_phase_consistency(*phase_function_sigma_s, interaction_centrenadir)


def test_invalid_query_raises(phase_invalid_query, interaction_centrenadir):
    phase, match = phase_invalid_query
    ctx = mi.PhaseFunctionContext(None)
    with pytest.raises(Exception, match=match):
        phase.eval_pdf(ctx, interaction_centrenadir, mi.Vector3f([0, 0, -1]))


def test_sample_distribution(variant_scalar_mono_polarized_double):
    rng_local = np.random.default_rng(42)

    r_bounds, v_bounds = (4.0, 25.0), (0.11, 0.22)
    grid = build_grid(22, 2, r_bounds=r_bounds, v_bounds=v_bounds, min_pts=700, max_pts=1000)

    n_query = 100
    rv_volumes, meis = make_query_interactions(n_query, r_bounds, v_bounds)
    ctx = mi.PhaseFunctionContext(None)

    n_samples_per_query = 2000
    methods = ["bisect", "search", "stochastic", "tabulate"]
    cdf_methods = ["bisect", "search", "tabulate"]

    s1_all = rng_local.random((n_query, n_samples_per_query))
    s2_all = rng_local.random((n_query, n_samples_per_query, 2))

    cos_theta = np.empty((len(methods), n_query, n_samples_per_query))

    for m, method in enumerate(methods):
        phase, *_ = build_phase_function(grid, rv_volumes, method)

        t0 = time.perf_counter()
        for q, mei in enumerate(meis):
            for k in range(n_samples_per_query):
                s1 = float(s1_all[q, k])
                s2 = mi.Point2f(float(s2_all[q, k, 0]), float(s2_all[q, k, 1]))
                wo, weight, pdf = phase.sample(ctx, mei, s1, s2)
                cos_theta[m, q, k] = float(-dr.dot(wo, mei.wi))
        elapsed = time.perf_counter() - t0

        n_total = n_query * n_samples_per_query
        print(f"{method}: {n_total} samples in {elapsed*1e3:.2f} ms "
              f"({elapsed/n_total*1e6:.2f} us/sample)")

    failures = []

    ref_idx = methods.index(cdf_methods[0])
    ref_cos = cos_theta[ref_idx]
    for method in cdf_methods[1:]:
        idx = methods.index(method)
        close = np.isclose(cos_theta[idx], ref_cos, atol=1e-4, rtol=1e-4)
        n_bad = int((~close).sum())
        print(f"{method}: {int(close.sum())}/{close.size} samples matched {cdf_methods[0]} "
              f"to numerical precision ({n_bad} mismatches)")

        if n_bad:
            q_idx, k_idx = np.unravel_index(np.flatnonzero(~close)[:5], close.shape)
            failures.append(
                f"{method}: {n_bad}/{close.size} samples disagreed with {cdf_methods[0]} "
                f"beyond numerical precision; first mismatches at (query,sample)="
                f"{list(zip(q_idx.tolist(), k_idx.tolist()))[0]}"
            )

    stochastic_idx = methods.index("stochastic")
    cur_mean, cur_std = cos_theta[stochastic_idx].mean(axis=1), cos_theta[stochastic_idx].std(axis=1)
    ref_mean, ref_std = ref_cos.mean(axis=1), ref_cos.std(axis=1)

    good = np.isclose(cur_mean, ref_mean, atol=0.05) & np.isclose(cur_std, ref_std, rtol=0.15)
    bad = ~good
    print(f"stochastic: {int(good.sum())}/{n_query} interactions matched {cdf_methods[0]} in "
          f"distribution (up to 10 first bad query indices: {np.flatnonzero(bad).tolist()[:10]})")

    if bad.any():
        failures.append(
            f"stochastic: {int(bad.sum())}/{n_query} interactions mismatched vs {cdf_methods[0]} "
        )

    assert not failures, "\n".join(failures)


@pytest.mark.parametrize("n_r,n_v", [
    (n_r, n_v) for n_r in range(10, 60, 10) for n_v in range(10, 60, 10)
])
def test_grid_sizes(n_r, n_v, variant_scalar_mono_polarized_double):
    grid = build_grid(n_r, n_v)
    nodes, grid_start, grid_len, phase_mueller, r_eff_grid, v_eff_grid, sigma_s_weight = grid

    r_eff_grid_np = np.array(r_eff_grid).ravel()
    v_eff_grid_np = np.array(v_eff_grid).ravel()
    r_eff_val = float(r_eff_grid_np[len(r_eff_grid_np) // 2])
    v_eff_val = float(v_eff_grid_np[len(v_eff_grid_np) // 2])

    rv_volumes = make_rv_volumes((r_eff_val, r_eff_val), (v_eff_val, v_eff_val))
    mei = make_interaction(wi=(0, 0, 1))
    ctx = mi.PhaseFunctionContext(None)

    for method in ["bisect", "search", "stochastic", "reconstruct", "tabulate"]:
        phase, *_ = build_phase_function(grid, rv_volumes, method)
        val, pdf = phase.eval_pdf(ctx, mei, mi.Vector3f([0, 0, -1]))
        assert np.isfinite(float(pdf))
        assert float(pdf) > 0.0


def test_large_grid_sweep(variant_scalar_mono_polarized_double, interaction_directions):
    grid = build_grid(4, 4, r_bounds=(10.0, 14.0), v_bounds=(0.11, 0.30))
    nodes, grid_start, grid_len, phase_mueller, r_eff_grid, v_eff_grid, sigma_s_weight = grid

    r_eff_grid_np = np.array(r_eff_grid).ravel()
    v_eff_grid_np = np.array(v_eff_grid).ravel()
    query_points = cell_query_points(r_eff_grid_np, v_eff_grid_np)

    for blending_method in ["bisect", "search", "stochastic", "reconstruct", "tabulate"]:
        for r_val, v_val in query_points:
            rv_volumes = make_rv_volumes((r_val, r_val), (v_val, v_val))
            phase_bundle = build_phase_function(grid, rv_volumes, blending_method)
            check_phase_consistency(*phase_bundle, interaction_directions)


def test_get_nodes_and_eval_max(variant_scalar_mono_polarized_double):
    n_r, n_v = 3, 3
    grid = build_grid(n_r, n_v)
    nodes, grid_start, grid_len, phase_mueller, r_eff_grid, v_eff_grid, sigma_s_weight = grid
    rv_volumes = make_rv_volumes((10.5, 10.5), (0.165, 0.165))
    phase, *_ = build_phase_function(grid, rv_volumes, "bisect")

    nodes_np         = np.array(nodes)
    grid_start_np    = np.array(grid_start)
    grid_len_np      = np.array(grid_len)
    phase_mueller_np = np.array(phase_mueller).reshape(-1, 6)

    nds = phase.get_nodes()
    nds_np = np.array(nds)
    assert np.allclose(nds_np, np.unique(nodes_np))

    values = dr.scalar.ArrayXf64(np.zeros(len(nds_np)))
    phase.eval_max(nds, values)

    ref_values = reference_eval_max(nds_np, grid_start_np, grid_len_np, phase_mueller_np, nodes_np, n_r * n_v)

    assert np.allclose(np.array(values), ref_values, rtol=1e-4, atol=1e-8)


def test_degenerate_directions(variant_scalar_mono_polarized_double, interaction_directions):
    mei = interaction_directions
    ctx = mi.PhaseFunctionContext(None)
    grid = build_grid(2, 2)
    rv_volumes = make_rv_volumes((10.5, 10.5), (0.165, 0.165))

    for method in ["bisect", "search", "stochastic", "reconstruct", "tabulate"]:
        phase, r_eff, v_eff, r_eff_grid, v_eff_grid, nodes, phase_mueller, grid_start, grid_len, sigma_s_weight = \
            build_phase_function(grid, rv_volumes, method)
        mock = mock_phase(
            r_eff, v_eff, r_eff_grid, v_eff_grid,
            nodes, phase_mueller, grid_start, grid_len, sigma_s_weight, mei)

        for wo in (mi.Vector3f(mei.wi), -mi.Vector3f(mei.wi)):
            val, pdf = phase.eval_pdf(ctx, mei, wo)
            assert np.all(np.isfinite(np.array(val)))
            assert np.isfinite(float(pdf))
            assert float(pdf) >= 0.0

            val_mock, pdf_mock = mock.eval_pdf(ctx, mei, wo)
            assert np.isclose(float(pdf), float(pdf_mock), rtol=1e-4)
            assert np.allclose(np.array(val), np.array(val_mock), rtol=1e-3, atol=1e-6)


def test_llvm_variant(variant_llvm_ad_spectral_polarized):
    rng_local = np.random.default_rng(0)

    n_r, n_v, min_pts, max_pts = 2, 2, 10, 15
    n_entries = n_r * n_v
    grid_len   = rng_local.integers(min_pts, max_pts, n_entries).astype(np.uint32)
    grid_start = np.concatenate([[0], np.cumsum(grid_len[:-1])]).astype(np.uint32)
    total_pts  = int(grid_len.sum())

    nodes         = np.empty(total_pts)
    phase_mueller = np.empty((total_pts, 6))
    for i in range(n_entries):
        s, l = int(grid_start[i]), int(grid_len[i])
        cos_vals = np.sort(rng_local.uniform(-1.0, 1.0, l))
        cos_vals[0], cos_vals[-1] = -1.0, 1.0
        nodes[s : s + l] = cos_vals

        g = rng_local.uniform(-0.9, 0.9)
        g2 = g * g
        denom = (1.0 + g2 - 2.0 * g * cos_vals) ** 1.5
        m11 = (1.0 - g2) / denom / (4.0 * np.pi)
        phase_mueller[s : s + l, 0] = m11
        phase_mueller[s : s + l, 1] = m11 * 0.1
        phase_mueller[s : s + l, 2] = m11 * 0.9
        phase_mueller[s : s + l, 3] = m11 * 0.8
        phase_mueller[s : s + l, 4] = m11 * 0.05
        phase_mueller[s : s + l, 5] = m11 * 0.8

    r_bounds, v_bounds = (4.0, 25.0), (0.11, 0.22)
    r_eff_grid = mi.VolumeGrid(np.linspace(*r_bounds, n_r).reshape(1, 1, -1, 1))
    v_eff_grid = mi.VolumeGrid(np.linspace(*v_bounds, n_v).reshape(1, 1, -1, 1))
    sigma_s_weight = mi.VolumeGrid(np.ones(n_entries).reshape(1, 1, -1, 1))

    r_eff_val = np.full((1, 1, 1), 0.5 * (r_bounds[0] + r_bounds[1]))
    v_eff_val = np.full((1, 1, 1), 0.5 * (v_bounds[0] + v_bounds[1]))
    r_eff_volume = {"type": "gridvolume", "grid": mi.VolumeGrid(r_eff_val), "filter_type": "nearest"}
    v_eff_volume = {"type": "gridvolume", "grid": mi.VolumeGrid(v_eff_val), "filter_type": "nearest"}

    mei = mi.MediumInteraction3f()
    mei.wi = mi.Vector3f(0, 0, 1)
    mei.p = mi.Point3f(0.5, 0.5, 0.5)
    mei.sh_frame = mi.Frame3f(mei.wi)
    ctx = mi.PhaseFunctionContext(None)

    cdf_methods = ["bisect", "search", "tabulate"]
    cos_theta = {}

    for method in cdf_methods:
        phase = mi.load_dict({
            "type": "particlephase",
            "r_eff_volume": r_eff_volume,
            "v_eff_volume": v_eff_volume,
            "r_eff_grid": r_eff_grid,
            "v_eff_grid": v_eff_grid,
            "n_r": n_r,
            "n_v": n_v,
            "nodes": mi.Float(nodes),
            "phase_mueller": mi.Float(phase_mueller.flatten()),
            "grid_start": mi.UInt32(grid_start),
            "grid_len": mi.UInt32(grid_len),
            "blending_method": method,
            "sigma_s_weight": sigma_s_weight,
        })
        s1 = mi.Float(0.3)
        s2 = mi.Point2f(0.4, 0.6)
        wo, weight, pdf = phase.sample(ctx, mei, s1, s2)
        dr.eval(wo, weight, pdf)
        assert np.all(np.isfinite(np.array(dr.ravel(weight))))
        assert np.isfinite(float(pdf[0])) and float(pdf[0]) >= 0.0
        cos_theta[method] = float(-dr.dot(wo, mei.wi)[0])

    ref = cos_theta[cdf_methods[0]]
    for method in cdf_methods[1:]:
        assert np.isclose(cos_theta[method], ref, atol=1e-4, rtol=1e-4), (
            f"{method} disagreed with {cdf_methods[0]} under llvm_ad_spectral_polarized: "
            f"{cos_theta[method]} vs {ref}"
        )
