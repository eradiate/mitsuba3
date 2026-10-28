"""Extremum structure traversal (``begin``/``next``) and medium tracking."""

import drjit as dr
import mitsuba as mi
import numpy as np
import pytest

from .utils import FAR, walk_iter


def _grid_volume(resolution, seed=0):
    rng = np.random.default_rng(seed)
    data = (rng.random(resolution) + 0.1).astype(np.float32)
    return mi.load_dict(
        {
            "type": "gridvolume",
            "grid": mi.VolumeGrid(data.transpose(2, 1, 0)),
            "filter_type": "nearest",
            "accel": False,
        }
    )


def _spherical_volume(resolution, rmin=0.2, rmax=1.0, seed=0):
    rng = np.random.default_rng(seed)
    data = (rng.random(resolution) + 0.1).astype(np.float32)
    return mi.load_dict(
        {
            "type": "sphericalcoordsvolume",
            "volume": {
                "type": "gridvolume",
                "grid": mi.VolumeGrid(data.transpose(2, 1, 0)),
                "filter_type": "nearest",
                "accel": False,
            },
            "rmin": rmin,
            "rmax": rmax,
            "fillmin": 0.7,
            "fillmax": 0.2,
        }
    )


def _extremum(volume, kind, extremum_res):
    d = {"type": kind}
    if extremum_res is not None:
        d["resolution"] = mi.ScalarVector3i(*extremum_res)
    extremum = mi.load_dict(d)
    extremum.update_extremum(volume.bbox(), volume)
    return extremum


def _rays():
    """Oblique rays plus axis-aligned, inside-start and missing ones."""
    rng = np.random.default_rng(1)
    out = []
    for _ in range(8):
        o = rng.random(3) * 2.0 - 0.5
        d = rng.random(3) * 2.0 - 1.0
        out.append((o, d / np.linalg.norm(d)))
    out += [
        ([-0.5, 0.3, 0.7], [1.0, 0.0, 0.0]),
        ([1.5, 0.3, 0.7], [-1.0, 0.0, 0.0]),
        ([0.25, 0.25, 0.25], [0.0, 0.0, 1.0]),
        ([0.5, -1.0, 0.5], [0.0, 1.0, 0.0]),
        ([2.0, 2.0, 2.0], [1.0, 1.0, 1.0]),
    ]
    return [
        mi.Ray3f(o=mi.Point3f(*map(float, o)), d=mi.Vector3f(*map(float, d)))
        for o, d in out
    ]


def _spherical_rays():
    """Rays through a domain centred on the origin, plus the polar axis, the
    equatorial plane and a near-axis azimuth jump."""
    rng = np.random.default_rng(2)
    out = []
    for _ in range(8):
        o = rng.random(3) * 3.0 - 1.5
        d = rng.random(3) * 2.0 - 1.0
        out.append((o, d / np.linalg.norm(d)))
    out += [
        ([0.0, 0.0, 1.5], [0.0, 0.0, -1.0]),
        ([0.0, 0.0, -1.5], [0.0, 0.0, 1.0]),
        ([-1.5, 0.5, 0.0], [1.0, 0.0, 0.0]),
        ([-1.5, 0.05, 0.3], [1.0, 0.0, 0.0]),
        ([1.5, 1.5, 1.5], [1.0, 1.0, 1.0]),
    ]
    return [
        mi.Ray3f(o=mi.Point3f(*map(float, o)), d=mi.Vector3f(*map(float, d)))
        for o, d in out
    ]


def _structures():
    """``(name, structure, rays, bbox)`` for the three structure plugins."""
    return [
        (
            "global",
            _extremum(_grid_volume((4, 4, 4)), "globalextremum", None),
            _rays(),
            (0.0, 1.0),
        ),
        (
            "grid",
            _extremum(_grid_volume((4, 5, 3)), "gridextremum", (2, 3, 3)),
            _rays(),
            (0.0, 1.0),
        ),
        (
            "radial",
            _extremum(
                _spherical_volume((4, 4, 4)), "sphericalcoordsextremum", (4, 1, 1)
            ),
            _spherical_rays(),
            (-1.0, 1.0),
        ),
        (
            "spherical_3d",
            _extremum(
                _spherical_volume((4, 4, 4)), "sphericalcoordsextremum", (3, 4, 4)
            ),
            _spherical_rays(),
            (-1.0, 1.0),
        ),
    ]


def test_structure_segments_tile(variant_scalar_rgb):
    """
    Tests that all structures segments that tile perfectly.
    """
    for name, extremum, rays, _ in _structures():
        for ray in rays:
            segments = walk_iter(extremum, ray, 0.0, 4.0)
            if not segments:
                continue

            assert segments[0][0] >= 0.0, name
            assert segments[-1][1] <= 4.0, name
            for before, after in zip(segments[:-1], segments[1:]):
                assert before[1] == after[0], (name, ray)


def _overlap(components):
    return mi.load_dict(
        {"type": "overlap", **{f"m{i}": c for i, c in enumerate(components)}}
    )


@pytest.mark.parametrize(
    "o,d,maxt,length",
    [
        ([0.5, 0.5, 0.5], [1, 0, 0], dr.inf, 0.5),
        ([-1.0, 0.5, 0.5], [1, 0, 0], 1.25, 0.25),
        ([2.0, 0.5, 0.5], [1, 0, 0], dr.inf, 0.0),
        ([-2.0, 0.5, 0.5], [1, 0, 0], 1.0, 0.0),
        ([-1.0, 5.0, 0.5], [1, 0, 0], dr.inf, 0.0),
    ],
    ids=["inside", "clipped", "behind", "beyond_maxt", "miss"],
)
def test_track_range(variant_scalar_rgb, o, d, maxt, length):
    """Test that the tracking function properly covers the medium's box range.
    Uses transmittance as a proxy of distance, since it is exact with constant
    media."""
    medium = _constant(2.0, 0.5)
    ray = mi.Ray3f(mi.Ray3f(o=mi.Point3f(o), d=mi.Vector3f(d)), maxt)
    sample = _tracking_sample(0)
    tr = medium.sample_transmittance(ray, sample, 0)

    assert dr.allclose(tr, dr.exp(-2.0 * length))


@pytest.mark.parametrize(
    "o,d,maxt,ot",
    [
        ([-1.0, 0.5, 0.5], [1, 0, 0], dr.inf, 1.0 + 2.0 + 2.0),
        ([-1.0, 0.5, 0.5], [1, 0, 0], 1.75, 0.75 + 0.5),
        ([2.0, 0.5, 0.5], [1, 0, 0], dr.inf, 2.0),
        ([2.75, 0.5, -1.0], [0, 0, 1], dr.inf, 4.0),
        ([-1.0, 1.5, 0.5], [1, 0, 0], dr.inf, 0.0),
    ],
    ids=["all", "clipped", "from_gap", "one_component", "miss"],
)
def test_track_overlap_range(variant_scalar_rgb, o, d, maxt, ot):
    """Joint traversal of ``overlap`` components that partially overlap and
    leave a gap: sigma_t 1 on x in [0, 1], 2 on [0.5, 1.5] and 4 on [2.5, 3]
    (y, z in [0, 1]). As in ``test_track_range``, residual ratio tracking is
    exact in constant media: ``exp(-ot)``."""
    medium = _overlap(
        [
            _constant(1.0, 0.5),
            _constant(2.0, 0.5, ((0.5, 0.0, 0.0), (1.5, 1.0, 1.0))),
            _constant(4.0, 0.5, ((2.5, 0.0, 0.0), (3.0, 1.0, 1.0))),
        ]
    )

    ray = mi.Ray3f(mi.Ray3f(o=mi.Point3f(o), d=mi.Vector3f(d)), maxt)
    sample = _tracking_sample(0)
    tr = medium.sample_transmittance(ray, sample, 0)

    assert dr.allclose(tr, dr.exp(-ot))


def _tracking_sample(seed):
    """First free-flight sample, drawn from the generator of stream ``seed``."""
    rng = mi.PCG32()
    rng.seed(initseq=mi.UInt64(seed))
    return rng.next_float32()


def _track_test(medium, ray, seed, tracking):
    """``(distance, component, weight, transmittance, pdf_ratio)`` of the
    medium method that ``tracking`` names, ``None`` where it does not produce
    a field: ``delta`` is ``sample_real_interaction``, ``ratio`` is
    ``sample_transmittance_pdf`` (``weight`` is its throughput) and ``rrt`` is
    ``sample_transmittance``."""
    sample = _tracking_sample(seed)
    if tracking == "delta":
        ms, weight = medium.sample_real_interaction(ray, sample, 0)
        return ms.mei.t, ms.component, weight, None, ms.pdf_ratio
    if tracking == "ratio":
        tr, throughput, pdf_ratio = medium.sample_transmittance_pdf(ray, sample, 0)
        return None, None, throughput, tr, pdf_ratio
    tr = medium.sample_transmittance(ray, sample, 0)
    return None, None, None, tr, None


TRACKING = ["delta", "ratio", "rrt"]


def _constant(sigma_t, albedo, aabb=((0.0,) * 3, (1.0,) * 3)):
    return mi.load_dict(
        {
            "type": "eoheterogeneous",
            "sigma_t": {"type": "constvolume", "value": sigma_t},
            "albedo": albedo,
            "has_spectral_extinction": False,
            "aabb_min": mi.ScalarPoint3f(*aabb[0]),
            "aabb_max": mi.ScalarPoint3f(*aabb[1]),
        }
    )


def test_track_component_frequencies(variants_vec_rgb):
    """Real collisions in ``overlap`` of two constant media pick each
    component in proportion to its extinction, and weigh it by its albedo
    alone: chi-squared on the counts, means of the throughput."""

    n = 100_000
    sigma_t, albedo = (1.0, 3.0), (0.9, 0.3)
    medium = _overlap([_constant(s, a) for s, a in zip(sigma_t, albedo)])

    ray = mi.Ray3f(o=mi.Point3f(0.5, 0.5, -1.0), d=mi.Vector3f(0.0, 0.0, 1.0))
    ray = dr.gather(mi.Ray3f, ray, dr.zeros(mi.UInt32, n))

    distance, component, throughput, *_ = _track_test(
        medium, ray, dr.arange(mi.UInt32, n), "delta"
    )

    real = dr.isfinite(distance)
    counts = [dr.count(real & (component == i))[0] for i in range(2)]
    expected = [sum(counts) * s / sum(sigma_t) for s in sigma_t]
    chi2 = sum((c - e) ** 2 / e for c, e in zip(counts, expected))

    assert chi2 < 15.0, counts
    for i in range(2):
        picked = real & (component == i)
        weight = dr.select(picked, throughput.x, albedo[i])
        assert dr.allclose(weight, albedo[i], rtol=1e-5)


def _repeat(inner, lattice=None, aabb=FAR):
    d = {"type": "repeat", "inner": inner}
    if lattice is not None:
        d["lattice"] = mi.ScalarVector3f(*lattice)
    if aabb is not None:
        d["aabb_min"] = mi.ScalarPoint3f(*aabb[0])
        d["aabb_max"] = mi.ScalarPoint3f(*aabb[1])
    return mi.load_dict(d)


def _step_grid(**props):
    """sigma_t 1 on x in [0, 0.5] and 2 on [0.5, 1] (y, z in [0, 1]), under
    a global majorant of 2: null collisions on the first half. ``props`` are
    added to the medium's properties."""
    return mi.load_dict(
        {
            "type": "eoheterogeneous",
            "sigma_t": {
                "type": "gridvolume",
                "grid": mi.VolumeGrid(mi.TensorXf([1.0, 2.0], shape=(1, 1, 2, 1))),
                "filter_type": "nearest",
            },
            "albedo": 0.8,
            "has_spectral_extinction": False,
            "extremum": {"type": "globalextremum"},
            **props,
        }
    )


def _shells():
    """sigma_t 1 for r < 0.5, 2 for r in [0.5, 0.75] and 0.5 for r in
    [0.75, 1], under a single radial majorant of 2 on [0.5, 1]: null
    collisions in the outer shell."""
    return mi.load_dict(
        {
            "type": "eoheterogeneous",
            "sigma_t": {
                "type": "sphericalcoordsvolume",
                "volume": {
                    "type": "gridvolume",
                    "grid": mi.VolumeGrid(mi.TensorXf([2.0, 0.5], shape=(1, 1, 2, 1))),
                    "filter_type": "nearest",
                },
                "rmin": 0.5,
                "rmax": 1.0,
                "fillmin": 1.0,
            },
            "albedo": 0.8,
            "has_spectral_extinction": False,
            "extremum": {"type": "sphericalcoordsextremum"},
        }
    )


# name: (medium, ray origin, [(t, optical depth over [0, t])]) for a ray
# along +x; the last checkpoint is past the medium. Every component has null
# collisions where it overlaps another, so that a weight computed with the
# wrong component biases the estimates.
TRACKING_CASES = {
    "grid": (_step_grid, (-1.0, 0.5, 0.5), [(1.5, 0.5), (2.0, 1.5)]),
    "overlap": (
        lambda: _overlap(
            [_step_grid(), _constant(0.5, 0.8, ((0.25, 0.0, 0.0), (1.25, 1.0, 1.0)))]
        ),
        (-1.0, 0.5, 0.5),
        [(1.25, 0.25), (1.5, 0.625), (2.0, 1.875), (2.25, 2.0)],
    ),
    "repeat": (
        lambda: _repeat(
            _step_grid(), (1.5,) * 3, ((-0.5, -1.0, -1.0), (2.75, 2.0, 2.0))
        ),
        (-1.0, 0.5, 0.5),
        [(1.5, 0.5), (2.0, 1.5), (3.0, 2.0), (3.5, 3.0)],
    ),
    "shells": (
        _shells,
        (-2.0, 0.0, 0.0),
        [(1.25, 0.125), (1.5, 0.625), (2.5, 1.625), (2.75, 2.125), (3.0, 2.25)],
    ),
}


def _tracking_case(name):
    """``(medium, ray, checkpoints)`` of ``TRACKING_CASES[name]``."""
    medium, o, checkpoints = TRACKING_CASES[name]
    ray = mi.Ray3f(o=mi.Point3f(o), d=mi.Vector3f(1.0, 0.0, 0.0))
    return medium(), ray, checkpoints


def _track_seeds(medium, ray, n, tracking):
    """``(distance, weight, transmittance, pdf_ratio)`` of ``_track_test``
    over seeds ``0..n-1`` as arrays, spectra reduced to their first channel,
    ``None`` where ``_track_test`` gives ``None``: one wide call, or one call
    per seed in scalar variants."""
    if dr.is_jit_v(mi.Float):
        ray = dr.gather(mi.Ray3f, ray, dr.zeros(mi.UInt32, n))
        rec = _track_test(medium, ray, dr.arange(mi.UInt32, n), tracking)
        return [np.array(rec[0])] + [
            None if rec[i] is None else np.array(rec[i][0]) for i in (2, 3, 4)
        ]
    recs = [_track_test(medium, ray, seed, tracking) for seed in range(n)]
    return [np.array([r[0] for r in recs])] + [
        None if recs[0][i] is None else np.array([r[i][0] for r in recs])
        for i in (2, 3, 4)
    ]


@pytest.mark.parametrize("case", TRACKING_CASES)
@pytest.mark.parametrize("tracking", TRACKING)
def test_track_unbiased(variants_all_rgb, case, tracking):
    """Ratio tracking estimates the transmittance ``exp(-ot)``; delta tracking
    samples real collisions with CDF ``1 - exp(-ot(t))`` and weighs them by
    the albedo. Means agree within 5 standard errors over 100k lanes. Grey
    majorants leave ``pdf_ratio`` at one."""
    n, z = 100_000, 5.0
    medium, ray, checkpoints = _tracking_case(case)
    distance, weight, tr, pdf_ratio = _track_seeds(medium, ray, n, tracking)
    if tracking != "rrt":
        assert np.allclose(pdf_ratio, 1.0, rtol=1e-5)

    if tracking != "delta":
        expected = np.exp(-checkpoints[-1][1])
        assert abs(tr.mean() - expected) < z * tr.std() / np.sqrt(n) + 1e-6
    else:
        for t, ot in checkpoints:
            p = 1.0 - np.exp(-ot)
            assert abs(np.mean(distance < t) - p) < z * np.sqrt(p * (1 - p) / n), t
        real = np.isfinite(distance)
        assert np.allclose(weight[real], 0.8, rtol=1e-5)
        assert np.allclose(weight[~real], 1.0, rtol=1e-5)


@pytest.mark.parametrize("tracking", ["ratio", "rrt"])
def test_ratio_track_russian_roulette(variants_all_rgb, tracking):
    """With an ``rr_threshold`` above the transmittance exp(-1.5) of the step
    grid, Russian roulette kills more estimates than ratio tracking alone,
    and the mean still agrees with exp(-1.5) within 5 standard errors."""
    n, z = 100_000, 5.0
    ray = mi.Ray3f(o=mi.Point3f(-1.0, 0.5, 0.5), d=mi.Vector3f(1.0, 0.0, 0.0))
    tr = _track_seeds(_step_grid(rr_threshold=0.5), ray, n, tracking)[2]
    tr_ref = _track_seeds(_step_grid(), ray, n, tracking)[2]

    assert np.mean(tr == 0.0) > np.mean(tr_ref == 0.0)
    assert abs(tr.mean() - np.exp(-1.5)) < z * tr.std() / np.sqrt(n)


@pytest.mark.parametrize("case", TRACKING_CASES)
@pytest.mark.parametrize("tracking", TRACKING)
def test_track_vcall_matches_host(variants_vec_rgb, case, tracking):
    """Tracking through a ``MediumPtr`` vcall, lanes diverging, against the
    host pointer: same draws, same result. Tracing through the vcall also
    guards against loop re-recording errors (a radial structure stepping its
    unused indices raised one)."""
    n = 1024
    medium, ray, _ = _tracking_case(case)
    ray = dr.gather(mi.Ray3f, ray, dr.zeros(mi.UInt32, n))
    seed = dr.arange(mi.UInt32, n)
    expected = _track_test(medium, ray, seed, tracking)
    actual = _track_test(dr.full(mi.MediumPtr, medium, n), ray, seed, tracking)
    for e, a in zip(expected, actual):
        assert dr.all(e == a, axis=None)
