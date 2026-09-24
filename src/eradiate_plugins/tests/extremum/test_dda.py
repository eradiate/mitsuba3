"""DDA route (``dda_init``/``dda_next``) on structures and media."""

import drjit as dr
import mitsuba as mi
import numpy as np
import pytest

MAX_STEPS = 512


def _f(x):
    return float(np.array(x).reshape(-1)[0])


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


def _structure_dict(kind, extremum_res):
    d = {"type": kind}
    if extremum_res is not None:
        d["resolution"] = mi.ScalarVector3i(*extremum_res)
    return d


def _extremum(volume, kind, extremum_res):
    extremum = mi.load_dict(_structure_dict(kind, extremum_res))
    extremum.update_extremum(volume.bbox(), volume)
    return extremum


def _component(volume, kind, extremum_res):
    """A medium whose extremum structure matches ``_extremum(volume, ...)``."""
    return mi.load_dict(
        {
            "type": "eoheterogeneous",
            "sigma_t": volume,
            "albedo": 0.8,
            "has_spectral_extinction": False,
            "extremum_structure": _structure_dict(kind, extremum_res),
        }
    )


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


def _clip(ray, bbox=(0.0, 1.0), maxt=4.0):
    """Entry/exit distances through the domain ``bbox``, or None if missed.

    ``sample_test`` does not clip to the structure's bbox, so the routes only
    agree over a range already clipped to the domain, as an integrator passes.
    """
    o = np.array(ray.o).reshape(-1)
    d = np.array(ray.d).reshape(-1)
    lo, hi = bbox
    with np.errstate(divide="ignore", invalid="ignore"):
        t0 = (lo - o) / d
        t1 = (hi - o) / d
    near = max(np.nanmax(np.minimum(t0, t1)), 0.0)
    far = min(np.nanmin(np.maximum(t0, t1)), maxt)
    return None if near >= far else (near, far)


def _walk_dda(extremum, ray, mint, maxt):
    state = extremum.dda_init(ray, mint, maxt)
    out = []
    for _ in range(MAX_STEPS):
        if not _f(state.mint) < _f(state.maxt):
            return out
        segment, state = extremum.dda_next(state)
        out.append(
            (
                _f(segment.mint),
                _f(segment.maxt),
                _f(segment.minorant()),
                _f(segment.majorant()),
            )
        )
    pytest.fail("dda traversal did not terminate")


def _sample(segments, mint, maxt, target_ot):
    """Delta-tracking oracle over a piecewise-constant majorant: returns
    ``(distance, leftover_ot)`` like ``sample_test``."""
    for lo, hi, _, majorant in segments:
        lo, hi = max(lo, mint), min(hi, maxt)
        if hi <= lo:
            continue
        ot = (hi - lo) * majorant
        if target_ot < ot:
            return lo + target_ot / max(majorant, 1e-7), target_ot
        target_ot -= ot
    return np.inf, target_ot


def _sum(walks, mint, maxt):
    """Pointwise sum of several piecewise-constant walks, as segments."""
    bounds = sorted(
        {mint, maxt} | {t for w in walks for s in w for t in s[:2] if mint < t < maxt}
    )
    out = []
    for lo, hi in zip(bounds[:-1], bounds[1:]):
        t = 0.5 * (lo + hi)
        majorant = sum(s[3] for w in walks for s in w if s[0] <= t < s[1])
        out.append((lo, hi, 0.0, majorant))
    return out


def _assert_close(expected, actual, ray, atol=0.0):
    for e, a in zip(expected, actual):
        assert np.allclose(_f(e), _f(a), rtol=1e-5, atol=atol), ray


def _structures():
    """``(name, structure, rays, bbox)`` for the three structure plugins."""
    return [
        (
            "global",
            _extremum(_grid_volume((4, 4, 4)), "extremum_global", None),
            _rays(),
            (0.0, 1.0),
        ),
        (
            "grid",
            _extremum(_grid_volume((4, 5, 3)), "extremum_grid", (2, 3, 3)),
            _rays(),
            (0.0, 1.0),
        ),
        (
            "radial",
            _extremum(_spherical_volume((4, 4, 4)), "extremum_spherical", (4, 1, 1)),
            _spherical_rays(),
            (-1.0, 1.0),
        ),
        (
            "spherical_3d",
            _extremum(_spherical_volume((4, 4, 4)), "extremum_spherical", (3, 4, 4)),
            _spherical_rays(),
            (-1.0, 1.0),
        ),
    ]


def test_structure_segments_tile(variant_scalar_rgb):
    for name, extremum, rays, _ in _structures():
        for ray in rays:
            segments = _walk_dda(extremum, ray, 0.0, 4.0)
            if not segments:
                continue
            assert segments[0][0] >= 0.0, name
            assert segments[-1][1] <= 4.0, name
            for before, after in zip(segments[:-1], segments[1:]):
                assert before[1] == after[0], (name, ray)


@pytest.mark.parametrize("target_ot", [0.0, 0.05, 0.5, 2.0, 1e6])
def test_structure_matches_traverse_extremum(variant_scalar_rgb, target_ot):
    for name, extremum, rays, bbox in _structures():
        for ray in rays:
            clipped = _clip(ray, bbox)
            if clipped is None:
                continue
            mint, maxt = clipped
            expected = extremum.sample_test(ray, mint, maxt, target_ot)
            actual = _sample(
                _walk_dda(extremum, ray, mint, maxt), mint, maxt, target_ot
            )
            # atol: `dda_init` re-derives the entry distance from its own bbox
            # test, which can differ from `mint` by an ulp.
            _assert_close(expected, actual, (name, ray), atol=1e-6)


@pytest.mark.parametrize("target_ot", [0.0, 0.5, 2.0])
def test_plain_medium_matches_structure(variant_scalar_rgb, target_ot):
    volume = _grid_volume((4, 5, 3))
    extremum = _extremum(volume, "extremum_grid", (2, 3, 3))
    medium = _component(volume, "extremum_grid", (2, 3, 3))
    for ray in _rays():
        clipped = _clip(ray)
        if clipped is None:
            continue
        mint, maxt = clipped
        expected = extremum.sample_test(ray, mint, maxt, target_ot)
        actual = medium.sample_test_dda(ray, mint, maxt, target_ot)
        _assert_close(expected, actual, ray)


def _overlap_parts():
    """Grid in [0, 1]^3, spherical and global in [-1, 1]^3: gaps, partial
    overlap and never-entered components all occur along one ray fan."""
    parts = [
        (_grid_volume((4, 5, 3)), "extremum_grid", (2, 3, 3)),
        (_spherical_volume((4, 4, 4)), "extremum_spherical", (3, 4, 4)),
        (_grid_volume((4, 4, 4), seed=3), "extremum_global", None),
    ]
    structures = [_extremum(*p) for p in parts]
    components = [_component(*p) for p in parts]
    return structures, components


def _overlap(components):
    return mi.load_dict(
        {"type": "overlap", **{f"m{i}": c for i, c in enumerate(components)}}
    )


@pytest.mark.parametrize("n", [1, 2, 3])
@pytest.mark.parametrize("target_ot", [0.0, 0.05, 0.5, 2.0, 1e6])
def test_overlap_matches_summed_components(variant_scalar_rgb, n, target_ot):
    structures, components = _overlap_parts()
    structures, medium = structures[:n], _overlap(components[:n])
    bbox = (0.0, 1.0) if n == 1 else (-1.0, 1.0)
    for ray in _spherical_rays():
        clipped = _clip(ray, bbox)
        if clipped is None:
            continue
        mint, maxt = clipped
        walks = [_walk_dda(s, ray, mint, maxt) for s in structures]
        expected = _sample(_sum(walks, mint, maxt), mint, maxt, target_ot)
        actual = medium.sample_test_dda(ray, mint, maxt, target_ot)
        _assert_close(expected, actual, (n, ray), atol=1e-6)


def test_overlap_vcall_matches_direct(variants_vec_rgb):
    """The state list survives a ``MediumPtr`` vcall per segment."""
    _, components = _overlap_parts()
    medium = _overlap(components)
    rays = [r for r in _spherical_rays() if _clip(r, (-1.0, 1.0)) is not None]
    bounds = np.array([_clip(r, (-1.0, 1.0)) for r in rays])
    o = mi.Point3f(np.array([np.array(r.o).reshape(-1) for r in rays]).T)
    d = mi.Vector3f(np.array([np.array(r.d).reshape(-1) for r in rays]).T)
    ray = mi.Ray3f(o=o, d=d)
    mint, maxt = mi.Float(bounds[:, 0]), mi.Float(bounds[:, 1])
    ptr = dr.full(mi.MediumPtr, medium, dr.width(mint))

    for target_ot in (0.0, 0.5, 2.0):
        expected = medium.sample_test_dda(ray, mint, maxt, target_ot)
        actual = ptr.sample_test_dda(ray, mint, maxt, target_ot)
        for e, a in zip(expected, actual):
            assert np.allclose(np.array(e), np.array(a), rtol=1e-5)
