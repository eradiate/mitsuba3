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
    """Entry/exit distances through the domain ``bbox``, or None if missed:
    the range an integrator passes."""
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
    ``(distance, leftover_ot)`` like ``sample_test_dda``."""
    for lo, hi, _, majorant in segments:
        lo, hi = max(lo, mint), min(hi, maxt)
        if hi <= lo:
            continue
        ot = (hi - lo) * majorant
        if target_ot < ot:
            return lo + target_ot / max(majorant, 1e-7), target_ot
        target_ot -= ot
    return np.inf, target_ot


def _sample_dda(structure, ray, mint, maxt, target_ot):
    return _sample(_walk_dda(structure, ray, mint, maxt), mint, maxt, target_ot)


def _sum(walks, mint, maxt):
    """Pointwise sum of several piecewise-constant walks, as segments."""
    bounds = sorted(
        {mint, maxt} | {t for w in walks for s in w for t in s[:2] if mint < t < maxt}
    )
    out = []
    for lo, hi in zip(bounds[:-1], bounds[1:]):
        t = 0.5 * (lo + hi)
        inside = [s for w in walks for s in w if s[0] <= t < s[1]]
        out.append((lo, hi, sum(s[2] for s in inside), sum(s[3] for s in inside)))
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
        expected = _sample_dda(extremum, ray, mint, maxt, target_ot)
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


def _track_media():
    return [
        _component(_grid_volume((4, 5, 3)), "extremum_grid", (2, 3, 3)),
        _component(_spherical_volume((4, 4, 4)), "extremum_spherical", (4, 1, 1)),
        _component(_spherical_volume((4, 4, 4)), "extremum_spherical", (3, 4, 4)),
    ]


@pytest.mark.parametrize(
    "o,d,maxt,expected",
    [
        ([0.5, 0.5, 0.5], [1, 0, 0], np.inf, (0.0, 0.5)),
        ([-1.0, 0.5, 0.5], [1, 0, 0], 1.25, (1.0, 1.25)),
        ([2.0, 0.5, 0.5], [1, 0, 0], np.inf, (0.0, np.inf)),
        ([-2.0, 0.5, 0.5], [1, 0, 0], 1.0, (0.0, np.inf)),
        ([-1.0, 5.0, 0.5], [1, 0, 0], np.inf, (0.0, np.inf)),
    ],
    ids=["inside", "clipped", "behind", "beyond_maxt", "miss"],
)
def test_prepare_medium_traversal(variant_scalar_rgb, o, d, maxt, expected):
    """The range is clipped to ``[0, ray.maxt]``; an empty one is ``[0, inf]``,
    so the lane counts as escaped."""
    medium = _component(_grid_volume((2, 2, 2)), "extremum_global", None)
    ray = mi.Ray3f(mi.Ray3f(o=mi.Point3f(o), d=mi.Vector3f(d)), maxt)
    _, mint, maxt = medium.prepare_medium_traversal(ray)
    assert np.allclose((_f(mint), _f(maxt)), expected)


def _track(medium, ray, seed, ratio):
    """Python delta / residual ratio tracking over ``_walk_medium``, with the
    draws of ``track_test``: ``(distance, throughput)``. Non-spectral media,
    ``use_rrt`` on."""
    hit, mint, maxt = medium.intersect_aabb(ray)
    mint, maxt = max(_f(mint), 0.0), _f(maxt)
    rng = mi.PCG32(initseq=seed)
    if not (hit and mint < maxt < np.inf):
        return np.inf, np.ones(3)
    target = -np.log(1.0 - rng.next_float32())
    throughput = np.ones(3)
    for lo, hi, minorant, majorant in _walk_medium(medium, ray, mint, maxt):
        control = minorant if ratio else 0.0
        residual = majorant - control
        start = lo
        while target < (hi - start) * residual:
            t = start + target / max(residual, 2.0**-24)
            throughput *= np.exp(-(t - start) * control)
            mei = dr.zeros(mi.MediumInteraction3f)
            mei.t, mei.p = t, ray(t)
            sigma_s, _, sigma_t = (
                np.array(x).reshape(-1) for x in medium.get_scattering_coefficients(mei)
            )
            if ratio:
                throughput *= np.maximum(1.0 - (sigma_t - control) / residual, 0.0)
            else:
                rng.next_float32()
                if not rng.next_float32() < np.mean((majorant - sigma_t) / majorant):
                    return t, throughput * sigma_s / sigma_t
            target = -np.log(1.0 - rng.next_float32())
            start = t
        throughput *= np.exp(-(hi - start) * control)
        target -= (hi - start) * residual
    return np.inf, throughput


@pytest.mark.parametrize("ratio", [False, True])
def test_dda_track_matches_python_tracking(variant_scalar_rgb, ratio):
    """Same draws: ``dda_track`` samples the collisions of the Python tracker,
    including lanes that stay in a segment across null ones."""
    for medium in _track_media():
        for ray in _spherical_rays() + _rays():
            for seed in range(16):
                expected = _track(medium, ray, seed, ratio)
                actual = medium.track_test(ray, seed, ratio)
                for e, a in zip(expected, actual):
                    assert np.allclose(
                        np.array(e), np.array(a).reshape(-1), rtol=1e-4, atol=1e-5
                    ), (ray, seed)


@pytest.mark.parametrize("ratio", [False, True])
def test_dda_track_vcall_matches_scalar(variants_vec_rgb, ratio):
    """``dda_track`` through a ``MediumPtr`` vcall, lanes diverging, against
    the same rays one at a time in ``scalar_rgb``: same draws, same result."""
    rays = _spherical_rays() + _rays()
    n_seeds = 16
    variant = mi.variant()
    mi.set_variant("scalar_rgb")
    expected = [
        [
            [np.array(x).reshape(-1) for x in medium.track_test(ray, seed, ratio)]
            for seed in range(n_seeds)
            for ray in _spherical_rays() + _rays()
        ]
        for medium in _track_media()
    ]
    mi.set_variant(variant)

    o = np.array([np.array(r.o).reshape(-1) for r in rays] * n_seeds).T
    d = np.array([np.array(r.d).reshape(-1) for r in rays] * n_seeds).T
    ray = mi.Ray3f(o=mi.Point3f(o), d=mi.Vector3f(d))
    seed = mi.UInt32(np.repeat(np.arange(n_seeds), len(rays)))

    for medium, lanes in zip(_track_media(), expected):
        ptr = dr.full(mi.MediumPtr, medium, dr.width(seed))
        distance, throughput = ptr.track_test(ray, seed, ratio)
        distance, throughput = np.array(distance), np.array(throughput)
        for i, (e_distance, e_throughput) in enumerate(lanes):
            assert np.allclose(e_distance, distance[i], rtol=1e-5, atol=1e-6), i
            assert np.allclose(e_throughput, throughput[:, i], rtol=1e-5, atol=1e-6), i


def test_dda_track_in_symbolic_loop(variants_vec_rgb):
    """``dda_track`` nested in an outer symbolic loop, as in an integrator:
    re-recording the outer loop must not change what the inner loop presumed
    constant (a radial structure stepping its unused indices did)."""
    rays = _spherical_rays()
    o = np.array([np.array(r.o).reshape(-1) for r in rays]).T
    d = np.array([np.array(r.d).reshape(-1) for r in rays]).T
    ray = mi.Ray3f(o=mi.Point3f(o), d=mi.Vector3f(d))

    for medium in _track_media():
        ptr = dr.full(mi.MediumPtr, medium, len(rays))
        i, total = dr.while_loop(
            (dr.zeros(mi.UInt32, len(rays)), dr.zeros(mi.Float, len(rays))),
            lambda i, total: i < 2,
            lambda i, total, ptr=ptr: (
                i + 1,
                total + ptr.track_test(ray, i, False)[1].x,
            ),
        )
        dr.eval(total)


FAR = ((-20.0,) * 3, (20.0,) * 3)


def _repeat(inner, lattice=None, aabb=FAR):
    d = {"type": "repeat", "inner": inner}
    if lattice is not None:
        d["lattice"] = mi.ScalarVector3f(*lattice)
    if aabb is not None:
        d["aabb_min"] = mi.ScalarPoint3f(*aabb[0])
        d["aabb_max"] = mi.ScalarPoint3f(*aabb[1])
    return mi.load_dict(d)


def _bbox(structure):
    bbox = structure.bbox()
    return np.array(bbox.min).reshape(-1), np.array(bbox.extents()).reshape(-1)


def _repeat_cases():
    """``(name, medium, entries, domain)``. ``entries`` are the medium's
    structures as ``(structure, tiling)``, ``tiling`` being ``(cell_min,
    lattice, region)`` or None; ``domain`` is the box the rays aim at. Every
    region is finite; ``FAR`` holds every ray, so the abutting grid has no
    gap."""
    grid = (_grid_volume((4, 5, 3)), "extremum_grid", (2, 3, 3))
    radial = (_spherical_volume((4, 4, 4)), "extremum_spherical", (4, 1, 1))
    glob = (_grid_volume((4, 4, 4), seed=3), "extremum_global", None)
    aabb = ((-1.2, -0.3, -2.0), (2.7, 3.1, 1.6))

    def tiled(part, lattice, aabb):
        structure = _extremum(*part)
        cell_min, extents = _bbox(structure)
        lattice = extents if lattice is None else np.array(lattice)
        medium = _repeat(_component(*part), lattice, aabb)
        return medium, (structure, (cell_min, lattice, aabb))

    cases = []
    for name, part, lattice, region in [
        ("grid", grid, None, FAR),
        ("grid_gaps", grid, (1.5, 1.25, 1.0), FAR),
        ("grid_aabb", grid, (1.5, 1.0, 2.0), aabb),
        ("radial_gaps", radial, (2.5, 3.0, 2.0), FAR),
        ("global_gaps", glob, (1.25, 1.5, 1.75), FAR),
    ]:
        medium, entry = tiled(part, lattice, region)
        cases.append((name, medium, [entry], _bbox(entry[0])))

    medium, entry = tiled(grid, (1.5, 1.25, 1.0), aabb)
    cases.append(
        (
            "overlap_repeat",
            _overlap([medium, _component(*glob)]),
            [entry, (_extremum(*glob), None)],
            _bbox(entry[0]),
        )
    )

    lattice = np.array([2.5, 2.25, 2.0])
    structures = [_extremum(*grid), _extremum(*radial)]
    cases.append(
        (
            "repeat_overlap",
            _repeat(_overlap([_component(*grid), _component(*radial)]), lattice),
            [(s, (_bbox(s)[0], lattice, FAR)) for s in structures],
            _bbox(structures[1]),
        )
    )
    return cases


def _repeat_rays(domain):
    """Long rays crossing several cells; axis-aligned ones through a tile,
    through a gap and beside a finite region."""
    tile_min, extents = domain
    rng = np.random.default_rng(4)
    out = []
    for _ in range(8):
        o = tile_min + extents * (rng.random(3) * 3.0 - 1.0)
        d = rng.random(3) * 2.0 - 1.0
        out.append((o, d / np.linalg.norm(d)))
    for axis in range(3):
        d = np.zeros(3)
        d[axis] = 1.0 if axis != 1 else -1.0
        for offset in ([0.3, 0.6, 0.4], [1.1, 1.05, 1.02], [0.3, -0.8, 0.4]):
            out.append((tile_min + extents * np.array(offset), d))
    return [
        mi.Ray3f(o=mi.Point3f(*map(float, o)), d=mi.Vector3f(*map(float, d)))
        for o, d in out
    ]


def _slab(o, d, lo, hi):
    lo, hi = np.array(lo), np.array(hi)
    inside = (lo <= o) & (o <= hi)
    with np.errstate(divide="ignore", invalid="ignore"):
        t0, t1 = (lo - o) / d, (hi - o) / d
    near = np.where(d == 0.0, np.where(inside, -np.inf, np.inf), np.minimum(t0, t1))
    far = np.where(d == 0.0, np.where(inside, np.inf, -np.inf), np.maximum(t0, t1))
    return near.max(), far.min()


def _walk_entry(structure, tiling, ray, mint, maxt):
    """Oracle for one entry: split its part of ``[mint, maxt]`` at the world
    lattice planes, walk the structure over each cell on the ray shifted into
    the tile, and fill each stretch without a segment with one zero segment."""
    o = np.array(ray.o).reshape(-1).astype(np.float64)
    d = np.array(ray.d).reshape(-1).astype(np.float64)
    ts, cell_min, lattice = {mint, maxt}, np.zeros(3), np.full(3, np.inf)
    if tiling is not None:
        cell_min, lattice, region = tiling
        r0, r1 = _slab(o, d, *region)
        start, end = max(mint, r0), min(maxt, r1)
        if start >= end:
            return [(mint, maxt, 0.0, 0.0)]
        ts = {start, end}
        for i in range(3):
            if d[i] == 0.0:
                continue
            n0, n1 = sorted(
                (o[i] + d[i] * t - cell_min[i]) / lattice[i] for t in (start, end)
            )
            for n in range(int(np.ceil(n0)), int(np.floor(n1)) + 1):
                t = (cell_min[i] + n * lattice[i] - o[i]) / d[i]
                if start < t < end:
                    ts.add(float(t))

    out, cursor = [], mint
    for lo, hi in zip(sorted(ts)[:-1], sorted(ts)[1:]):
        cell = np.floor((o + d * 0.5 * (lo + hi) - cell_min) / lattice)
        shifted = o - np.nan_to_num(cell) * np.where(np.isinf(lattice), 0, lattice)
        shifted = mi.Ray3f(o=mi.Point3f(*map(float, shifted)), d=ray.d)
        for segment in _walk_dda(structure, shifted, lo, hi):
            if segment[0] > cursor:
                out.append((cursor, segment[0], 0.0, 0.0))
            out.append(segment)
            cursor = segment[1]
        if cursor < hi:
            out.append((cursor, hi, 0.0, 0.0))
            cursor = hi
    if cursor < maxt:
        out.append((cursor, maxt, 0.0, 0.0))
    return out


def _walk_medium(medium, ray, mint, maxt):
    state = medium.dda_init(ray, mint, maxt)
    out = []
    for _ in range(MAX_STEPS):
        if not _f(state.mint) < _f(state.maxt):
            return out
        segment, state = medium.dda_step(state, ray)
        out.append(
            (
                _f(segment.mint),
                _f(segment.maxt),
                _f(segment.minorant()),
                _f(segment.majorant()),
            )
        )
    pytest.fail("dda traversal did not terminate")


def _t_max(domain):
    return 6.0 * float(domain[1].max())


def _solid(segments, eps=1e-5):
    """Drop the float-noise slivers both routes leave at cell boundaries."""
    return [s for s in segments if s[1] - s[0] > eps]


def test_repeat_segments_tile(variant_scalar_rgb):
    for name, medium, _, domain in _repeat_cases():
        t_max = _t_max(domain)
        for ray in _repeat_rays(domain):
            segments = _walk_medium(medium, ray, 0.0, t_max)
            assert segments[0][0] == 0.0, (name, ray)
            assert segments[-1][1] == t_max, (name, ray)
            for before, after in zip(segments[:-1], segments[1:]):
                assert before[0] <= before[1] == after[0], (name, ray)


def test_repeat_matches_shifted_structures(variant_scalar_rgb):
    gaps = {}
    for name, medium, entries, domain in _repeat_cases():
        t_max = _t_max(domain)
        for ray in _repeat_rays(domain):
            walks = [_walk_entry(s, tiling, ray, 0.0, t_max) for s, tiling in entries]
            expected = _solid(_sum(walks, 0.0, t_max))
            actual = _solid(_walk_medium(medium, ray, 0.0, t_max))
            assert len(expected) == len(actual), (name, ray)
            assert np.allclose(expected, actual, rtol=1e-5, atol=1e-5), (name, ray)
            gaps[name] = gaps.get(name, False) or any(s[3] == 0.0 for s in actual)
    assert gaps == {
        "grid": False,
        "grid_gaps": True,
        "grid_aabb": True,
        "radial_gaps": True,
        "global_gaps": True,
        "overlap_repeat": True,
        "repeat_overlap": True,
    }


def test_repeat_rejects_invalid_tiling(variant_scalar_rgb):
    grid = (_grid_volume((4, 5, 3)), "extremum_grid", (2, 3, 3))
    with pytest.raises(RuntimeError, match="already tiled"):
        _repeat(_repeat(_component(*grid)))
    with pytest.raises(RuntimeError, match="does not fit"):
        _repeat(_component(*grid), (1.5, 0.5, 1.0))
    with pytest.raises(RuntimeError, match="aabb_min"):
        _repeat(_component(*grid), aabb=None)


def test_repeat_wide_matches_scalar(variants_vec_backends_once_rgb):
    """All rays of a case in one wide walk, lanes crossing cells and gaps out
    of step, against the same walks one ray at a time in ``scalar_rgb``."""
    variant = mi.variant()
    mi.set_variant("scalar_rgb")
    expected = {
        name: [
            _solid(_walk_medium(medium, ray, 0.0, _t_max(domain)))
            for ray in _repeat_rays(domain)
        ]
        for name, medium, _, domain in _repeat_cases()
    }
    mi.set_variant(variant)

    for name, medium, _, domain in _repeat_cases():
        rays = _repeat_rays(domain)
        o = np.array([np.array(r.o).reshape(-1) for r in rays]).T
        d = np.array([np.array(r.d).reshape(-1) for r in rays]).T
        ray = mi.Ray3f(o=mi.Point3f(o), d=mi.Vector3f(d))
        n = len(rays)
        state = medium.dda_init(
            ray, dr.zeros(mi.Float, n), dr.full(mi.Float, _t_max(domain), n)
        )
        actual = [[] for _ in rays]
        for _ in range(MAX_STEPS):
            active = state.mint < state.maxt
            if not dr.any(active):
                break
            segment, state = medium.dda_step(state, ray, active)
            columns = [
                np.array(x)
                for x in (
                    segment.mint,
                    segment.maxt,
                    segment.minorant(),
                    segment.majorant(),
                )
            ]
            for lane in np.nonzero(np.array(active))[0]:
                actual[lane].append(tuple(float(c[lane]) for c in columns))
        else:
            pytest.fail("dda traversal did not terminate")

        for e, a in zip(expected[name], actual):
            for before, after in zip(a[:-1], a[1:]):
                assert before[1] == after[0], name
            a = _solid(a)
            assert len(e) == len(a), name
            assert np.allclose(e, a, rtol=1e-5, atol=1e-5), name
