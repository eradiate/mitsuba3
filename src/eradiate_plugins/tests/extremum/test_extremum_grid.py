import mitsuba as mi
import numpy as np
import pytest


def make_gridvolume(data, filter_type="nearest", to_world=None, wrap_mode="clamp"):
    wrap = wrap_mode != "none"
    wrap_mode = wrap_mode if wrap_mode != "none" else "clamp"

    d = {
        "type": "gridvolume",
        "grid": mi.VolumeGrid(data.astype(np.float32)),
        "filter_type": filter_type,
        "accel": False,
        "wrap_mode": wrap_mode,
        "wrap": wrap,
    }
    if to_world is not None:
        d["to_world"] = to_world
    return mi.load_dict(d)


def build_extremum_grid(volume, resolution, domain=None, scale=1.0):
    struct = mi.load_dict({"type": "extremum_grid", "resolution": list(resolution)})
    if domain is None:
        domain = volume.bbox()
    struct.build(mi.BoundingBox3f(domain.min, domain.max), volume, scale)
    return struct


def extremum_data(struct, resolution):
    data = mi.traverse(struct)["data"].numpy()
    data = data.reshape(resolution[2], resolution[1], resolution[0], 2)
    return data.transpose(2, 1, 0, 3)  # -> [x, y, z, {minorant, majorant}]


def march(struct, ray, t=0.0, max_steps=200):
    """March next_segment until the segment extends to infinity."""
    segments = []
    for _ in range(max_steps):
        seg = struct.next_segment(ray, t)
        assert seg.mint == t
        assert seg.maxt > seg.mint, "no forward progress"
        segments.append((seg.mint, seg.maxt, seg.minorant(), seg.majorant()))
        t = seg.maxt
        if not np.isfinite(t):
            return segments
    pytest.fail("next_segment did not terminate")


# ------------------------------------------------------------------------
# Build: stored bounds
# ------------------------------------------------------------------------


def test_build_high_res(variant_scalar_mono_double):
    n_x, n_y, n_z = 4, 8, 3
    data = np.linspace(1, n_x * n_y * n_z, n_x * n_y * n_z).reshape(n_x, n_y, n_z)
    volume = make_gridvolume(data.transpose(2, 1, 0))

    struct = build_extremum_grid(volume, (n_x, n_y, n_z))
    grid = extremum_data(struct, (n_x, n_y, n_z))

    # One extremum cell per data texel: minorant == majorant == texel value
    assert np.allclose(data, grid[:, :, :, 0])
    assert np.allclose(data, grid[:, :, :, 1])


def test_build_half_res(variant_scalar_mono_double):
    n_x, n_y, n_z = 4, 8, 4
    data = np.linspace(1, n_x * n_y * n_z, n_x * n_y * n_z).reshape(n_x, n_y, n_z)
    volume = make_gridvolume(data.transpose(2, 1, 0))

    struct = build_extremum_grid(volume, (2, 4, 2))
    grid = extremum_data(struct, (2, 4, 2))

    assert np.allclose(data[::2, ::2, ::2], grid[:, :, :, 0])
    assert np.allclose(data[1::2, 1::2, 1::2], grid[:, :, :, 1])


def test_build_scale(variant_scalar_mono_double):
    data = np.ones((2, 2, 2))
    volume = make_gridvolume(data)
    struct = build_extremum_grid(volume, (2, 2, 2), scale=3.0)
    grid = extremum_data(struct, (2, 2, 2))
    assert np.allclose(grid, 3.0)


def test_double_build_guard(variant_scalar_mono_double):
    # Invariant: one structure belongs to one medium — rebuilding with a
    # different domain must fail loudly.
    volume = make_gridvolume(np.ones((2, 2, 2)))
    struct = build_extremum_grid(volume, (2, 2, 2))
    with pytest.raises(RuntimeError, match="different domain"):
        struct.build(mi.BoundingBox3f([0, 0, 0], [2, 2, 2]), volume, 1.0)


# ------------------------------------------------------------------------
# next_segment
# ------------------------------------------------------------------------


def test_next_segment_march(variant_scalar_mono_double):
    # Volume over [0,1]^3 (values varying along x), clamp-extended domain
    # [-1,2]^3. Exercises: empty segment outside the domain, semi-infinite
    # border cells on both sides, exact tiling of interior cells.
    values = np.array([1.0, 2.0, 3.0, 4.0])
    data = values.reshape(1, 1, 4)  # (z, y, x)
    volume = make_gridvolume(data)
    domain = mi.BoundingBox3f([-1, -1, -1], [2, 2, 2])
    struct = build_extremum_grid(volume, (4, 1, 1), domain=domain)

    ray = mi.Ray3f(o=[-2.0, 0.5, 0.5], d=[1.0, 0.0, 0.0])
    segments = march(struct, ray)

    reference = [
        (0.0, 1.0, 0.0, 0.0),  # before the domain
        (1.0, 2.0, 1.0, 1.0),  # border region x in [-1, 0]: clamped to cell 0
        (2.0, 2.25, 1.0, 1.0),  # interior cells
        (2.25, 2.5, 2.0, 2.0),
        (2.5, 2.75, 3.0, 3.0),
        (2.75, 3.0, 4.0, 4.0),
        (3.0, 4.0, 4.0, 4.0),  # border region x in [1, 2]: clamped to cell 3
        (4.0, np.inf, 0.0, 0.0),  # past the domain
    ]
    assert len(segments) == len(reference)
    for seg, ref in zip(segments, reference):
        assert np.allclose(seg, ref, rtol=1e-5)


def test_next_segment_outside_miss(variant_scalar_mono_double):
    volume = make_gridvolume(np.ones((2, 2, 2)))
    struct = build_extremum_grid(volume, (2, 2, 2))

    # Ray that never enters the domain: one empty segment to infinity
    ray = mi.Ray3f(o=[-1.0, 5.0, 0.5], d=[1.0, 0.0, 0.0])
    seg = struct.next_segment(ray, 0.0)
    assert seg.maxt == np.inf
    assert seg.minorant() == 0.0 and seg.majorant() == 0.0


@pytest.mark.parametrize(
    "origin, direction",
    [
        # Along the x=0.5, y=0.5 cell face/edge line
        ([0.5, 0.5, -1.0], [0.0, 0.0, 1.0]),
        # Through cell corners on the main diagonal
        ([-0.5, -0.5, -0.5], [1.0, 1.0, 1.0]),
        # Starting exactly on a face
        ([0.25, 0.5, 0.5], [1.0, 0.0, 0.0]),
    ],
)
def test_next_segment_stall(variant_scalar_mono_double, origin, direction):
    # Regression for the forward-nudge: rays exactly along cell boundaries
    # must make forward progress and terminate (march() asserts both)
    rng = np.random.default_rng(11)
    volume = make_gridvolume(rng.uniform(0.1, 2.0, (4, 4, 4)))
    struct = build_extremum_grid(volume, (4, 4, 4))

    d = np.array(direction, dtype=np.float64)
    ray = mi.Ray3f(o=origin, d=(d / np.linalg.norm(d)).tolist())
    march(struct, ray)
