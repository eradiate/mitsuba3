import mitsuba as mi
import numpy as np
import pytest


def make_spherical_volume(rmin=0.5, rmax=1.0, fillmin=1.0, fillmax=0.0):
    # Four radial layers with values 0.4 (innermost) .. 0.1 (outermost)
    data = np.ascontiguousarray(
        (np.ones((4, 1, 1)) * np.linspace(1, 4, 4).reshape(-1, 1, 1) * 0.1)[::-1]
    ).transpose(2, 1, 0)
    return mi.load_dict(
        {
            "type": "sphericalcoordsvolume",
            "volume": {
                "type": "gridvolume",
                "grid": mi.VolumeGrid(data.copy()),
                "filter_type": "nearest",
                "accel": False,
            },
            "rmin": rmin,
            "rmax": rmax,
            "fillmin": fillmin,
            "fillmax": fillmax,
        }
    )


def build_spherical_struct(volume, resolution=(2, 1, 1)):
    struct = mi.load_dict(
        {"type": "extremum_spherical", "resolution": list(resolution)}
    )
    struct.build(volume.bbox(), volume, 1.0)
    return struct


def march(struct, ray, t=0.0, max_steps=50):
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


def assert_segments(segments, reference):
    assert len(segments) == len(reference)
    for seg, ref in zip(segments, reference):
        assert np.allclose(seg, ref, rtol=1e-5, atol=1e-6), f"{seg} != {ref}"


def test_gate_non_spherical_volume(variant_scalar_mono_double):
    # The gate guards correctness, not just tightness: extremum_local axis
    # meanings depend on the spherical frame
    volume = mi.load_dict(
        {
            "type": "gridvolume",
            "grid": mi.VolumeGrid(np.ones((2, 2, 2), dtype=np.float32)),
            "accel": False,
        }
    )
    struct = mi.load_dict({"type": "extremum_spherical", "resolution": [2, 1, 1]})
    with pytest.raises(RuntimeError, match="spherical frame"):
        struct.build(volume.bbox(), volume, 1.0)


def test_radial_march(variant_scalar_mono_double):
    # Downward ray through the sphere center: shell crossings at r = 1.0,
    # 0.75, 0.5 (two shells over [rmin, rmax], reversed data). The [0, rmin]
    # fill region includes fillmin = 1.0 in its majorant.
    struct = build_spherical_struct(make_spherical_volume())
    ray = mi.Ray3f(o=[0, 0, 1.0], d=[0, 0, -1.0])
    assert_segments(
        march(struct, ray),
        [
            (0.0, 0.25, 0.1, 0.2),
            (0.25, 0.5, 0.3, 0.4),
            (0.5, 1.5, 0.4, 1.0),
            (1.5, 1.75, 0.3, 0.4),
            (1.75, 2.0, 0.1, 0.2),
            (2.0, np.inf, 0.0, 0.0),
        ],
    )


def test_radial_march_origin_outside(variant_scalar_mono_double):
    # Ray origin outside the domain: leading empty segment up to the domain
    struct = build_spherical_struct(make_spherical_volume())
    ray = mi.Ray3f(o=[0, 0, 1.5], d=[0, 0, -1.0])
    assert_segments(
        march(struct, ray),
        [
            (0.0, 0.5, 0.0, 0.0),
            (0.5, 0.75, 0.1, 0.2),
            (0.75, 1.0, 0.3, 0.4),
            (1.0, 2.0, 0.4, 1.0),
            (2.0, 2.25, 0.3, 0.4),
            (2.25, 2.5, 0.1, 0.2),
            (2.5, np.inf, 0.0, 0.0),
        ],
    )


def test_radial_march_tangent(variant_scalar_mono_double):
    # Impact parameter 0.75: the ray is tangent to the inner shell boundary
    # at closest approach (t = 1); only the outer shell is traversed. The
    # corner region outside rmax carries the conservative fill majorant.
    struct = build_spherical_struct(make_spherical_volume())
    ray = mi.Ray3f(o=[0.75, 0, 1.0], d=[0, 0, -1.0])
    t_enter = 1.0 - np.sqrt(1.0 - 0.75**2)
    assert_segments(
        march(struct, ray),
        [
            (0.0, t_enter, 0.0, 0.1),
            (t_enter, 1.0, 0.1, 0.2),
            (1.0, 2.0 - t_enter, 0.1, 0.2),
            (2.0 - t_enter, 2.0, 0.0, 0.1),
            (2.0, np.inf, 0.0, 0.0),
        ],
    )
