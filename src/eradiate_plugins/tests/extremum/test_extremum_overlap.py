import mitsuba as mi
import numpy as np
import pytest


def grid_struct(value, box_min, box_max, resolution=(2, 2, 2)):
    extents = np.subtract(box_max, box_min)
    volume = mi.load_dict(
        {
            "type": "gridvolume",
            "grid": mi.VolumeGrid(np.full((2, 2, 2), value, dtype=np.float32)),
            "filter_type": "nearest",
            "accel": False,
            "to_world": mi.ScalarTransform4f()
            .translate(list(box_min))
            .scale(extents.tolist()),
        }
    )
    struct = mi.load_dict({"type": "extremum_grid", "resolution": list(resolution)})
    struct.build(mi.BoundingBox3f(list(box_min), list(box_max)), volume, 1.0)
    return struct


def march(struct, ray, t=0.0, max_steps=100):
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


def test_overlap_disjoint_gap(variant_scalar_mono_double):
    # Disjoint children: entries, exits, and the gap between them all fall
    # out of the same earliest-exit minimum; the gap is an empty segment
    a = grid_struct(1.0, [0, 0, 0], [1, 1, 1])
    b = grid_struct(2.0, [2, 0, 0], [3, 1, 1])
    overlap = mi.load_dict({"type": "extremum_overlap", "a": a, "b": b})

    ray = mi.Ray3f(o=[-1.0, 0.5, 0.5], d=[1.0, 0.0, 0.0])
    assert_segments(
        march(overlap, ray),
        [
            (0.0, 1.0, 0.0, 0.0),
            (1.0, 1.5, 1.0, 1.0),  # child a, internal cell boundary at 1.5
            (1.5, 2.0, 1.0, 1.0),
            (2.0, 3.0, 0.0, 0.0),  # gap
            (3.0, 3.5, 2.0, 2.0),  # child b
            (3.5, 4.0, 2.0, 2.0),
            (4.0, np.inf, 0.0, 0.0),
        ],
    )


def test_overlap_summed_bounds(variant_scalar_mono_double):
    # Overlapping region [1, 2]: bounds are the sums of the children's bounds.
    # (Resolution >= 2: a single-cell grid stores the loose (0, max) global
    # pair by design.)
    a = grid_struct(1.0, [0, 0, 0], [2, 1, 1], resolution=(2, 1, 1))
    b = grid_struct(2.0, [1, 0, 0], [3, 1, 1], resolution=(2, 1, 1))
    overlap = mi.load_dict({"type": "extremum_overlap", "a": a, "b": b})

    ray = mi.Ray3f(o=[0.0, 0.5, 0.5], d=[1.0, 0.0, 0.0])
    assert_segments(
        march(overlap, ray),
        [
            (0.0, 1.0, 1.0, 1.0),
            (1.0, 2.0, 3.0, 3.0),
            (2.0, 3.0, 2.0, 2.0),
            (3.0, np.inf, 0.0, 0.0),
        ],
    )
