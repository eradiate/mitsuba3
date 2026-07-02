import drjit as dr
import mitsuba as mi
import numpy as np


def generate_extremum_irregular2(
    volume_grid,
    resolution,
    filter_type="nearest",
    transform=None,
    **kwargs,
):
    if transform is None:
        transform = mi.ScalarAffineTransform4f()

    volume = mi.load_dict(
        {
            "type": "gridvolume",
            "grid": volume_grid,
            "filter_type": filter_type,
            "accel": False,
            "to_world": transform,
        }
    )
    extremum_struct = mi.load_dict(
        {
            "type": "extremum_irregular2",
            "volume": volume,
            "resolution": resolution,
            **kwargs,
        }
    )
    return volume, extremum_struct


def get_boxes(extremum_struct):
    """Return (box_lo, box_hi, box_extrema, index_grid) as numpy arrays."""
    params = mi.traverse(extremum_struct)
    n_boxes = len(params["box_extrema"]) // 2
    return (
        params["box_lo"].numpy().reshape(n_boxes, 3).copy(),
        params["box_hi"].numpy().reshape(n_boxes, 3).copy(),
        params["box_extrema"].numpy().reshape(n_boxes, 2).copy(),
        params["index_grid"].numpy().copy(),
    )


def make_step_volume_x(n_x=8, n_y=4, n_z=4, low=0.2, high=1.0):
    """Step function along x: low for the first half, high for the second."""
    return np.where(np.arange(n_x)[:, None, None] < n_x // 2, low, high) * np.ones(
        (n_x, n_y, n_z)
    )


def assert_compare_segment(ref, other):
    assert dr.allclose(ref.mint, other.mint)
    assert dr.allclose(ref.maxt, other.maxt)
    assert dr.allclose(ref.minorant(), other.minorant())
    assert dr.allclose(ref.majorant(), other.majorant())


# ------------------------------------------------------------------------
# Build tests
# ------------------------------------------------------------------------


def test_build_homogeneous(variant_scalar_rgb):
    # A homogeneous field has a zero extremum span, so the null-collision cost
    # term vanishes and the surface-area term drives every cell into one box.
    data = np.full((4, 6, 8), 0.5)
    _, extremum_struct = generate_extremum_irregular2(
        mi.VolumeGrid(data.transpose(2, 1, 0)), mi.ScalarVector3i(4, 6, 8)
    )
    lo, hi, extrema, index = get_boxes(extremum_struct)

    assert lo.shape[0] == 1
    assert (lo[0] == [0, 0, 0]).all() and (hi[0] == [4, 6, 8]).all()
    assert np.allclose(extrema[0], [0.5, 0.5], rtol=1e-5)
    assert (index == 0).all()


def test_build_step_null_dominated(variant_scalar_rgb):
    # A large cn makes the null-collision term dominate: merging across the
    # x = 4 discontinuity (span 0.8) is rejected, so no box straddles it.
    data = make_step_volume_x()
    _, extremum_struct = generate_extremum_irregular2(
        mi.VolumeGrid(data.transpose(2, 1, 0)),
        mi.ScalarVector3i(8, 4, 4),
        cn=1e6,
    )
    lo, hi, extrema, _ = get_boxes(extremum_struct)
    for l_, h in zip(lo, hi):
        assert h[0] <= 4 or l_[0] >= 4


def test_build_step_traversal_dominated(variant_scalar_rgb):
    # With cn = 0 only the surface-area term remains, which always decreases on
    # merge, so the whole step volume collapses into a single box.
    data = make_step_volume_x()
    _, extremum_struct = generate_extremum_irregular2(
        mi.VolumeGrid(data.transpose(2, 1, 0)),
        mi.ScalarVector3i(8, 4, 4),
        cn=0.0,
    )
    lo, hi, extrema, _ = get_boxes(extremum_struct)
    assert lo.shape[0] == 1
    assert (lo[0] == [0, 0, 0]).all() and (hi[0] == [8, 4, 4]).all()
    assert np.allclose(extrema[0], [0.2, 1.0], rtol=1e-4)


def test_build_partition_invariants(variant_scalar_rgb):
    rng = np.random.default_rng(1)
    n_x, n_y, n_z = 8, 6, 4
    data = rng.uniform(0.0, 1.0, size=(n_x, n_y, n_z))
    data[data < 0.3] = 0.0  # include empty regions
    _, extremum_struct = generate_extremum_irregular2(
        mi.VolumeGrid(data.transpose(2, 1, 0)),
        mi.ScalarVector3i(n_x, n_y, n_z),
    )
    lo, hi, extrema, index = get_boxes(extremum_struct)

    # Every cell is covered by exactly one box
    coverage = np.zeros((n_x, n_y, n_z), dtype=int)
    for l_, h in zip(lo, hi):
        coverage[l_[0] : h[0], l_[1] : h[1], l_[2] : h[2]] += 1
    assert (coverage == 1).all()

    # The cell index matches the box bounds (x-fastest linearization)
    index = index.reshape(n_z, n_y, n_x).transpose(2, 1, 0)
    for b, (l_, h) in enumerate(zip(lo, hi)):
        assert (index[l_[0] : h[0], l_[1] : h[1], l_[2] : h[2]] == b).all()

    # Box extrema bound the covered data
    for b, (l_, h) in enumerate(zip(lo, hi)):
        block = data[l_[0] : h[0], l_[1] : h[1], l_[2] : h[2]]
        assert extrema[b, 0] <= block.min() + 1e-6
        assert extrema[b, 1] >= block.max() - 1e-6


def test_build_scaled(variant_scalar_rgb):
    data = np.full((4, 4, 4), 0.5)
    _, extremum_struct = generate_extremum_irregular2(
        mi.VolumeGrid(data.transpose(2, 1, 0)), mi.ScalarVector3i(4, 4, 4), scale=3.0
    )
    _, _, extrema, _ = get_boxes(extremum_struct)
    assert np.allclose(extrema[0], [1.5, 1.5], rtol=1e-5)


# ------------------------------------------------------------------------
# Point evaluation
# ------------------------------------------------------------------------


def test_eval_1(variants_any_scalar, variants_any_llvm):
    data = make_step_volume_x()
    _, extremum_struct = generate_extremum_irregular2(
        mi.VolumeGrid(data.transpose(2, 1, 0)),
        mi.ScalarVector3i(8, 4, 4),
        cn=1e6,
    )

    for p, expected in [
        ((0.25, 0.5, 0.5), 0.2),
        ((0.75, 0.5, 0.5), 1.0),
        # Out-of-bounds points exercise the clamp policy
        ((-0.5, 0.5, 0.5), 0.2),
        ((1.5, 0.5, 0.5), 1.0),
    ]:
        it = dr.zeros(mi.Interaction3f)
        it.p = mi.Point3f(*p)
        minorant, majorant = extremum_struct.eval_1(it)
        assert dr.allclose(minorant, expected, rtol=1e-4)
        assert dr.allclose(majorant, expected, rtol=1e-4)


# ------------------------------------------------------------------------
# Traversal tests
# ------------------------------------------------------------------------


def test_sample_homogeneous(variants_any_scalar, variants_any_llvm):
    # A homogeneous volume merges into a single box: one segment spans the
    # entire traversal interval
    data = np.full((4, 4, 4), 0.5)
    _, extremum_struct = generate_extremum_irregular2(
        mi.VolumeGrid(data.transpose(2, 1, 0)), mi.ScalarVector3i(4, 4, 4)
    )

    ray = mi.Ray3f(
        o=mi.ScalarVector3f(0.0, 0.5, 0.5),
        d=mi.ScalarVector3f(1.0, 0.0, 0.0),
    )
    res, ot_acc = extremum_struct.sample_segment(ray, 0.0, 1.0, 0.3, True)
    ref_segment = mi.ExtremumSegment(mint=0.0, maxt=1.0, majorant=0.5, minorant=0.5)
    assert dr.allclose(ot_acc, 0.0)
    assert_compare_segment(ref_segment, res)

    # Target beyond the total optical thickness: no valid segment
    res, ot_acc = extremum_struct.sample_segment(ray, 0.0, 1.0, 0.8, True)
    assert not dr.any(res.valid())
    assert dr.allclose(ot_acc, 0.5, rtol=1e-5)


def test_sample_parity_with_grid(variants_any_scalar, variants_any_llvm):
    """Cross-check against extremum_grid on random rays.

    A large cn blocks every merge between cells with distinct values, so the
    irregular grid degenerates to per-cell boxes and both structures bound the
    field with identical piecewise extrema. Their segmentations may still
    differ (the irregular grid widens boundary boxes in clamp mode), so the
    comparison uses segmentation-invariant quantities: target validity and the
    derived distance at which the target optical thickness is reached.
    """
    rng = np.random.default_rng(42)
    n_x, n_y, n_z = 8, 6, 4
    data = rng.uniform(0.05, 1.0, size=(n_x, n_y, n_z))
    resolution = mi.ScalarVector3i(n_x, n_y, n_z)

    for wrap_mode in ["clamp", "repeat", "mirror"]:
        volume_grid = mi.VolumeGrid(data.transpose(2, 1, 0))
        volume = mi.load_dict(
            {
                "type": "gridvolume",
                "grid": volume_grid,
                "filter_type": "nearest",
                "accel": False,
            }
        )
        grid_struct = mi.load_dict(
            {
                "type": "extremum_grid",
                "volume": volume,
                "resolution": resolution,
                "wrap_mode": wrap_mode,
            }
        )
        _, irregular_struct = generate_extremum_irregular2(
            volume_grid, resolution, cn=1e6, wrap_mode=wrap_mode
        )

        # No merge occurs across distinct-valued cells: per-cell boxes
        lo, _, _, _ = get_boxes(irregular_struct)
        assert lo.shape[0] == n_x * n_y * n_z

        for trial in range(100):
            o = rng.uniform(-0.5, 1.5, 3)
            d = rng.normal(size=3)
            if trial % 7 == 0:
                d[trial % 3] = 0.0  # direction with a zero component
            if trial % 11 == 0:
                o = np.round(o * n_x) / n_x  # origin exactly on a cell face
            norm = np.linalg.norm(d)
            if norm == 0.0:
                continue
            d = d / norm

            ray = mi.Ray3f(
                mi.Point3f(*[float(v) for v in o]),
                mi.Vector3f(*[float(v) for v in d]),
            )
            maxt = float(rng.uniform(0.5, 3.0))
            tau = float(rng.uniform(0.01, 1.5))

            seg_ref, ot_ref = grid_struct.sample_segment(ray, 0.0, maxt, tau, True)
            seg, ot = irregular_struct.sample_segment(ray, 0.0, maxt, tau, True)

            valid_ref, valid = dr.all(seg_ref.valid()), dr.all(seg.valid())
            assert valid_ref == valid, (wrap_mode, trial)
            if valid_ref:
                # Distance at which the target optical thickness is reached
                ts_ref = seg_ref.mint + (tau - ot_ref) / seg_ref.majorant()
                ts = seg.mint + (tau - ot) / seg.majorant()
                assert dr.allclose(ts_ref, ts, atol=1e-5), (wrap_mode, trial)
                assert dr.allclose(seg_ref.minorant(), seg.minorant())
                assert dr.allclose(seg_ref.majorant(), seg.majorant())
                # The grid cell segment is contained in the irregular box
                eps = 1e-4
                assert dr.all(seg.mint <= seg_ref.mint + eps)
                assert dr.all(seg.maxt >= seg_ref.maxt - eps)
            else:
                # Total accumulated optical thickness must agree
                assert dr.allclose(ot_ref, ot, atol=1e-5), (wrap_mode, trial)


# ------------------------------------------------------------------------
# Volume update tests
# ------------------------------------------------------------------------


def test_refit_keeps_topology(variants_any_scalar, variants_any_llvm):
    data = make_step_volume_x()
    volume, extremum_struct = generate_extremum_irregular2(
        mi.VolumeGrid(data.transpose(2, 1, 0)),
        mi.ScalarVector3i(8, 4, 4),
        cn=1e6,
        rebuild_threshold=-1.0,
    )
    lo0, hi0, _, _ = get_boxes(extremum_struct)

    # Rescale the volume data: the refit must update the extrema in place
    params = mi.traverse(volume)
    new_data = make_step_volume_x(low=0.4, high=2.0)
    params["data"] = mi.TensorXf(
        new_data.transpose(2, 1, 0)[..., None].astype(np.float32)
    )
    params.update()

    lo1, hi1, extrema1, _ = get_boxes(extremum_struct)
    assert (lo0 == lo1).all() and (hi0 == hi1).all()
    # The refit bounds the new field: minorant 0.4 and majorant 2.0 appear
    assert np.isclose(extrema1[:, 0].min(), 0.4, rtol=1e-4)
    assert np.isclose(extrema1[:, 1].max(), 2.0, rtol=1e-4)
