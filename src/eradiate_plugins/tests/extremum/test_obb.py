import mitsuba as mi
import pytest


def make_obb(bbox_min, bbox_max, to_world):
    bbox = mi.BoundingBox3f(
        mi.Point3f(*bbox_min),
        mi.Point3f(*bbox_max),
    )
    return mi.OrientedBoundingBox3f(bbox, to_world)

def test_identical(variant_scalar_rgb, variants_any_llvm):
    """Same box should overlap itself"""
    trafo = mi.Transform4f.rotate([0, 0, 1], 30)
    obb = make_obb([-1, -1, -1], [1, 1, 1], trafo)
    assert obb.overlaps(obb)

def test_identity_both(variant_scalar_rgb, variants_any_llvm):
    """Two identity OBBs"""
    obb_a = make_obb([0, 0, 0], [2, 2, 2], mi.Transform4f())
    obb_b = make_obb([1, 1, 1], [3, 3, 3], mi.Transform4f())
    assert obb_a.overlaps(obb_b)

def test_identity_separated(variant_scalar_rgb, variants_any_llvm):
    """Two identity OBBs, separated"""
    obb_a = make_obb([0, 0, 0], [1, 1, 1], mi.Transform4f())
    obb_b = make_obb([2, 2, 2], [3, 3, 3], mi.Transform4f())
    assert not obb_a.overlaps(obb_b)

def test_rotated_overlapping(variant_scalar_rgb, variants_any_llvm):
    """Two rotated boxes that overlap"""
    trafo_a = mi.Transform4f.rotate([0, 0, 1], 30)
    trafo_b = mi.Transform4f.rotate([0, 0, 1], -30)
    obb_a = make_obb([-1, -1, -1], [1, 1, 1], trafo_a)
    obb_b = make_obb([-1, -1, -1], [1, 1, 1], trafo_b)
    assert obb_a.overlaps(obb_b)

def test_rotated_separated(variant_scalar_rgb):
    """Two rotated boxes that are separated"""
    trafo_a = mi.Transform4f.rotate([0, 0, 1], 30)
    trafo_b = (
        mi.Transform4f.translate([5, 0, 0])
        @ mi.Transform4f.rotate([0, 0, 1], -30)
    )
    obb_a = make_obb([-1, -1, -1], [1, 1, 1], trafo_a)
    obb_b = make_obb([-1, -1, -1], [1, 1, 1], trafo_b)
    assert not obb_a.overlaps(obb_b)

def test_symmetry(variant_scalar_rgb):
    """a.overlaps(b) == b.overlaps(a)"""
    trafo_a = (
        mi.Transform4f.translate([1, 0, 0])
        @ mi.Transform4f.rotate([0, 1, 0], 45)
    )
    trafo_b = (
        mi.Transform4f.translate([2, 0.5, 0])
        @ mi.Transform4f.rotate([1, 0, 0], 60)
    )
    obb_a = make_obb([-1, -1, -1], [1, 1, 1], trafo_a)
    obb_b = make_obb([-0.5, -0.5, -0.5], [0.5, 0.5, 0.5], trafo_b)
    assert obb_a.overlaps(obb_b) == obb_b.overlaps(obb_a)

def test_symmetry_separated(variant_scalar_rgb):
    """Symmetry holds for separated boxes too"""
    trafo_a = mi.Transform4f.rotate([1, 1, 0], 37)
    trafo_b = (
        mi.Transform4f.translate([10, 0, 0])
        @ mi.Transform4f.rotate([0, 0, 1], 73)
    )
    obb_a = make_obb([-1, -1, -1], [1, 1, 1], trafo_a)
    obb_b = make_obb([-1, -1, -1], [1, 1, 1], trafo_b)
    assert obb_a.overlaps(obb_b) == obb_b.overlaps(obb_a)
    assert not obb_a.overlaps(obb_b)

def test_arbitrary_rotation_known_overlap(variant_scalar_rgb):
    """Arbitrary rotations with geometrically verified overlap"""
    # Two unit boxes at origin, differently rotated: must overlap
    trafo_a = mi.Transform4f.rotate([1, 0, 0], 20)
    trafo_b = mi.Transform4f.rotate([0, 1, 0], 35)
    obb_a = make_obb([-1, -1, -1], [1, 1, 1], trafo_a)
    obb_b = make_obb([-1, -1, -1], [1, 1, 1], trafo_b)
    assert obb_a.overlaps(obb_b)

def test_non_uniform_scale(variant_scalar_rgb):
    """OBB-OBB with non-uniform scale"""
    trafo_a = mi.Transform4f.scale([3, 1, 1])
    trafo_b = (
        mi.Transform4f.translate([4, 0, 0])
        @ mi.Transform4f.scale([1, 3, 1])
    )
    obb_a = make_obb([-1, -1, -1], [1, 1, 1], trafo_a)
    obb_b = make_obb([-1, -1, -1], [1, 1, 1], trafo_b)
    # obb_a extends [-3,3] in X, obb_b centered at 4, extends [3,5] in X
    assert obb_a.overlaps(obb_b)

def test_non_uniform_scale_separated(variant_scalar_rgb):
    """OBB-OBB with non-uniform scale, separated"""
    trafo_a = mi.Transform4f.scale([3, 1, 1])
    trafo_b = (
        mi.Transform4f.translate([5, 0, 0])
        @ mi.Transform4f.scale([1, 3, 1])
    )
    obb_a = make_obb([-1, -1, -1], [1, 1, 1], trafo_a)
    obb_b = make_obb([-1, -1, -1], [1, 1, 1], trafo_b)
    # obb_a extends [-3,3] in X, obb_b centered at 5, extends [4,6] in X
    assert not obb_a.overlaps(obb_b)

