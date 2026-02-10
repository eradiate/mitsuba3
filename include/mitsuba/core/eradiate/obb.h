#pragma once

#include <mitsuba/core/bbox.h>
#include <mitsuba/core/transform.h>

NAMESPACE_BEGIN(mitsuba)

/**
 * \brief Oriented Bounding Box for SAT-based overlap testing
 *
 * Stores an OBB in decomposed form (center, rotation, half-extents)
 * optimized for repeated overlap tests using the Separating Axis Theorem.
 */
template <typename Float_>
struct OrientedBoundingBox {
    using Float       = Float_;
    using Mask        = dr::mask_t<Float>;
    using Point3f     = Point<Float, 3>;
    using Vector3f    = Vector<Float, 3>;
    using Matrix3f    = dr::Matrix<Float, 3>;
    using BBox3f      = BoundingBox<Point3f>;
    using Transform4f = Transform<Point<Float, 4>, true>;

    /// OBB center in world space
    Point3f  center;
    /// Orthonormal 3x3 rotation matrix (columns = OBB local axes in world space)
    Matrix3f rotation;
    /// Half-extents with scale absorbed (local half-extents * per-axis scale)
    Vector3f half_extents;

    /**
     * \brief Construct an OBB from a local-space AABB and a to-world transform
     *
     * Extracts rotation and scale from the transform's 3x3 submatrix.
     * Scale is absorbed into half_extents so the SAT formulas work directly.
     */
    OrientedBoundingBox(const BBox3f &bbox, const Transform4f &to_world) {
        Vector3f local_half = bbox.extents() * Float(0.5);

        center = to_world * Point3f(bbox.center());

        // Extract 3x3 columns, normalize to get rotation, absorb scale
        for (size_t j = 0; j < 3; ++j) {
            Vector3f col(
                to_world.matrix.entry(0, j),
                to_world.matrix.entry(1, j),
                to_world.matrix.entry(2, j)
            );
            Float len = dr::norm(col);
            rotation.entry(0, j) = col.x() / len;
            rotation.entry(1, j) = col.y() / len;
            rotation.entry(2, j) = col.z() / len;
            half_extents.entry(j) = local_half.entry(j) * len;
        }
    }

    /**
     * \brief Test overlap between two oriented bounding boxes
     *
     * General OBB-OBB overlap test using the Separating Axis Theorem.
     * Computes the relative rotation between the two OBBs and tests
     * 15 potential separating axes.
     *
     * \return true if the boxes overlap
     */
    Mask overlaps(const OrientedBoundingBox &other) const {

        // Translation in world space
        Vector3f t_world = other.center - center;
        // Express translation in this box's local frame
        Matrix3f Rt = dr::transpose(rotation);
        Vector3f t = Rt*t_world;

        // Relative rotation: C = R_this^T * R_other
        Matrix3f C = Rt * other.rotation;
        Matrix3f absC = dr::abs(C) + dr::Epsilon<Float>;

        const Vector3f &a = half_extents;
        const Vector3f &b = other.half_extents;
        Mask result(true);
        Float ra, rb;
        Vector3f ra_v, rb_v, proj;

        // Test axes: this box's face normals (axes of A, expressed in A's frame)
        ra_v = a;
        rb_v = absC * b;
        result &= dr::all(dr::abs(t) <= ra_v + rb_v);

        // Test axes: other box's face normals (axes of B, expressed in A's frame)
        ra_v = absC * a;
        rb_v = b;
        proj = C*t;
        result &= dr::all(proj <= ra_v + rb_v);

        // When both OBBs are axis-aligned, C ≈ I and cross-product axes are
        // degenerate (zero or duplicate of face normals). Skip them.
        if constexpr (!dr::is_array_v<Float>) {
            Matrix3f identity = dr::identity<Matrix3f>();
            if ( dr::all_nested( dr::abs(C - identity) < 1e-6f ) )
                return result;
        }

        // Test axes: 9 edge-edge cross products (A_i x B_j)
        // A0 x B0
        ra = a.entry(1) * absC.entry(2, 0) + a.entry(2) * absC.entry(1, 0);
        rb = b.entry(1) * absC.entry(0, 2) + b.entry(2) * absC.entry(0, 1);
        result &= dr::abs(t.entry(2) * C.entry(1, 0) - t.entry(1) * C.entry(2, 0)) <= ra + rb;

        // A0 x B1
        ra = a.entry(1) * absC.entry(2, 1) + a.entry(2) * absC.entry(1, 1);
        rb = b.entry(0) * absC.entry(0, 2) + b.entry(2) * absC.entry(0, 0);
        result &= dr::abs(t.entry(2) * C.entry(1, 1) - t.entry(1) * C.entry(2, 1)) <= ra + rb;

        // A0 x B2
        ra = a.entry(1) * absC.entry(2, 2) + a.entry(2) * absC.entry(1, 2);
        rb = b.entry(0) * absC.entry(0, 1) + b.entry(1) * absC.entry(0, 0);
        result &= dr::abs(t.entry(2) * C.entry(1, 2) - t.entry(1) * C.entry(2, 2)) <= ra + rb;

        // A1 x B0
        ra = a.entry(0) * absC.entry(2, 0) + a.entry(2) * absC.entry(0, 0);
        rb = b.entry(1) * absC.entry(1, 2) + b.entry(2) * absC.entry(1, 1);
        result &= dr::abs(t.entry(0) * C.entry(2, 0) - t.entry(2) * C.entry(0, 0)) <= ra + rb;

        // A1 x B1
        ra = a.entry(0) * absC.entry(2, 1) + a.entry(2) * absC.entry(0, 1);
        rb = b.entry(0) * absC.entry(1, 2) + b.entry(2) * absC.entry(1, 0);
        result &= dr::abs(t.entry(0) * C.entry(2, 1) - t.entry(2) * C.entry(0, 1)) <= ra + rb;

        // A1 x B2
        ra = a.entry(0) * absC.entry(2, 2) + a.entry(2) * absC.entry(0, 2);
        rb = b.entry(0) * absC.entry(1, 1) + b.entry(1) * absC.entry(1, 0);
        result &= dr::abs(t.entry(0) * C.entry(2, 2) - t.entry(2) * C.entry(0, 2)) <= ra + rb;

        // A2 x B0
        ra = a.entry(0) * absC.entry(1, 0) + a.entry(1) * absC.entry(0, 0);
        rb = b.entry(1) * absC.entry(2, 2) + b.entry(2) * absC.entry(2, 1);
        result &= dr::abs(t.entry(1) * C.entry(0, 0) - t.entry(0) * C.entry(1, 0)) <= ra + rb;

        // A2 x B1
        ra = a.entry(0) * absC.entry(1, 1) + a.entry(1) * absC.entry(0, 1);
        rb = b.entry(0) * absC.entry(2, 2) + b.entry(2) * absC.entry(2, 0);
        result &= dr::abs(t.entry(1) * C.entry(0, 1) - t.entry(0) * C.entry(1, 1)) <= ra + rb;

        // A2 x B2
        ra = a.entry(0) * absC.entry(1, 2) + a.entry(1) * absC.entry(0, 2);
        rb = b.entry(0) * absC.entry(2, 1) + b.entry(1) * absC.entry(2, 0);
        result &= dr::abs(t.entry(1) * C.entry(0, 2) - t.entry(0) * C.entry(1, 2)) <= ra + rb;

        return result;
    }

    DRJIT_STRUCT_NODEF(OrientedBoundingBox, center, rotation, half_extents)
};

NAMESPACE_END(mitsuba)
