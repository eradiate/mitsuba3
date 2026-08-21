#pragma once

#include <mitsuba/core/fwd.h>
#include <mitsuba/core/platform.h>
#include <mitsuba/core/vector.h>
#include <drjit/texture.h>

NAMESPACE_BEGIN(mitsuba)

/**
 * \brief Volume's coordinate type.
 */
enum class VolumeCoordFlag : uint32_t {
    Grid,
    Spherical
};


/**
 * \brief Frame parameters of radially parameterized volumes in world space.
 */
template<typename Float>
struct VolumeParametrization {

    using Point3f = Point<Float, 3>;
    using BoundingBox3f = BoundingBox<Point3f>;
    using AffineTransform4f = Transform<Point<Float, 4>, true>;

    /// The volume transform
    AffineTransform4f to_world;

    /// The local space range spanned by the volume e.g. dim 0 -> [rmin, rmax] in spherical coords.
    BoundingBox3f uv_range;

    /// The volume coordinate flag, default to grid.
    VolumeCoordFlag flag;

    /// Specifies if the volume wraps.
    bool wrap;

    /// The wrapping behaviour when exiting the volume boundaries.
    dr::WrapMode wrap_mode;

    VolumeParametrization():
        to_world(AffineTransform4f()),
        uv_range(BoundingBox3f(Point3f(0.f), Point3f(1.f))),
        flag(VolumeCoordFlag::Grid),
        wrap(false),
        wrap_mode(dr::WrapMode::Clamp) {}

    VolumeParametrization(
        AffineTransform4f to_world, BoundingBox3f uv_range,
        VolumeCoordFlag flag, bool wrap = false,
        dr::WrapMode wrap_mode = dr::WrapMode::Clamp):
        to_world(to_world),
        uv_range(uv_range),
        flag(flag),
        wrap(wrap),
        wrap_mode(wrap_mode) {}
};

/**
 * \brief Applies the configured texture wrapping mode to an integer
 * position
 */
template <typename T> T wrap(
    const T &pos, const T &res,
    const dr::divisor<dr::scalar_t<T>> *inv_res,
    dr::WrapMode wrap_mode) {
    using Scalar = dr::scalar_t<T>;

    constexpr size_t Dimension = dr::size_v<T>;

    static_assert(
        std::is_integral_v<Scalar> &&
        std::is_signed_v<Scalar>
    );

    if (wrap_mode == dr::WrapMode::Clamp) {
        return clip(pos, 0, res - 1);
    } else {
        T value_shift_neg = select(pos < 0, pos + 1, pos);

        T div;
        for (size_t i = 0; i < Dimension; ++i)
            div[i] = inv_res[i](value_shift_neg[i]);

        T mod = pos - div * res;
        mod[mod < 0] += T(res);

        if (wrap_mode == dr::WrapMode::Mirror)
            // Starting at 0, flip the texture every other repetition
            // (flip when: even number of repetitions in negative direction,
            // or odd number of repetitions in positive direction)
            mod = select(((div & 1) == 0) ^ (pos < 0), mod, res - 1 - mod);

        return mod;
    }
}

NAMESPACE_END(mitsuba)
