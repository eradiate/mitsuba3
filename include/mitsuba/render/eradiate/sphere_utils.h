#pragma once

#include <mitsuba/core/fwd.h>
#include <mitsuba/core/platform.h>
#include <mitsuba/core/vector.h>

NAMESPACE_BEGIN(mitsuba)


/**
 * \brief Frame parameters of radially parameterized volumes in world space.
 */
template<typename Float>
struct SphericalParameters {

    using Point3f = Point<Float, 3>;

    Point3f center;
    Float rmin, rmax;

    SphericalParameters():
        center(Point3f(0.f)),
        rmin(dr::Infinity<Float>),
        rmax(0.f) {}

    SphericalParameters(Point3f center, Float rmin, Float rmax):
        center(center),
        rmin(rmin),
        rmax(rmax) {}

    MI_INLINE bool valid() {
        return rmax > rmin;
    }
};

NAMESPACE_END(mitsuba)
