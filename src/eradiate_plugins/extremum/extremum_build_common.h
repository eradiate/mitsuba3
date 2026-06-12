#pragma once

#include <mitsuba/core/math.h>
#include <mitsuba/render/volume.h>
#include <nanothread/nanothread.h>

NAMESPACE_BEGIN(mitsuba)

/**
 * \brief Build a fine extremum grid from a volume.
 *
 * Constructs a regular grid of the given resolution in the volume's local
 * space, where each cell stores the [minorant, majorant] of the extinction
 * volume over the corresponding region. The returned buffer is interleaved
 * (minorant, majorant) per cell, with X-fastest / Z-slowest linearization.
 *
 * A safety factor of (1 - eps, 1 + eps) is applied to the values to guarantee
 * conservative bounds under floating point rounding.
 *
 * Shared between extremum structure plugins (e.g. ``extremum_grid``,
 * ``extremum_irregular``).
 *
 * \param volume      Extinction volume queried via \ref Volume::extremum()
 * \param resolution  Grid resolution along XYZ
 * \param scale       Scale factor applied to the extremum values
 *
 * \return The interleaved extremum grid of size ``2 * prod(resolution)``
 */
template <typename Float, typename Spectrum>
DynamicBuffer<Float> build_fine_extremum_grid(
    const Volume<Float, Spectrum> *volume,
    const Vector<int32_t, 3> &resolution,
    dr::scalar_t<Float> scale
) {
    MI_IMPORT_CORE_TYPES()
    using FloatStorage = DynamicBuffer<Float>;

    // local space supergrid cell size
    const ScalarVector3f cell_size = dr::rcp(ScalarVector3f(resolution));

    ScalarVector2f safety_factor(1.f - dr::Epsilon<Float>,
                                 1.f + dr::Epsilon<Float>);

    // Allocate extremum grid data
    size_t n = dr::prod(resolution);

    size_t n_threads  = pool_size() + 1;
    size_t grain_size = std::max(n / (4 * n_threads), (size_t) 1);

    FloatStorage extremum_grid = dr::empty<FloatStorage>(n * 2);

    // Early return if using the global majorant.
    if (n == 1) {
        ScalarFloat max = volume->max();
        dr::scatter(extremum_grid,
                    scale * Vector2f(0.f, max) * safety_factor,
                    UInt32(0));
        return extremum_grid;
    }

    if constexpr (!dr::is_jit_v<Float>) {
        auto guard = volume->pin();

        dr::parallel_for(
            dr::blocked_range<size_t>(0, n, grain_size),
            [&](const dr::blocked_range<size_t> &range) {
                // Recover x, y, z from block start (one-time div/mod per
                // block)

                for (auto idx = range.begin(); idx != range.end(); ++idx) {
                    // Store in linear array (Z-slowest, X-fastest)
                    int32_t x = idx % resolution.x();
                    int32_t y = (idx / resolution.x()) % resolution.y();
                    int32_t z = idx / (resolution.x() * resolution.y());

                    ScalarPoint3f cell_min =
                        ScalarVector3f(x, y, z) * cell_size;
                    ScalarPoint3f cell_max = cell_min + cell_size;
                    ScalarBoundingBox3f cell_bounds(
                        cell_min + math::RayEpsilon<Float>,
                        cell_max - math::RayEpsilon<Float>);

                    // Query volume for local extremum, currently assume
                    // local bounds.
                    auto [min, maj] = volume->extremum(cell_bounds);

                    dr::scatter(extremum_grid,
                                scale * Vector2f(min, maj) * safety_factor,
                                UInt32(idx));
                }
            });
    } else {

        UInt32 idx = dr::arange<UInt32>((uint32_t) n);

        UInt32 x = idx % resolution.x() ;
        UInt32 y = (idx / resolution.x())  % resolution.y();
        UInt32 z =  idx / (resolution.x() * resolution.y());

        Point3f cell_min = Vector3f(x, y, z) * cell_size;
        Point3f cell_max = cell_min + cell_size;
        BoundingBox3f cell_bounds(
                        cell_min + math::RayEpsilon<Float>,
                        cell_max - math::RayEpsilon<Float>
                    );

        auto [min, maj] = volume->extremum(cell_bounds);

        dr::scatter(extremum_grid, scale * min * safety_factor.x(), idx*2);
        dr::scatter(extremum_grid, scale * maj * safety_factor.y(), idx*2+1);
        dr::sync_thread();
    }

    return extremum_grid;
}

NAMESPACE_END(mitsuba)
