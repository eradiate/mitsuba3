# Eradiate Mitsuba — Changelog

Changes made to the Eradiate fork of Mitsuba 3 since the latest release
(**v0.4.3**, 2026-03-18).

---

## Unreleased

### New features

- **`eovolpath` integrator**: New volumetric path integrator that replicates and
  extends the upstream `volpath` integrator. Introduces the *DDIS* (Detector Directional Importance Sampling) method for improved variance
  reduction in participating media
  (`src/eradiate_plugins/integrators/eovolpath.cpp`).

- **Extremum plugins** (`extremum_grid`, `extremum_spherical`): New plugin
  type that computes local majorant/minorant bounds for a medium, enabling
  distance sampling using local majorants rather than a global majorant. Two
  implementations are provided:
  - `extremum_grid` — grid-based extremum tracking.
  - `extremum_spherical` — spherical-coordinates-based extremum tracking.

  These plugins expose a Python API and are accompanied by a comprehensive unit
  test suite.

### Improvements

- **`ExtremumSegment` refactor**: Separated extremum traversal algorithms from
  `sample_segment` into dedicated free functions, improving modularity and
  preparing the codebase for future tracking algorithms. `ExtremumSegment` now
  stores a `Vector2f` and exposes `majorant` / `minorant` as getters.
  (`include/mitsuba/render/eradiate/extremum_segment.h`)

- **Optical thickness rename**: Renamed the `tau` symbol to `ot` (optical
  thickness) throughout the extremum plugins and medium code for consistency
  with Eradiate's notation conventions.

- **`volume.h` additions**: Extended `Volume` header with helpers required by
  the extremum plugins.

### Bug fixes

- **`tabphase_polarized`**: Fixed out-of-sync node handling that could produce
  incorrect results for polarized tabulated phase functions.
  (`src/eradiate_plugins/phase/tabphase_polarized.cpp`)
