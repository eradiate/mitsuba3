# Eradiate Mitsuba — Changelog

Similar to the upstream Mitsuba 3 codebase, this fork does not strictly follow
the [Semantic Versioning](https://semver.org/) convention. However, we also
strive to document breaking API changes in the release notes below.

---

## v0.6.0 (upcoming)

### Extremum Refactoring

- **Renaming**: the `ExtremumStructure` plugin has been renamed to `Extremum`.
  Folders were renamed from `extremum` to `extrema`. Extremum structures plugin
  ids and files were renamed from `extremum_*` to `*extremum`
  (PR [57](https://github.com/eradiate/mitsuba3/pull/57)).
- **Interface**: the `Extremum` interface has been simplified to avoid
  duplication of volumes, transform, and other parameters that are required by
  both the medium and the extrema. A new `build` function has been added that
  is called by the medium. The extrema are undefined before `build` is called
  (PR [41](https://github.com/eradiate/mitsuba3/pull/41)).
- **Traversal**: dropped `traverse_extremum` and replaced it by a iterator-style
  API that exposes `begin` and `next`. Note that those are not DRJIT methods and
  can only be used on scalar pointers. The functions now have Python bindings.
  The medium now also exposes three new virtual methods `sample_real_interaction`,
  `sample_transmittance`, and `sample_transmittance_pdf`, which by default
  implement delta tracking, ratio tracking and residual ratio tracking. Those
  can be consumed by integrators e.g. `eovolpath`. A russian roulette mechanism
  was added to conditionally terminate path with low transmittance and the
  threshold can be controlled with `rr_threshold`
  (PR [52](https://github.com/eradiate/mitsuba3/pull/52)).
- **Overlapping**: support for overlapping media with a new `overlap` medium and
  and traversal support. By default, supports 4 overlapping media at a time
  which can be modified as a compiler option in `MI_MAX_OVERLAPPING_MEDIA`
  (PR [52](https://github.com/eradiate/mitsuba3/pull/41)).
- **Repeat**: support for periodic media with a new `repeat` media that can
  wrap any other media except another repeat
  (PR [53](https://github.com/eradiate/mitsuba3/pull/53)).
- **Tests**: new render tests that compare `volpath` reference of various scenes
  that include heterogeneous media. The reference are added to the
  `data_eradiate` submodule. Add `extremum` to the xml parser
  (PR [54](https://github.com/eradiate/mitsuba3/pull/54)).

### New features

- **Analytical cube**: add a new cube plugin which defines an analytical cube
  intersection. Can improve performance by up to 10% for atmospheric scenes with
  little geometry in scalar mode. Does not work in cuda
  (PR [60](https://github.com/eradiate/mitsuba3/pull/60)).

### Bug fixes

- **DDIS phase**: Clip the envelope range to avoid out of bound errors in CUDA
  (PR [51](https://github.com/eradiate/mitsuba3/pull/51)).
- **Bbox traversal**: fix edge case of ray parallel with the z-axis report false
  positive intersection (PR [51](https://github.com/eradiate/mitsuba3/pull/51)).
- **Extremum spherical**: correct fillmax evaluation position
  (PR [51](https://github.com/eradiate/mitsuba3/pull/51))

### Improvements

- **DDIS refactor**: refactor the DDIS phase function generation mechanism so
  that only media visible from the scene create a DDIS phase function. Note
  that the phase functions are not create if there are no associated scenes.
  Media with nested media can either get the maximum `ddis_threshold` value of
  its nested components or override it with its own value
  (PR [59](https://github.com/eradiate/mitsuba3/pull/59)).
- **Piecewise folding**: refactored the piecewise medium so that it uses the
  new `sample_real_interaction` and `sample_transmittance` methods, removing
  the need for `sample_analytical_interaction` and
  `transmittance_eval_analytical`. `piecewise_volpath` was entirely removed and
  superseded by `eovolpath`
  (PR [58](https://github.com/eradiate/mitsuba3/pull/58)).

## v0.6.0 (3rd November 2026)

### New features

- **Upstream Mitsuba rebase**: The base Mitsuba version is bumped from 3.7.1 to
  3.9.1. The two-version leap bring major upgrades listed in the Mitsuba
  [release notes](https://mitsuba.readthedocs.io/en/stable/release_notes.html).
  This, in particular, unlocks support for Python 3.14 and higher
  (PR [#56](https://github.com/eradiate/mitsuba3/pull/56)).
- **`eoheterogeneous` medium**: A new heterogeneous that replicates and extends
  the upstream `heterogeneous` medium. Introduces a separate bounding box
  parametrization that allows the medium to have a different bounding box to
  the underlying `sigma_t` volume bbox. In turns allows to exploit the wrap
  mechanism of the volume and implement periodic boundaries
  (PR [#32](https://github.com/eradiate/mitsuba3/pull/32)).
- **`particlephase` phase function**: A polarized phase function
  tabulated on an irregular grid of effective radius and effective variance,
  each grid entry carrying its own irregular `cos(theta)` discretization. At
  a query point, one of the 4 bracketing entries is picked, weighted by its
  bilinear blend weight, and sampled directly.
  (PR [#35](https://github.com/eradiate/mitsuba3/pull/35)).
- **`fisheye` sensor**: new circular fisheye camera. A film sample is mapped
  to its distance from the optical centre in image-circle radii, and the selected
  radial projection model turns that normalized radius into an angle `theta`
  from the optical axis. Supports using real calibration data for faithful
  emulation of real lenses.
  (PR [#49](https://github.com/eradiate/mitsuba3/pull/49)).

### Bug fixes

- **`multiphase` phase**: Guard the plugin from evaluating phase function
  with a weight equal to zero at the interaction point. Handle the case where
  all phase weights are zero at the interaction point. Use the blendphase
  parameter format convention.
- **`ocean_legacy` and `ocean_grasp` bsdf**: Add bound checks to avoid nan
  values due to floating point precision errors.
  (PR [#45](https://github.com/eradiate/mitsuba3/pull/45)).

### Improvements

- Improve the *DDIS* (Detector Directional Importance Sampling) variance
  reduction technique by generating DDIS phase functions, envelopes of the
  max phase function used in a medium. Implemented by `piecewise` and
  `eoheterogeneous` media
  (PR [#32](https://github.com/eradiate/mitsuba3/pull/32)).
- Implement prediction-based splitting (PBS) and N-tuple local estimate (NLE)
  variance reduction methods, completing work on VROOM method support in the
  `eovolpath` integrator
  (PR [#33](https://github.com/eradiate/mitsuba3/pull/33)).
- Add an experimental implementation for periodic boundaries of extremum grid
  structures
  (PR [#34](https://github.com/eradiate/mitsuba3/pull/34)).
- Create a new custom scene traversal function which traverses all parents of
  a node, not just the first one. Combine it with Eradiate's traversal function
  (PR [#37](https://github.com/eradiate/mitsuba3/pull/37)).
- Reduce divergence with upstream `microfacet.h` by remove the angle parameter
  and correlated computations. This can simply be achieved by rotating the input
  and output vectors. Optimize the LUT generation in `ocean_legacy`.
  (PR [#47](https://github.com/eradiate/mitsuba3/pull/47)).

## v0.5.0 (9th June 2026)

This release is the first to feature a full list of changes. It also contains
significant upgrades on our specific tooling and documentation, including a new,
Eradiate-specific docs site deployed to Read The Docs.

We introduce a first, experimental implementation of variance reduction methods
targeted toward cloudy scenes. All of these changes are subject to change in the
few next versions.

### New features

- **`eovolpath` integrator**: A new volumetric path tracer that replicates and
  extends the upstream `volpath` integrator. Introduces the *DDIS* (Detector
  Directional Importance Sampling) method for improved variance reduction in
  participating media, as well as residual ratio tracking for transmittance
  estimation during emitter sampling.
  (PR [#25](https://github.com/eradiate/mitsuba3/pull/25),
  [#27](https://github.com/eradiate/mitsuba3/pull/27))

- **Extremum plugin family**: This new plugin family that provides
  majorant/minorant bounds for a medium, enabling distance sampling using local
  majorants rather than a global majorant. Three implementations are provided:

  - `extremum_global` — spatially invariant extremum.
  - `extremum_grid` — grid-based extremum tracking.
  - `extremum_spherical` — spherical-coordinates-based extremum tracking.

  (PR [#24](https://github.com/eradiate/mitsuba3/pull/24))

- **Stokes-moment integrator**: This new integrator plugin combines the `moment`
  and `stokes` integrators into a single implementation which retrieves both the
  second moment and Stokes vector components of radiance samples simultaneously.
  This addresses a limitation that would constrain the nesting `moment` and
  `stokes` to a configuration that would not allow getting the second moment of
  all Stokes vector components (only the intensity).
  (PR [#29](https://github.com/eradiate/mitsuba3/pull/29))

### Improvements

- **`Volume` interface**: Added a `min()` method to the `Volume` interface,
  propagated to `VolumeGrid` and `Texture`, as required by the extremum plugins.
  (PR [#27](https://github.com/eradiate/mitsuba3/pull/27))

- **`Medium` analytical transmittance API**: Renamed `sample_interaction_real`
  to `sample_interaction_analytical` and `eval_transmittance_pdf_real` to
  `eval_analytical_transmittance`, which now return the full unpolarized
  spectrum. Reworked the `piecewise` medium to compute transmittance and PDF
  correctly in spectral mode, with new spectral unit tests.
  (PR [#31](https://github.com/eradiate/mitsuba3/pull/31))

### Bug fixes

- **`tabphase_polarized`**: Fixed out-of-sync node handling that could produce
  incorrect results for polarized tabulated phase functions.
  (PR [#26](https://github.com/eradiate/mitsuba3/pull/26))

- **`piecewise_volpath`**: Fixed the integrator to work with non-spectral media.
  (PR [#31](https://github.com/eradiate/mitsuba3/pull/31))
