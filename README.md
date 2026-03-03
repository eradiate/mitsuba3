![Eradiate logo](https://raw.githubusercontent.com/eradiate/eradiate/main/docs/_static/eradiate-logo.svg "Eradiate — A new-generation radiative transfer simulation package")

# Eradiate Radiometric Kernel ­— Mitsuba 3


[![pypi][pypi-badge]][pypi-url]

[pypi-badge]: https://img.shields.io/pypi/v/eradiate-mitsuba?style=flat-square
[pypi-url]: https://pypi.org/project/eradiate-mitsuba/

## Introduction

This is the radiometric kernel of the
[Eradiate radiative transfer model](https://github.com/eradiate/eradiate), based
on [Mitsuba 3](https://github.com/mitsuba-renderer/mitsuba3), the rendering
system developed at [EPFL](https://www.epfl.ch) in Switzerland. This modified
version builds on Mitsuba 3's feature set to provide a modern technical
foundation for the implementation of high-performance Monte Carlo ray tracing
techniques to solve atmospheric radiative transfer problems.

## How does this project relate to Mitsuba 3?

The Eradiate kernel is a fork of Mitsuba 3. It inherits many of its distinctive
features, such as retargetability (the ability to apply compile-time
transformations to the C++ codebase to change fundamental aspects of the
simulation such as spectral representation or polarization), and is powered by
the same just-in-time compiler [Dr.Jit](https://github.com/mitsuba-renderer/drjit).

However, the Eradiate kernel is more specialized than upstream Mitsuba 3:

- **Monochromatic computation**: Currently, Eradiate hands the management of the
  spectral dimension to pre- and post-processing components external to the
  radiometric engine. Consequently, all radiometric computations are
  monochromatic.

- **Double precision**: Earth observation problems involve a wide range of
  length scales: objects of a size of ~1 cm, such as leaves, coexist with
  objects of a size of ~10.000 km, such as planets. In such conditions,
  performing ray tracing safely is challenging. For safety, Eradiate currently
  relies only on double-precision variants, which has two consequences: (1) the
  high-performance BVH supplied by Embree cannot be used (it supports only
  single-precision), and we must fall back to the built-in kd-tree; (2) using
  JIT-compiled variants targeting the GPU is not yet possible. Addressing these
  limitations is planned.

- **Scalar variants only**: For now, only scalar variants are used. Supporting
  JIT compilation is planned.

## About

The Eradiate kernel is maintained by [Rayference](https://www.rayference.eu/).

Mitsuba was created by [Wenzel Jakob](https://rgl.epfl.ch/people/wjakob), and
the Eradiate team thanks all Mitsuba contributors, in particular
[EPFL's Realistic Graphics Lab](https://rgl.epfl.ch/), for their work and support.
