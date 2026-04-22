.. _sec-eradiate-integrators:

Integrators
===========

The Eradiate kernel provides volumetric path tracing integrators that extend
Mitsuba's ``volpath`` to address Earth observation use cases:

- ``piecewise_volpath`` — handles piecewise-homogeneous plane-parallel
  atmospheres with a fast regular tracking algorithm.
- ``eovolpath`` — a full-featured Earth observation volumetric path tracer. It
  extends ``volpath`` with EO-specific variance reduction methods and tracking
  estimators.

  .. versionadded:: 0.5.0
