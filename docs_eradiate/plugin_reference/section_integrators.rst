.. _sec-eradiate-integrators:

Integrators
===========

The Eradiate fork provides two volumetric path integrators that extend
Mitsuba's :monosp:`volpath` to address Earth observation use cases:

- :monosp:`piecewise_volpath` — handles piecewise-homogeneous layered
  atmospheres by sampling free paths within each homogeneous layer separately.
- :monosp:`eovolpath` — a full-featured Earth observation volumetric path
  integrator. It replicates and extends :monosp:`volpath` with the *DDIS*
  (Decomposition Direct Importance Sampling) phase-sampling method, which
  decomposes the phase function into a direct and diffuse component to reduce
  variance in anisotropic scattering regimes.

.. note::

   :monosp:`eovolpath` was introduced after v0.4.3.
