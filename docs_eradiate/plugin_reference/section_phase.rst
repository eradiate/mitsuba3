.. _sec-eradiate-phase:

Phase functions
===============

The Eradiate fork provides phase function plugins suited for atmospheric
radiative transfer:

- :monosp:`rayleigh_polarized` — Rayleigh scattering phase matrix with full
  polarimetric support, appropriate for molecular scattering in the atmosphere.
- :monosp:`tabphase_irregular` — tabulated phase function defined on an
  irregular angular grid, suitable for aerosol or cloud particle scattering.
- :monosp:`tabphase_polarized` — tabulated polarimetric phase matrix (Müller
  matrix form) for polarized light transport through anisotropic media.
- :monosp:`multiphase` — linearly combines multiple phase functions, enabling
  composite media (e.g., a mixture of molecular and aerosol scattering).
