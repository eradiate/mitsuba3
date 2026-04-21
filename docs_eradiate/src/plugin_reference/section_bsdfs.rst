.. _sec-eradiate-bsdfs:

BSDFs
=====

The Eradiate fork provides a collection of BSDFs (Bidirectional Scattering
Distribution Functions) tailored for Earth observation applications. Unlike the
generic material models in upstream Mitsuba 3, these are one-sided analytical
reflectance models designed for large-scale atmospheric and surface radiative
transfer simulations.

The following models are available:

- **Lambertian variants**: :monosp:`bilambertian`, :monosp:`mqdiffuse` — simple
  diffuse models with bidirectional or anisotropic corrections.
- **Empirical surface models**: :monosp:`rpv` (Rahman–Pinty–Verstraete),
  :monosp:`rtls` (Ross-Thick Li-Sparse), :monosp:`hapke`, :monosp:`maignan` —
  semi-empirical bidirectional reflectance distribution functions calibrated
  against satellite or in-situ measurements.
- **Ocean surface models**: :monosp:`ocean_grasp`, :monosp:`ocean_legacy`,
  :monosp:`ocean_mishchenko` — ocean BRDF parameterizations accounting for
  wind-driven surface roughness, foam, and subsurface scattering.
- **Measured reflectance**: :monosp:`measured_mono` — tabulated reflectance
  data for monochromatic simulations.
- **Utility**: :monosp:`selectbsdf` — conditionally selects among multiple
  BSDFs based on a scene parameter.

