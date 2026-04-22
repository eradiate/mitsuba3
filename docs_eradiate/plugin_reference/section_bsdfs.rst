.. _sec-eradiate-bsdfs:

BSDFs
=====

The Eradiate kernel provides a collection of BSDFs (Bidirectional Scattering
Distribution Functions) targeted at Earth observation applications. The following
models are available:

- **Empirical surface models**: ``rpv`` (Rahman-Pinty-Verstraete),
  ``rtls`` (Ross-Thick Li-Sparse), ``hapke``, ``maignan`` —
  semi-empirical bidirectional reflectance distribution functions calibrated
  against satellite or in-situ measurements.
- **Water surface models**: ``ocean_grasp``, ``ocean_legacy``,
  ``ocean_mishchenko`` — ocean BRDF parametrizations accounting for
  wind-driven surface roughness, foam, and subsurface scattering.
- **Measured reflectance**: ``mqdiffuse``, ``measured_mono`` — special-purpose data-driven material models.
- **Material models**: ``bilambertian`` — simple
  diffuse material model supporting reflection and transmission.
- **Utility**: ``selectbsdf`` — conditionally selects among multiple
  BSDFs based on a material index texture.
