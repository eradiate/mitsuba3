.. _sec-eradiate-sensors:

Sensors
=======

The Eradiate kernel provides a suite of distant measurement sensors designed for
top-of-atmosphere and surface-level radiometric measurements.

- ``hdistant`` — hemispherical distant sensor collecting radiance over a
  hemisphere, mapping zenith and azimuth sampling coordinates to film coordinates.
- ``mdistant`` — multi-directional distant sensor recording radiance at a
  discrete set of viewing directions simultaneously.
- ``mpdistant`` — multi-pixel variant of the ``distant`` plugin that maps the
  target shape's texture coordinates to film coordinates.
- ``mradiancemeter`` — a convenience plugin that groups together an arbitrary
  number of radiancemeter sensors.
- ``distantflux`` — integrates the upwelling radiant flux (irradiance) over the
  entire hemisphere above the scene, mapping zenith and azimuth sampling
  coordinates to film coordinates.
