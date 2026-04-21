.. _sec-eradiate-sensors:

Sensors
=======

The Eradiate fork provides a suite of distant measurement sensors designed for
top-of-atmosphere and surface-level radiometric measurements. Unlike Mitsuba's
perspective or orthographic cameras, these sensors record upward-leaving or
downward-arriving radiance integrated over specified angular configurations.

- :monosp:`hdistant` — hemispherical distant sensor collecting radiance over a
  hemisphere, parameterized by zenith and azimuth sampling.
- :monosp:`mdistant` — multi-directional distant sensor recording radiance at a
  discrete set of viewing directions simultaneously.
- :monosp:`mpdistant` — polarimetric variant of :monosp:`mdistant` that records
  full Stokes vector radiance.
- :monosp:`mradiancemeter` — a distant radiance meter measuring the specific
  intensity along a single fixed direction.
- :monosp:`distantflux` — integrates the upwelling radiant flux (irradiance)
  over the entire hemisphere above the scene.

