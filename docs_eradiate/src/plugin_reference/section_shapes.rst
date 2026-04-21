.. _sec-eradiate-shapes:

Shapes
======

The Eradiate fork provides two shape plugins:

- :monosp:`arectangle` — an axis-aligned rectangle shape, equivalent to
  Mitsuba's :monosp:`rectangle` but with additional support for Eradiate's
  atmospheric scene geometry conventions.
- :monosp:`instancelist` — efficiently instantiates a large collection of
  identical sub-scenes (e.g., a canopy of identical trees) by referencing a
  shared shape group with per-instance transforms.

