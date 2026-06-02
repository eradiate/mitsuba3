.. _sec-eradiate-shapes:

Shapes
======

The Eradiate kernel provides two shape plugins:

- ``arectangle`` — an axis-aligned rectangle shape, equivalent to
  Mitsuba's ``rectangle`` but with additional support for Eradiate's
  atmospheric scene geometry conventions.
- ``instancelist`` — efficiently instantiates a large collection of
  identical subsets of a scene geometry (*e.g.* a canopy of identical trees) by
  referencing a shared shape group with many per-instance transforms.
