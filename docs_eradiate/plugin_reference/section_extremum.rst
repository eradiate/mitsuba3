.. _sec-eradiate-extremum:

Extremum structures
===================

Extremum structures are spatial data structures that store local majorant and
minorant bounds of a heterogeneous medium's extinction coefficient. They are
used to implement *delta tracking* with locally adaptive majorants, which
significantly reduces wasted null-collision samples compared to using a global
majorant.

The abstract base class :py:class:`mitsuba.ExtremumStructure` defines the
interface. Two concrete implementations are provided:

- :monosp:`extremum_grid` — stores bounds on a regular 3-D grid aligned with
  the medium's bounding box.
- :monosp:`extremum_spherical` — stores bounds in a spherical-coordinate grid,
  suitable for planetary-scale atmospheric media.

Both plugins expose a Python API and can be traversed with
:py:func:`mitsuba.traverse`. The associated :py:class:`mitsuba.ExtremumSegment`
struct represents a single traversal segment and is documented in the
:doc:`../api_reference/index`.

.. note::

   This plugin category was introduced after v0.4.3 and is not present in the
   corresponding PyPI release of ``eradiate-mitsuba``.
