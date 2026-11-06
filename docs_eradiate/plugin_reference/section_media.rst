.. _sec-eradiate-media:

Participating media
===================

DDIS
----

The :ref:`eovolpath <integrator-eovolpath>` integrator implements the VROOM
variance reduction methods, which include directional importance sampling
(DDIS). All media, Eradiate's (e.g.
:ref:`eoheterogeneous <medium-eoheterogeneous>`) as well as Mitsuba's (e.g.
``homogeneous``), accept a ``ddis_threshold`` parameter: the probability of
sampling the scattering direction from the medium's DDIS phase function instead
of its phase function. DDIS is active where ``ddis_threshold`` is positive and
the integrator has ``use_ddis`` enabled. The default, 0, disables it.

The DDIS phase function is not specified by the user. It is computed as the
envelope (the maximum) of the medium's phase functions, and is exposed as the
``ddis_phase_function`` scene parameter.

**Nested media.** :ref:`overlap <medium-overlap>` and
:ref:`repeat <medium-repeat>` contain other media and form a tree. Only the
root medium of the tree is used for DDIS:

- A ``ddis_threshold`` set on the root medium is used, and the values of the
  nested media are ignored.
- Otherwise, the root medium uses the largest value among the media it
  contains, at any depth.

**Which media are built.** The scene builds the DDIS phase function of each
root medium it reaches directly: shape interior and exterior media, and sensor
media. It is refreshed when the medium's phase functions change. Nested media
never build one. Media the scene does not reach directly (for instance inside an
``instance`` or ``shapegroup``) have no DDIS phase function and render without
DDIS.
