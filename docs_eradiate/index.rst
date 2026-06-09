:hide-toc:
:layout: landing

Eradiate Mitsuba
================

**Date**: |today| |
**Version**: |release|

This is the documentation for the Eradiate kernel, the radiometric core of the
`Eradiate radiative transfer model <https://www.eradiate.eu/>`__, based on the
`Mitsuba 3 <https://github.com/mitsuba-renderer/mitsuba3>`__ renderer. It
provides extensions tailored for Earth observation radiative transfer
simulations.

.. grid:: 1 1 2 3
   :gutter: 2
   :padding: 0

   .. grid-item-card:: :iconify:`material-symbols:book-2 height=1.5em` User Guide
      :link: user_guide/index
      :link-type: doc

      Installation, basic usage, and worked examples.

   .. grid-item-card:: :iconify:`mdi:puzzle height=1.5em` Plugin Reference
      :link: plugin_reference/index
      :link-type: doc

      Documentation for all Eradiate-specific Mitsuba plugins.

   .. grid-item-card:: :iconify:`material-symbols:description height=1.5em` API Reference
      :link: api_reference/index
      :link-type: doc

      Python API for Eradiate-specific classes and utilities.

   .. grid-item-card:: :iconify:`material-symbols:code height=1.5em` Developer Guide
      :link: dev_guide/index
      :link-type: doc

      Build instructions, contribution guidelines, and internals.

   .. grid-item-card:: :iconify:`mdi:clock height=1.5em` Changelog
      :link: changelog
      :link-type: doc

      Release history and migration notes.

.. toctree::
   :hidden:
   :caption: User Guide

   user_guide/index

.. toctree::
   :hidden:
   :caption: Reference

   plugin_reference/index
   api_reference/index

.. toctree::
   :hidden:
   :caption: Developer Guide

   dev_guide/index

.. toctree::
   :hidden:
   :caption: About

   changelog
   bibliography
