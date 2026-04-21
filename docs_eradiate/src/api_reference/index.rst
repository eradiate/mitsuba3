API Reference
=============

This page documents the Python-facing API introduced by the Eradiate fork of
Mitsuba 3. Only additions specific to the fork are covered here; the base
Mitsuba 3 API is documented at https://mitsuba.readthedocs.io.

.. warning::

   This page requires two steps before it renders correctly:

   1. Install ``eradiate-mitsuba`` (or build it from source with
      ``pixi run build-kernel``).
   2. Build the API subsite to generate the intermediate RST::

         sphinx-build -b html docs_api _build/html_api

   Then rebuild this site normally.

.. include:: ../../generated/eradiate_api.rst
