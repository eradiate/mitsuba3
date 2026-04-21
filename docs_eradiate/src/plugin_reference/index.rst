.. _sec-eradiate-plugin-reference:

Plugin Reference
================

This section documents all Mitsuba plugins provided by the Eradiate fork.
They are organized into 9 categories corresponding to the subdirectories of
``src/eradiate_plugins/``.

Plugin parameters are documented in tables with four columns: **Parameter**,
**Type**, **Description**, and **Flags**. The Flags column uses the following
abbreviations:

- |exposed| — the parameter is exposed as a scene parameter (traversable)
- |readonly| — exposed but not modifiable after construction
- |differentiable| — the parameter supports automatic differentiation
- |discontinuous| — differentiation requires special handling to avoid bias

.. toctree::
   :maxdepth: 1
   :glob:

   ../../generated/plugins_*
