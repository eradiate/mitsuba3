Introduction
============

.. todo::

   Write a user-facing introduction covering:

   - Installation (``pip install eradiate-mitsuba`` or from source)
   - Variant selection with ``mi.set_variant(...)``
   - Basic scene setup: loading plugins, setting parameters
   - Running a simple simulation and reading the output
   - Differences from upstream Mitsuba 3

A minimal starting point:

.. code-block:: python

   import mitsuba as mi

   mi.set_variant('scalar_mono')

   scene = mi.load_dict({
       'type': 'scene',
       # ... scene description
   })
