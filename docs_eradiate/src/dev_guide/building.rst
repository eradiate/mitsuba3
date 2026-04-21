Building from source
====================

.. todo::

   Document the full build process, including:

   - Prerequisites (CMake, ninja, Python ≥ 3.9, pixi)
   - Configuring the kernel: ``pixi run kernel-configure``
   - Building the kernel: ``pixi run build-kernel``
   - Building a Python wheel: ``pixi run wheel-build``
   - Running the test suite: ``pixi run test``
   - Building this documentation: ``sphinx-build -b html docs_eradiate _build/html``
   - Eradiate-specific CMake options

Minimal quickstart:

.. code-block:: bash

   # From the Eradiate monorepo root:
   pixi run build-kernel
   pixi run test-quick
