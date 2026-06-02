API documentation pipeline
==========================

The API reference is built in two stages before the main site renders it.

Overview
--------

.. code-block:: text

   docs_eradiate/_docs_api/conf.py  (api_doc_structure)
           │
           │  builder-inited: generate_list_api_callback
           ▼
   docs_eradiate/_docs_api/list_api.rst   ← generated; do not edit by hand
           │
           │  sphinx-build (docs-api task)
           │  autodoc + nanobind patches
           ▼
   docs_eradiate/generated/extracted_rst_api.rst   ← raw RST extracted from docstrings
   docs_eradiate/generated/eradiate_mitsuba_api.rst        ← index with ``.. include::`` ranges
           │
           │  sphinx-build (docs task)
           ▼
   docs_eradiate/_build/html/api_reference/index.html

Stage 1 — API subsite (``docs-api``)
------------------------------------

Source directory: ``docs_eradiate/_docs_api/``

The ``conf.py`` for this subsite does several things during the build:

1. **Loads a Mitsuba variant** (``scalar_rgb``) so that nanobind-bound
   classes are importable by autodoc.

2. **Generates** ``list_api.rst`` at ``builder-inited`` via
   ``generate_list_api_callback``.  The function walks the top-level
   ``mitsuba`` namespace and emits an ``.. autoclass::``,
   ``.. autofunction::``, or ``.. autodata::`` directive for every symbol
   whose fully-qualified name matches a pattern in ``api_doc_structure``.
   The file is sorted alphabetically and must not be edited by hand.

3. **Patches Sphinx's inspect utilities** (``mpatch_ismethod``,
   ``mpatch_isclassmethod``) so that nanobind method objects are classified
   correctly and not skipped.

4. **Intercepts autodoc callbacks** to reformat C++ / nanobind docstrings
   into valid RST:

   - ``process_signature_callback`` — strips and sanitizes C++ type
     annotations from signatures.
   - ``process_docstring_callback`` — converts Markdown code fences to RST
     ``.. code-block::`` directives, adds typed parameter / return entries,
     handles overloaded signatures, and emits each symbol as a
     ``.. py:class::`` / ``.. py:method::`` / … block into the in-memory
     ``extracted_rst`` list.

5. **Writes two generated files** in ``build-finished`` via
   ``write_rst_file_callback``:

   - ``generated/extracted_rst_api.rst`` — the full concatenated RST for all
     documented symbols, one block per symbol.
   - ``generated/eradiate_mitsuba_api.rst`` — an index file that uses
     ``.. include:: /generated/extracted_rst_api.rst`` with ``:start-line:``
     / ``:end-line:`` to splice each symbol's block under the section
     heading defined in ``api_doc_structure``.

.. note::

   On incremental builds Sphinx reuses its pickled environment and skips
   re-reading sources, so the autodoc callbacks do not fire. The
   ``build-finished`` handler detects this (``last_block_name is None``) and
   returns early without overwriting the generated files.  ``list_api.rst``
   is still regenerated (step 2 runs at ``builder-inited``, which always
   fires), but the HTML output and the ``generated/`` files are not updated.

Controlling what is documented
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Adding a symbol** — add a regex pattern to the appropriate section in
``api_doc_structure`` in ``docs_eradiate/_docs_api/conf.py``.
``list_api.rst`` is regenerated automatically on the next build; there is
no need to touch it.

**Organizing into sections** — each key of ``api_doc_structure`` becomes a
section heading in the rendered reference; each value is a list of regex
patterns matched with ``re.fullmatch`` against the fully-qualified Python
name.  Symbols that do not match any pattern are placed in an "Other"
section automatically.

.. code-block:: python

   api_doc_structure = {
       "Constants": [
           r"mitsuba.ERD_MI_([\w]+)",
       ],
       "Medium": [
           r"mitsuba\.Medium",
       ],
       "Extremum": [
           r"mitsuba\.ExtremumSegment([\w]*)",
           r"mitsuba\.ExtremumStructure([\w]*)",
       ],
   }

Stage 2 — main site (``docs``)
------------------------------

Source directory: ``docs_eradiate/``

``api_reference/index.rst`` includes the generated index:

.. code-block:: rst

   .. include:: ../generated/eradiate_mitsuba_api.rst

Sphinx then processes ``eradiate_mitsuba_api.rst``, which splices the relevant line
ranges from ``extracted_rst_api.rst`` into the rendered page. Both generated
files must exist before this build runs, so ``docs-api`` must be run first.

Build tasks
-----------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Task
     - Description
   * - ``pixi run docs-api``
     - Build only the API subsite. Run this after changing
       ``api_doc_structure``.
   * - ``pixi run docs``
     - Build the main site. Run ``pixi run docs-api`` first.
   * - ``pixi run docs-serve``
     - Live-reload server. Watches ``docs_eradiate/`` but ignores
       ``generated/`` to avoid rebuild loops caused by the file-generation
       step.
   * - ``pixi run docs-clean``
     - Remove ``docs_eradiate/_build/``. Run before ``docs`` to force a
       full rebuild.
