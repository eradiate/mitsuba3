#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Eradiate Mitsuba fork — Sphinx documentation configuration

import os
import re
import sys

sys.path.insert(0, os.path.abspath("_exts"))

# -- Version info sourced from mitsuba.h ----------------------------------


def _read_version():
    header = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..",
        "include",
        "mitsuba",
        "mitsuba.h",
    )
    text = open(header, encoding="utf-8").read()
    major, minor, patch = re.findall(
        r"#define ERD_MI_VERSION_(?:MAJOR|MINOR|PATCH)\s+(\d+)", text
    )[:3]
    return f"{major}.{minor}", f"{major}.{minor}.{patch}"


# -- General configuration ------------------------------------------------

source_suffix = [".rst", ".md"]
master_doc = "index"

project = "Eradiate Mitsuba"
copyright = "2024, Rayference"
author = "Rayference"

version, release = _read_version()

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.todo",
    "sphinx_design",
    "sphinx_copybutton",
    "sphinx_tabs.tabs",
    "sphinxcontrib.bibtex",
    "nbsphinx",
    "myst_parser",
    "pluginparameters",
    "subfig",
]

exclude_patterns = [
    "_build",
    "**.ipynb_checkpoints",
    "_docs_api/**",
    "plugin_reference/section_*.rst",
    "generated/extracted_rst_api.rst",
    "generated/eradiate_api.rst",
]

default_role = "any"
language = "en"
primary_domain = "py"
highlight_language = "python"

rst_prolog = r"""
.. role:: paramtype

.. role:: monosp

.. |spectrum| replace:: :paramtype:`spectrum`
.. |texture| replace:: :paramtype:`texture`
.. |float| replace:: :paramtype:`float`
.. |bool| replace:: :paramtype:`boolean`
.. |int| replace:: :paramtype:`integer`
.. |false| replace:: :monosp:`false`
.. |true| replace:: :monosp:`true`
.. |string| replace:: :paramtype:`string`
.. |bsdf| replace:: :paramtype:`bsdf`
.. |phase| replace:: :paramtype:`phase`
.. |point| replace:: :paramtype:`point`
.. |vector| replace:: :paramtype:`vector`
.. |transform| replace:: :paramtype:`transform`
.. |volume| replace:: :paramtype:`volume`
.. |tensor| replace:: :paramtype:`tensor`

.. |drjit| replace:: :monosp:`drjit`
.. |numpy| replace:: :monosp:`numpy`

.. |nbsp| unicode:: 0xA0
   :trim:

.. |exposed| replace:: :abbr:`P (This parameter will be exposed as a scene parameter)`
.. |readonly| replace:: :abbr:`R (This parameter will be exposed as a scene parameter, but cannot be modified.)`
.. |differentiable| replace:: :abbr:`∂ (This parameter is differentiable)`
.. |discontinuous| replace:: :abbr:`D (This parameter might introduce discontinuities. Therefore it requires special handling during differentiation to prevent bias))`

"""

# -- Theme ----------------------------------------------------------------

html_theme = "shibuya"
html_title = "Eradiate Mitsuba"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_show_sourcelink = False

html_theme_options = {
    "github_url": "https://github.com/eradiate/eradiate-mitsuba",
    "accent_color": "teal",
}

# -- nbsphinx -------------------------------------------------------------

nbsphinx_execute = "never"

# -- myst-parser ----------------------------------------------------------

myst_enable_extensions = ["colon_fence"]

# -- intersphinx ----------------------------------------------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "mitsuba": ("https://mitsuba.readthedocs.io/en/stable/", None),
}

# -- bibtex ---------------------------------------------------------------

bibtex_bibfiles = ["references.bib"]

# -- todo -----------------------------------------------------------------

todo_include_todos = True

# -- Plugin doc generation ------------------------------------------------

build_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "generated")


def custom_step(app):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import generate_plugin_doc

    if not os.path.exists(build_dir):
        os.mkdir(build_dir)
    generate_plugin_doc.generate(build_dir)
    # Create a placeholder eradiate_api.rst if the docs_api subsite hasn't been
    # built yet, so the include in src/api_reference/index.rst doesn't fail.
    api_rst = os.path.join(build_dir, "eradiate_api.rst")
    if not os.path.exists(api_rst):
        with open(api_rst, "w", encoding="utf-8") as f:
            f.write(
                ".. note::\n\n"
                "   Run ``sphinx-build -b html _docs_api _build/html_api`` first\n"
                "   to generate the API reference content.\n"
            )


def setup(app):
    app.connect("builder-inited", custom_step)
