# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join("..", "..")))

from traceratops._version import __version__  # noqa: E402

autodoc_mock_imports = [
    "numpy",
    "matplotlib",
    "pandas",
    "astropy",
    "tqdm",
    "scipy",
    "sklearn",
]


project = "traceratops"
copyright = "2025, Marcelo Nollmann, Xavier Devos"
author = "Marcelo Nollmann, Xavier Devos"
release = __version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",  # include documentation from docstring
    "sphinx.ext.napoleon",  # allow google or numpy docstring
    "myst_parser",  # parse markdown files to be understood by sphinx
    "sphinxarg.ext",
]

templates_path = ["_templates"]
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output


html_theme = "sphinx_rtd_theme"
html_context = {"default_mode": "light"}
html_static_path = ["_static"]
