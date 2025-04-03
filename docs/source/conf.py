# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import shutil

from traceratops._version import __version__


def inject_pr_template():
    def __add_title(file_path):
        with open(file_path, "r", encoding="utf-8") as f:
            content = f.read()
        new_content = "# Pull Request template (COPY)\n\n" + content
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(new_content)

    def __convert_checkboxes(file_path):
        updated_lines = []
        with open(file_path, "r", encoding="utf-8") as f:
            for line in f:
                stripped = line.strip()
                if stripped.startswith("- [ ]"):
                    label = stripped[5:].strip()  # le texte après '- [ ]'
                    html = f'<label><input type="checkbox"> {label}</label>'
                    updated_lines.append(html + "\n")
                else:
                    updated_lines.append(line)

        with open(file_path, "w", encoding="utf-8") as f:
            f.writelines(updated_lines)

    source_md = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "PULL_REQUEST_TEMPLATE.md")
    )
    destination_md = os.path.join(
        os.path.dirname(__file__),
        "contribute",
        "checklist",
        "pr_template_COPY.md",
    )
    shutil.copyfile(source_md, destination_md)

    __add_title(destination_md)
    __convert_checkboxes(destination_md)


inject_pr_template()

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
    "sphinx_panels",  # for creating panels like pandas or numpy main doc page
]

templates_path = ["_templates"]
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output


html_theme = "sphinx_rtd_theme"
html_context = {"default_mode": "light"}
html_static_path = ["_static"]

myst_heading_anchors = 2
