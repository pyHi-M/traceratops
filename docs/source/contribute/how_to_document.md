# How to document

## 0. Introduction to Sphinx

*Sphinx is a documentation generator. It is a tool for translating text comments, in reST ([reStructuredText](https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html)) by default or Md ([markdown](https://www.markdownguide.org/cheat-sheet/)) via [MyST extension](https://myst-parser.readthedocs.io/en/latest/), in source files into HTML or PDF files using LaTeX.*

### 0.1. Features

- Support for several output formats (HTML, PDF, EPUB, manual pages, etc.)
- Automatic generation of code-related documentation (with support for many languages)
- Easy reference between pages
- Its management of extensions to adapt it to all situations and languages

### 0.2. Comparison of reST and MyST 

Source : https://deusyss.developpez.com/tutoriels/Python/SphinxDoc/#LVI

ReST:

```reStructuredText
.. toctree::
caption: Episodes
maxdepth: 1

basics
creating-using-web
creating-using-desktop
contributing
doi
websites
```

MyST:

To use it, activate your dev env and `pip install myst-parser`

````markdown
```{toctree}
---
caption: Episodes
maxdepth: 1
---

basics
creating-using-web
creating-using-desktop
contributing
doi
websites
```
````

## 1. First steps

***To be done only once.***

### 1.1. Create your dev environment

```shell
conda create -n traceratops_dev python=3.11
conda activate traceratops_dev
pip install sphinx
```

### 1.2. Generate the "docs" folder structure

```shell
mkdir docs
cd docs
sphinx-quickstart
```

The last command ask you few question:

```
> Separate source and build directories: YES
> Project name:
> Author name(s):
> Project release []: 0.0.1
```

It create few files:

```shell
.
├── build (empty folder, used when you generate doc)
├── make.bat (don't modify, it's used to build doc)
├── Makefile (don't modify, it's used to build doc)
└── source (folder that you will fill with doc pages) 
    ├── conf.py (the place to configurate your doc)
    ├── index.rst (main doc page)
    ├── _static	(put your static data here, like picture)
    └── _templates (unused?)
```

### 1.3. Build locally your doc

- ```cd docs/```
- ```sphinx-build -b html source/ build/html```

### 1.4. Open your doc inside a web navigator

```firefox build/html/index.html```

## 2. generate documentation from code

### 2.1. Generate automatically doc from your code docstring

Activate your dev env and `pip install sphinx-apidoc`

Open file: `conf.py`
add:

```python
import os
import sys
sys.path.insert(0, os.path.abspath("../../<code_source_folder>/"))
extensions = [
    "sphinx.ext.autodoc",  # include documentation from docstring
    "sphinx.ext.napoleon",  # allow google or numpy docstring
    "myst_parser",  # parse markdown files to be understood by sphinx
]
```

Run inside a terminal:
```sphinx-apidoc -a -o docs/source/ <code_source_folder>/```

### 2.2. Generate doc from `ArgumentParser` descriptions

*For more details, see the [official documentation of sphinx-argparse](https://sphinx-argparse.readthedocs.io/en/latest/index.html).*

0. Activate your dev env and `pip install sphinx-argparse`

Enable the extension in your sphinx config:
```python
extensions = [
    ...,
    'sphinxarg.ext',
]
```


1. You need a python script with some arguments used for CLI.
`<code_source_folder>/example_script.py`:

```python
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(
        add_help=True,
        description="""You can use three double quote
        - To wirte a description 
          of this script
        - With differents lines""",
    )
    # add_argument_group() is optionnal
    # you can directly use parser.add_argument()
    parser_required = parser.add_argument_group("Required arguments")
    parser_required.add_argument("--input", help="Input data path")
    parser_opt = parser.add_argument_group("Optional arguments")
    parser_opt.add_argument("--output", help="Output filename")
    return parser # need to return the parser object
    
# How to collect argument value
if __name__ == "__main__":
    parser = parse_arguments()
    args = parser.parse_args()
    input_path = args.input
    output_path = args.output
```

2. Create a .rst (or .md) file inside your source folder:
```touch docs/source/example_script.rst```

3. Add title and keywords to generate automatically the doc:

```rst
example\_script
====================

.. argparse::
   :ref: <code_source_folder>.example_script.parse_arguments
   :prog: example_script
```

4. Link this file inside your `index.rst`:

   ```rst
   
   .. toctree::
      :maxdepth: 1
      example_script
   ```

5. Build locally your doc:

```sphinx-build -b html source/ build/html```

## 3. Tips

- Think to delete content of build folder, this can solve strange bug...

  You can use a bash script to build:

```bash
#!/bin/bash

rm -r docs/build/html
sphinx-build -b html docs/source docs/build/html
```

- If you have bug with image display, may be you can refresh cache on your web navigator.

- Custom the design with `html_theme` variable inside `conf.py`
  - For ReadTheDocs theme: 
    - `pip install sphinx-rtd-theme` 
    - `html_theme = "sphinx_rtd_theme"`

- On ReadTheDocs you can customise the branch used to build the doc:
  - readthedocs.org > dashboard > "your_project" > Settings > Default branch

## 4. Deploy on ReadTheDocs

1. Create an account on https://readthedocs.org/

2. On ReadTheDcos, connect with github:
  - Settings > Account > Connected services > Add new connection

3. On GitHub (Your personnal account) give access for readthedocs:
  - Profil > Settings > (Integrations) > Applications > Authorized OAuth Apps > Read the Docs Community (readthedocs.org)
  - You should see a list intitle "Organization access", clic on "Grant" button of the organization that host your repository.

4. On GitHub (code repository) configure your project for readthedocs:
  - On the default branch (`main` or `master` for old repo)
  - add a `.readthedocs.yaml` file at the root folder
  ```yaml
  # .readthedocs.yaml
  # Read the Docs configuration file
  # See https://docs.readthedocs.io/en/stable/config-file/v2.html for details

  # Required
  version: 2

  # Set the version of Python and other tools you might need
  build:
    os: ubuntu-24.04
    tools:
      python: "3.11"

  # Build documentation in the docs/ directory with Sphinx
  sphinx:
    builder: html
    configuration: docs/source/conf.py

  # Optionally declare the Python requirements required to build your docs
  python:
    install:
    - requirements: docs/requirements.txt
    - method: pip
      path: .
  ```

5. Add inside `docs` folder, a `requirements.txt` file with only the dependancies used to build doc:
  ```text
  sphinx
  numpydoc
  sphinx-rtd-theme
  myst-parser
  sphinx-argparse
  ```

6. inside the doc `conf.py` file, mock the code dependancies:

  *If you generate documentation via `sphinx-apidoc` or `sphinx-argparse`, the code will be compiled in order to extract the useful information for the doc in the code. The mock system, thanks to the `autodoc_mock_imports` variable, allows you to ignore external import instructions in the code when building the documentation.*

  ```python
  #...
  import os
  import sys

  autodoc_mock_imports = [
      "numpy",
      "matplotlib",
      "pandas",
      "astropy",
      "tqdm",
  ]

  sys.path.insert(0, os.path.abspath("../../traceratops/"))
  #...
  ```

  7. Go back to your [ReadTheDocs dashboard](https://app.readthedocs.org/dashboard/):

  - "Add project"
  - Follow instruction, use "Configure manually" if your repo isn't find automatically.
