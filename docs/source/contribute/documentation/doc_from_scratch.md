# Setup documentation (from scratch)

```{note}
**If everything is already set up**, simply run this script to build the documentation locally and verify your changes before committing them.

`bash build_doc.sh`
```


***To be done only once.***

## 1. First steps

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
    ├── conf.py (the place to configure your doc)
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

- [Generate automatically doc from your code docstring](generate_doc.md#from-docstring)

- [Generate doc from `ArgumentParser` descriptions](generate_doc.md#from-argumentparser)

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

3. On GitHub (Your personal account) give access for readthedocs:
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

5. Add inside `docs` folder, a `requirements.txt` file with only the dependencies used to build doc:
  ```text
  sphinx
  numpydoc
  sphinx-rtd-theme
  myst-parser
  sphinx-argparse
  ```

6. inside the doc `conf.py` file, mock the code dependencies:

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
