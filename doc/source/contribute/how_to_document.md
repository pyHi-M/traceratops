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

### 1.2. Generate the doc folder structure

```shell
mkdir doc
cd doc
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

- ```cd doc/```
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
```sphinx-apidoc -a -o doc/source/ <code_source_folder>/```

### 2.2. Generate doc from `ArgumentParser` descriptions

*For more details, see the [official documentation of sphinx-argparse](https://sphinx-argparse.readthedocs.io/en/latest/index.html).*

0. Activate your dev env and `pip install sphinx-argparse`

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
```touch doc/source/example_script.rst```

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

rm -r doc/build/html
sphinx-build -b html doc/source doc/build/html
```

- If you have bug with image display, may be you can refresh cache on your web navigator.

- Custom the design with `html_theme` variable inside `conf.py`
  - For ReadTheDocs: 
    - `pip install sphinx-rtd-theme` 
    - `html_theme = "sphinx_rtd_theme"`

