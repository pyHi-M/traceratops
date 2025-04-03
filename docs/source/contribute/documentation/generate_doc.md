# Generate documentation from code

## From `ArgumentParser`

*For more details, see the [official documentation of sphinx-argparse](https://sphinx-argparse.readthedocs.io/en/latest/index.html).*

- Activate your dev env and install:
    ```bash
    pip install sphinx-argparse
    ```

- Enable the extension in your sphinx config:
    ```python
    extensions = [
        ...,
        'sphinxarg.ext',
    ]
    ```


- You need a python script with some arguments used for CLI like:

    `<code_source_folder>/example_script.py`

    ```python
    import argparse

    def parse_arguments():
        parser = argparse.ArgumentParser(
            add_help=True,
            description="""You can use three double quote
            - To write a description
            of this script
            - With different lines""",
        )
        # add_argument_group() is optional
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

- Create a .rst (or .md) file inside your source folder:

    ```touch docs/source/example_script.rst```

- Add title and keywords to generate automatically the doc:

    ```rst
    example\_script
    ====================

    .. argparse::
    :ref: <code_source_folder>.example_script.parse_arguments
    :prog: example_script
    ```

- Link this file inside your `index.rst`:

   ```rst

   .. toctree::
      :maxdepth: 1
      example_script
   ```

## From """docstring"""

```{warning}
Not currently used for this project
```

*Generate automatically documentation from your module/class/function docstring.*

Activate your dev env and install:
```bash
pip install sphinx-apidoc
```

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
```bash
sphinx-apidoc -a -o docs/source/ <code_source_folder>/
```
