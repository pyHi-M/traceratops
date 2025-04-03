# Documentation

*This current documentation was created and configured using [this tutorial](documentation/doc_from_scratch.md). It's a set of text files in Markdown (MD) or reStructuredText (RST or ReST) language.*

*We use Sphinx tool to generate HTML pages and ReadTheDocs to host documentation online. On ReadTheDocs, the documentation is automatically build and update when your push code on `main` or `dev` branches to GitHub repository.*

## Shortcut

```{eval-rst}

.. panels::
    :card: + intro-card text-center
    :column: col-lg-4 col-md-6 col-sm-6 col-xs-12 p-2

    ---

    .. link-button:: documentation/doc_from_scratch
            :type: ref
            :text: Create a doc
            :classes: btn-block btn-info stretched-link text-white

    ---

    .. link-button:: documentation/generate_doc
            :type: ref
            :text: Doc from code
            :classes: btn-block btn-info stretched-link text-white

    ---

    .. link-button:: documentation/rest_with_myst
            :type: ref
            :text: ReST with MyST
            :classes: btn-block btn-info stretched-link text-white
```


```{toctree}
:maxdepth: 1
:hidden:

documentation/doc_from_scratch
documentation/generate_doc
documentation/rest_with_myst
```


## Build locally


```{eval-rst}

.. panels::
    :card: + intro-card

    ---
    :column: col-lg-4 col-md-12 p-2

    **Isolated virtual env**
    (similar with RTD env):

    .. code-block:: bash

       bash build_docs.sh

    ---
    :column: col-lg-8 col-md-12 p-2

    **Development env** (faster):

    .. code-block:: bash

       rm -r docs/build/html
       sphinx-build -b html docs/source docs/build/html

```

## New script

[Please, follow this tutorial.](documentation/generate_doc.md#from-argumentparser)

## Feature enhancement

### Python script

*Document inside your python script.*

- Update the `description` argument inside the `argparse.ArgumentParser(...)`
- If you create a new CLI option with `.add_argument(...)`, add a `help` argument like this:
```python
parser.add_argument("--z_min", help="Z minimum for a localization.", default=0, type=float)
```

### Reliability status

Choose the status of your main feature (can be non-exhaustive):
- `development` = no test available
- `beta-test` = only manual test(s)
- `stable` = automatic test(s)

Find the doc file inside `docs/source/scripts/` folder.

Add just after the main title: "**Reliability status:** (in bolt) `status value` (in backtick)"
```markdown
# figure_him_matrix

**Reliability status**: `stable`
```

### Markdown script (Optional)

*Document inside the markdown file linked to this script.*

- Find the doc file inside `docs/source/scripts/` folder
- Add general information or images after this block:
````rst
```{eval-rst}
.. argparse::
:ref: traceratops.figure_him_matrix.parse_arguments
:prog: figure_him_matrix
```
````
