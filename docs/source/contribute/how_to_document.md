# [How to] Document

*This current documentation was created and configured using [this tutorial](doc_from_scratch.md). It's a set of text files in Markdown (MD) or reStructuredText (RST or ReST) language.*

*We use Sphinx tool to generate HTML pages and ReadTheDocs to host documentation online. You can build localy HTML pages with the bash script `build_doc.sh`. On ReadTheDocs, the documentation is automaticaly build and update when your push code on `main` or `dev` branches to GitHub repository.*

## Document a feature from scratch

[Please, follow this tutorial.](generate_doc.md#from-argumentparser)

## Document a feature enhancement

### Inside your python script

- Update the `description` argument inside the `argparse.ArgumentParser(...)`
- If you create a new CLI option with `.add_argument(...)`, add a `help` argument like this:
```python
parser.add_argument("--z_min", help="Z minimum for a localization.", default=0, type=float)
```

### Update the reliability status

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

### [Optional] Inside the doc file of this script

- Find the doc file inside `docs/source/scripts/` folder
- Add general information or images after this block:
````rst
```{eval-rst}
.. argparse::
:ref: traceratops.figure_him_matrix.parse_arguments
:prog: figure_him_matrix
```
````


## Encapsulate ReST directive in Markdown with MyST

[*Source*](https://deusyss.developpez.com/tutoriels/Python/SphinxDoc/#LVI)

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

*To use it, activate your dev env and `pip install myst-parser`.*

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
