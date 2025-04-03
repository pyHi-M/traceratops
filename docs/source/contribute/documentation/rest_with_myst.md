# Encapsulate ReST directive in Markdown with MyST

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
