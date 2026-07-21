"""
New release versions are made using [Semantic Versioning](https://semver.org/).

> Given a version number MAJOR.MINOR.PATCH, increment the:

    1. MAJOR version when you make incompatible API changes,
    2. MINOR version when you add functionality in a backwards compatible manner, and
    3. PATCH version when you make backwards compatible bug fixes.

    Additional labels for pre-release and build metadata are available as extensions to the MAJOR.MINOR.PATCH format.

Examples:
"0.0.0.dev0"
"1.2.34.dev0"
"1.2.34a0"
"1.2.34"
"""

import importlib.metadata

__version__ = importlib.metadata.version("traceratops")
