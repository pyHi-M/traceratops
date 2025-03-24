#!/bin/bash

pip install -e ".[dev]"
rm -r docs/build/html
sphinx-build -b html docs/source docs/build/html
firefox docs/build/html/index.html
