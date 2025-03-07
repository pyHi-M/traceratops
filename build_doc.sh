#!/bin/bash

rm -r docs/build/html
sphinx-build -b html docs/source docs/build/html
