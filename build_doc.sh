#!/bin/bash

# >>> Conda init <<<
source ~/miniconda3/etc/profile.d/conda.sh
conda create -n test_build_doc python=3.11 --yes
conda activate test_build_doc
pip install --no-deps .
pip install -r docs/requirements.txt
rm -r docs/build/html
sphinx-build -b html docs/source docs/build/html
firefox docs/build/html/index.html
conda deactivate
conda env remove --name test_build_doc
