#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script will load a localizations table and analyze it.

It will specifically produce a plot with:
- the number of localizations per barcode
- the snr distribution per barcode
- scatterplot of the snr versus z
- scatterplot of roundness versus skew

"""

import argparse

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from traceratops.core.chromatin_trace_table import ChromatinTraceTable
from traceratops.core.localization_table import LocalizationTable

font = {"weight": "normal", "size": 12}
matplotlib.rc("font", **font)


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-L", "--localization", help="Localizations file path")
    return parser


def main():
    parser = parse_arguments()
    args = parser.parse_args()
    # [loops over lists of datafolders]
    process_intensities(args.localization, args.trace)
    print("Finished execution")


if __name__ == "__main__":
    main()
