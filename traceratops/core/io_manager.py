# -*- coding: utf-8 -*-
"""
Input/Output data manager module
"""

import json
import os

from astropy.table import Table


def create_folder(folder_path: str):
    """Create folder with `makedirs` from os module.
    It's a recursive directory creation function.

    Parameters
    ----------
    folder_path : str
        Path name of folder
    """

    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        print(f"$ Folder '{folder_path}' created successfully.")
    else:
        print(f"! [INFO] Folder '{folder_path}' already exists.")


def load_barcode_dict(file_name):
    """Loads a barcode type dict JSON

    Parameters
    ----------
    file_name : str
        JSON file name

    Returns
    -------
    dict
        Barcode dictionary
    """
    bc_dict = {}
    # Check if the file exists
    if not os.path.exists(file_name):
        print("File does not exist")
    else:
        # Opening JSON file
        with open(file_name, encoding="utf-8") as json_f:
            # returns JSON object as a dictionary
            barcode_type = json.load(json_f)
            print("$ {} barcode dictionary loaded")
            bc_dict = barcode_type

    return bc_dict


def save_json(data, file_name):
    """Save a python dict as a JSON file

    Parameters
    ----------
    data : dict
        Data to save
    file_name : str
        Output JSON file name
    """
    with open(file_name, mode="w", encoding="utf-8") as json_f:
        json.dump(data, json_f, ensure_ascii=False, sort_keys=True, indent=4)


def read_table_from_ecsv(path):
    """Read an astropy Table saved as an ``ecsv`` file."""
    # read ecsv file
    table = Table.read(path, format="ascii.ecsv")

    return table


def save_table_to_ecsv(data, path):
    """
    Save an astropy table into an ``ecsv`` file.
    The input should be an astropy Table object
    """
    data.write(
        path,
        format="ascii.ecsv",
        overwrite=True,
    )
