#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
contains functions and classes needed for the analysis and plotting of HiM matrices
"""


import csv
import itertools
import json
import os

import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as npl
from sklearn import manifold
from sklearn.model_selection import GridSearchCV, LeaveOneOut
from sklearn.neighbors import KernelDensity
from tqdm import trange


class AnalysisHiMMatrix:
    """
    Main use is to produce paper quality figures of HiM matrices, 3-way interaction matrices and HiM matrix ratios
    """

    def __init__(self, run_parameters, root_folder="."):
        self.data_folder = root_folder + os.sep + "scHiMmatrices"
        self.run_parameters = run_parameters
        self.root_folder = root_folder
        self.data = dict()
        self.data_files = dict()
        self.folders_to_load = []
        self.number_barcodes = 0

    def load_data(self):
        """
        loads dataset

        Returns
        -------
        self.foldes2Load contains the parameters used for the processing of HiM matrices.
        self.data_files dictionary containing the extensions needed to load data files
        self.data dictionary containing the datasets loaded
        """

        # loads datasets: parameter files
        filename_list_data_json = (
            self.root_folder + os.sep + self.run_parameters["parametersFileName"]
        )
        with open(filename_list_data_json, encoding="utf-8") as json_file:
            list_data = json.load(json_file)

        dataset_name = list(list_data.keys())[0]
        print(f"Dataset: {dataset_name}")

        output_filename = (
            self.data_folder
            + os.sep
            + dataset_name
            + "_label:"
            + self.run_parameters["label"]
            + "_action:"
            + self.run_parameters["action"]
        )

        filename_parameters_json = f"{output_filename}_parameters.json"
        with open(filename_parameters_json, encoding="utf-8") as json_file:
            folders_to_load = json.load(json_file)
        print(f"Loading parameter file:{filename_parameters_json}")

        # Creates filenames to be loaded
        data_files = {
            "ensembleContactProbability": "_ensembleContactProbability.npy",
            "SCmatrixCollated": "_SCmatrixCollated.npy",
            "SClabeledCollated": "_SClabeledCollated.npy",
        }

        if "3wayContacts_anchors" in list_data[dataset_name]:
            for i_anchor in list_data[dataset_name]["3wayContacts_anchors"]:
                new_key = f"anchor:{str(i_anchor - 1)}"
                data_files[new_key] = f"_{new_key}_ensemble3wayContacts.npy"
        else:
            print("No anchors found")

        # loads datasets: numpy matrices
        data = {}
        print(f"Loading datasets from: {output_filename}")
        for i_data_file, value in data_files.items():
            print(
                f"Loaded: {i_data_file}: <{os.path.basename(output_filename + value)}>"
            )
            data[i_data_file] = np.load(
                output_filename + data_files[i_data_file]
            ).squeeze()

        # loads datasets: lists
        run_name = load_list(f"{output_filename}_runName.csv")
        data["runName"] = run_name
        print(f"""Loaded runNames: {data["runName"]}""")

        data["uniqueBarcodes"] = load_list(f"{output_filename}_uniqueBarcodes.csv")
        print(f"""Loaded barcodes #: {data["uniqueBarcodes"]}""")
        self.number_barcodes = len(data["uniqueBarcodes"])

        print(f"""Total number of cells loaded: {data["SCmatrixCollated"].shape[2]}""")
        print(f"""Number Datasets loaded: {len(data["runName"])}""")

        # Exports data
        self.data = data
        self.data_files = data_files
        self.folders_to_load = folders_to_load
        self.list_data = list_data
        self.dataset_name = dataset_name

    # functions

    def plot_2d_matrix_simple(
        self,
        ifigure,
        matrix,
        unique_barcodes,
        yticks,
        xticks,
        cmtitle="probability",
        c_min=0,
        c_max=1,
        c_m="coolwarm",
        fontsize=12,
        colorbar=False,
        axis_ticks=False,
        n_cells=0,
        n_datasets=0,
        show_title=False,
        fig_title="",
    ):
        pos = ifigure.imshow(matrix, cmap=c_m)  # colormaps RdBu seismic

        if show_title:
            title_text = f"{fig_title} | N = {n_cells} | n = {n_datasets}"
            ifigure.title.set_text(title_text)

        # plots figure
        if xticks:
            ifigure.set_xlabel("barcode #", fontsize=fontsize)
            if not axis_ticks:
                ifigure.set_xticklabels(())
            else:
                print(f"barcodes:{unique_barcodes}")
                # ifigure.set_xticks(np.arange(matrix.shape[0]),unique_barcodes)
                ifigure.set_xticklabels(unique_barcodes)

        else:
            ifigure.set_xticklabels(())
        if yticks:
            ifigure.set_ylabel("barcode #", fontsize=fontsize)
            if not axis_ticks:
                ifigure.set_yticklabels(())
            else:
                # ifigure.set_yticks(np.arange(matrix.shape[0]), unique_barcodes)
                ifigure.set_yticklabels(unique_barcodes)
        else:
            ifigure.set_yticklabels(())

        for xtick, ytick in zip(
            ifigure.xaxis.get_majorticklabels(), ifigure.yaxis.get_majorticklabels()
        ):
            xtick.set_fontsize(fontsize)
            ytick.set_fontsize(fontsize)

        if colorbar:
            cbar = plt.colorbar(pos, ax=ifigure, fraction=0.046, pad=0.04)
            cbar.minorticks_on()
            cbar.set_label(cmtitle, fontsize=float(fontsize) * 0.85)
            pos.set_clim(vmin=c_min, vmax=c_max)

        pos.set_clim(vmin=c_min, vmax=c_max)

        return pos

    def n_cells_loaded(self):
        if self.run_parameters["action"] == "labeled":
            cells_with_label = [
                idx for idx, x in enumerate(self.data["SClabeledCollated"]) if x > 0
            ]
            n_cells = len(cells_with_label)
        elif self.run_parameters["action"] == "unlabeled":
            cells_with_label = [
                idx for idx, x in enumerate(self.data["SClabeledCollated"]) if x == 0
            ]
            n_cells = len(cells_with_label)
        else:
            n_cells = self.data["SCmatrixCollated"].shape[2]
        print(f"n_cells selected with label: {n_cells}")
        return n_cells

    def retrieve_sc_matrix(self):
        """
        retrieves single cells that have the label requested

        Returns
        -------
        self.sc_matrix_selected

        """
        n_cells = self.n_cells_loaded()
        sc_matrix_selected = np.zeros(
            (self.number_barcodes, self.number_barcodes, n_cells)
        )

        if self.run_parameters["action"] == "labeled":
            cells_with_label = [
                idx for idx, x in enumerate(self.data["SClabeledCollated"]) if x > 0
            ]
            for new_cell, i_cell in enumerate(cells_with_label):
                sc_matrix_selected[:, :, new_cell] = self.data["SCmatrixCollated"][
                    :, :, i_cell
                ]
        elif self.run_parameters["action"] == "unlabeled":
            cells_with_label = [
                idx for idx, x in enumerate(self.data["SClabeledCollated"]) if x == 0
            ]
            for new_cell, i_cell in enumerate(cells_with_label):
                sc_matrix_selected[:, :, new_cell] = self.data["SCmatrixCollated"][
                    :, :, i_cell
                ]
        else:
            sc_matrix_selected = self.data["SCmatrixCollated"]
        print(f"n_cells retrieved: {sc_matrix_selected.shape[2]}")
        self.sc_matrix_selected = sc_matrix_selected


def load_list(file_name):
    with open(file_name, newline="", encoding="utf-8") as csvfile:
        spamreader = csv.reader(csvfile, delimiter=" ", quotechar="|")
        run_name = []
        for row in spamreader:
            if len(run_name) > 0:
                run_name.append(row)
            else:
                run_name = row
    return run_name


def list_sc_to_keep(p, mask):
    if p["action"] == "all":
        try:
            cells_to_plot = range(len(mask))
        except TypeError:
            print(mask)
            cells_to_plot = range(mask.shape[0])
    elif p["action"] == "labeled":
        a = [i for i in range(len(mask)) if mask[i] == 1]
        cells_to_plot = a
    elif p["action"] == "unlabeled":
        a = [i for i in range(len(mask)) if mask[i] == 0]
        cells_to_plot = a

    print(
        f'>> label: {p["label"]}\t action:{p["action"]}\
            \t Ncells2plot:{max(cells_to_plot)}\t Ncells in sc_matrix:{len(mask)}'
    )

    return cells_to_plot


def plot_matrix(
    sc_matrix_collated,
    unique_barcodes,
    pixel_size,
    number_rois=1,
    output_filename="test",
    log_name_md="log.md",
    clim=1.4,
    c_m="seismic",
    figtitle="PWD matrix",
    cmtitle="distance, um",
    n_cells=0,
    mode="median",
    inverse_matrix=False,
    c_min=0,
    cells_to_plot=None,
    filename_ending="_HiMmatrix.png",
    font_size=22,
):
    if cells_to_plot is None:
        cells_to_plot = []
    ######################################################
    # Calculates ensemble matrix from single cell matrices
    ######################################################

    if len(sc_matrix_collated.shape) == 3:
        # matrix is 3D and needs combining SC matrices into an ensemble matrix
        if len(cells_to_plot) == 0:
            cells_to_plot = range(sc_matrix_collated.shape[2])

        mean_sc_matrix, keep_plotting = calculate_ensemble_pwd_matrix(
            sc_matrix_collated, pixel_size, cells_to_plot, mode=mode
        )

    else:
        # already an ensemble matrix --> no need for further treatment
        if mode == "counts":
            mean_sc_matrix = sc_matrix_collated
        else:
            mean_sc_matrix = pixel_size * sc_matrix_collated
        keep_plotting = True
    if keep_plotting:
        # no errors occurred

        # Calculates the inverse distance matrix if requested in the argument.
        if inverse_matrix:
            mean_sc_matrix = np.reciprocal(mean_sc_matrix)

        # plots figure
        plt.figure(figsize=(15, 15))
        pos = plt.imshow(mean_sc_matrix, cmap=c_m)  # colormaps RdBu seismic
        plt.xlabel("barcode #", fontsize=float(font_size) * 1.2)
        plt.ylabel("barcode #", fontsize=float(font_size) * 1.2)
        plt.title(
            f"{figtitle} | {str(mean_sc_matrix.shape[0])} barcodes | n={str(n_cells)}",
            fontsize=float(font_size) * 1.3,
        )

        plt.xticks(
            np.arange(sc_matrix_collated.shape[0]), unique_barcodes, fontsize=font_size
        )
        plt.yticks(
            np.arange(sc_matrix_collated.shape[0]), unique_barcodes, fontsize=font_size
        )
        cbar = plt.colorbar(pos, fraction=0.046, pad=0.04)
        cbar.ax.tick_params(labelsize=float(font_size) * 0.8)
        cbar.minorticks_on()
        cbar.set_label(cmtitle, fontsize=float(font_size) * 1.0)
        plt.clim(c_min, clim)

        if len(output_filename.split(".")) > 1:
            if output_filename.split(".")[1] == "png":
                out_fn = output_filename.split(".")[0] + filename_ending
            elif len(output_filename.split(".")[1]) == 3:
                # keeps original extension
                out_fn = output_filename
            else:
                # most likely the full filename contains other '.' in addition to that in the extension
                out_fn = output_filename + filename_ending
        else:
            out_fn = output_filename + filename_ending
        plt.savefig(out_fn)
        if not is_notebook():
            plt.close()
        if "png" not in out_fn:
            out_fn += ".png"
    else:
        # errors during pre-processing
        print("Error plotting figure. Not executing script to avoid crash.")

    return mean_sc_matrix


def shuffle_matrix(matrix, index):
    new_size = len(index)
    new_matrix = np.zeros((new_size, new_size, len(matrix[0][0])))

    if new_size > matrix.shape[0]:
        raise ValueError(
            f"Error: shuffle size {new_size} is larger than matrix dimensions {matrix.shape[0]}\nShuffle: {index}"
        )
    for i, j in itertools.product(range(new_size), range(new_size)):
        if index[i] < matrix.shape[0] and index[j] < matrix.shape[0]:
            new_matrix[i, j] = matrix[index[i], index[j]]
        else:
            raise ValueError(
                f"Out of index; matrix.shape[0]: {matrix.shape[0]} |i: {i} |index[i]: {index[i]} |j: {j} |index[j]: {index[j]}"
            )
    return new_matrix


def decodes_trace(single_trace):
    """
    from a trace entry, provides Numpy array with coordinates, barcode and trace names

    Parameters
    ----------
    single_trace : astropy table
        astropy table for a single trace.

    Returns
    -------
    list of barcodes
    x, y and z coordinates as numpy arrays,
    trace name as string

    """
    barcodes, X, Y, Z = (
        single_trace["Barcode #"],
        single_trace["x"],
        single_trace["y"],
        single_trace["z"],
    )
    trace_name = single_trace["Trace_ID"][0][:3]

    return barcodes, X, Y, Z, trace_name


def write_xyz_2_pdb(file_name, single_trace, barcode_type=dict()):
    # writes xyz coordinates to a PDB file with pseudoatoms
    # file_name : string of output file path, e.g. '/foo/bar/test2.pdb'
    # xyz      : n-by-3 numpy array with atom coordinates

    default_atom_name = "xxx"
    barcodes, X, Y, Z, trace_name = decodes_trace(single_trace)

    # builds NP array
    xyz = np.transpose(np.array([X, Y, Z]))

    # calculates center of mass
    center_of_mass = np.mean(X), np.mean(Y), np.mean(Z)

    # recenters and converts to A
    unit_conversion = 10.0  # converts from nm to Angstroms
    xyz = unit_conversion * (xyz - center_of_mass)

    # writes PDB file
    n_atoms = xyz.shape[0]

    # defines atom names from barcode properties
    if len(barcode_type) < 1:
        # all atoms have the same identity
        print("did not find barcode_type dictionary")
        for i, barcode in enumerate(barcodes):
            barcode_type[str(barcode)] = default_atom_name
    else:
        # adds missing keys
        for barcode in barcodes:
            if str(barcode) not in barcode_type.keys():
                barcode_type[str(barcode)] = default_atom_name
                print(f"$ fixing key {barcode} as not found in dict")

    """
        COLUMNS        DATA TYPE       CONTENTS
    --------------------------------------------------------------------------------
     1 -  6        Record name     "ATOM  "
     7 - 11        Integer         Atom serial number.
    13 - 16        Atom            Atom name.
    17             Character       Alternate location indicator.
    18 - 20        Residue name    Residue name.
    22             Character       Chain identifier.
    23 - 26        Integer         Residue sequence number.
    27             AChar           Code for insertion of residues.
    31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.
    39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.
    47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.
    55 - 60        Real(6.2)       Occupancy.
    61 - 66        Real(6.2)       Temperature factor (Default = 0.0).
    73 - 76        LString(4)      Segment identifier, left-justified.
    77 - 78        LString(2)      Element symbol, right-justified.
    79 - 80        LString(2)      Charge on the atom.
    """

    with open(file_name, mode="w+", encoding="utf-8") as fid:
        # atom coordinates
        # txt = "HETATM  {: 3d}  C{:02d} {} P   1      {: 5.3f}  {: 5.3f}  {: 5.3f}  0.00  0.00      PSDO C  \n"
        # txt = "HETATM  {: 3d}  {} {} P{: 3d}      {: 5.3f}  {: 5.3f}  {: 5.3f}  0.00  0.00      PSDO C  \n"
        # for i in range(n_atoms):
        #    atom_name = barcode_type[str(barcodes[i])]
        # fid.write(txt.format(i + 1, i + 1, trace_name, int(barcodes[i]), xyz[i, 0], xyz[i, 1], xyz[i, 2]))
        #    fid.write(txt.format(i + 1, atom_name, trace_name, int(barcodes[i]), xyz[i, 0], xyz[i, 1], xyz[i, 2]))

        # fills fields with correct spacing
        field_record = "HETATM"
        field_atom_number = " {:4d}"
        field_atom_name = "  {}"
        field_alternative_location_indicator = " "
        field_res_name = trace_name + " "
        field_chain_identifier = " "
        field_residue_seq_number = "{:4d}"
        field_code_insertion = "    "
        field_X = "{}"
        field_Y = "{}"
        field_Z = "{}"  # " {:0<7.3f}"
        field_occupancy = "   0.0"
        field_temp_factor = "   0.0"
        field_segment_identifier = "      " + "PSDO"
        field_element_symbol = " X"
        field_charge_atom = " X"
        txt = (
            field_record
            + field_atom_number
            + field_atom_name
            + field_alternative_location_indicator
            + field_res_name
            + field_chain_identifier
            + field_residue_seq_number
            + field_code_insertion
            + field_X
            + field_Y
            + field_Z
            + field_occupancy
            + field_temp_factor
            + field_segment_identifier
            + field_element_symbol
            + field_charge_atom
            + "\n"
        )

        # txt = "HETATM  {: 3d}  C{:02d} {} P   1      {: 5.3f}  {: 5.3f}  {: 5.3f}  0.00  0.00      PSDO C  \n"
        for i in range(n_atoms):
            atom_name = barcode_type[str(barcodes[i])]
            fid.write(
                txt.format(
                    i + 1,
                    atom_name,
                    int(barcodes[i]),
                    " {:0<7.3f}".format(xyz[i, 0])[:8],
                    " {:0<7.3f}".format(xyz[i, 1])[:8],
                    " {:0<7.3f}".format(xyz[i, 2])[:8],
                )
            )

        # connectivity
        txt1 = "CONNECT  {: 3d}  {: 3d}\n"
        txt2 = "CONNECT  {: 3d}  {: 3d}  {: 3d}\n"

        # first line of connectivity
        fid.write(txt1.format(1, 2))

        # consecutive lines
        for i in range(2, n_atoms):
            fid.write(txt2.format(i, i - 1, i + 1))

        # last line
        fid.write(txt1.format(i + 1, i))

        print(f"Done writing {file_name:s} with {n_atoms:d} atoms.")


def distances_2_coordinates(distances):
    """Infer coordinates from distances"""
    N = distances.shape[0]
    d_0 = []

    # pre-caching
    cache = {}
    for j in range(N):
        sumi = sum(distances[j, k] ** 2 for k in range(j + 1, N))
        cache[j] = sumi

    # compute distances from center of mass
    sum2 = sum(cache[j] for j in range(N))
    for i in range(N):
        sum1 = cache[i] + sum(distances[j, i] ** 2 for j in range(i + 1))

        val = 1 / N * sum1 - 1 / N**2 * sum2
        d_0.append(val)

    # generate gram matrix
    gram = np.zeros(distances.shape)
    for row in range(distances.shape[0]):
        for col in range(distances.shape[1]):
            dists = d_0[row] ** 2 + d_0[col] ** 2 - distances[row, col] ** 2
            gram[row, col] = 1 / 2 * dists

    # extract coordinates from gram matrix
    coordinates = []
    vals, vecs = npl.eigh(gram)

    vals = vals[N - 3 :]
    vecs = vecs.T[N - 3 :]

    # print('eigvals:', vals) # must all be positive for PSD (positive semidefinite) matrix

    # same eigenvalues might be small -> exact embedding does not exist
    # fix by replacing all but largest 3 eigvals by 0
    # better if three largest eigvals are separated by large spectral gap

    for val, vec in zip(vals, vecs):
        coord = vec * np.sqrt(val)
        coordinates.append(coord)

    return np.array(coordinates).T


def is_notebook():
    # sourcery skip: assign-if-exp, boolean-if-exp-identity, remove-unnecessary-cast, switch
    """
    This function detects if you are running on an ipython console or in the shell.
    It is used to either kill plots or leave them open.

    Returns
    -------
    TYPE Boolean
        true if running in Jupyter or Ipython consoles within spyder.
        false otherwise (terminal)

    """
    try:
        shell = get_ipython().__class__.__name__
        if shell == "ZMQInteractiveShell":
            return True  # Jupyter notebook or qtconsole
        elif shell == "TerminalInteractiveShell":
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False  # Probably standard Python interpreter


def plot_distance_histograms(
    sc_matrix_collated,
    pixel_size,
    output_filename="test",
    log_name_md="log.md",
    mode="hist",
    limit_n_plots=10,
    kernel_width=0.25,
    optimize_kernel_width=False,
    max_distance=4.0,
):
    if not is_notebook():
        n_plots_x = n_plots_y = sc_matrix_collated.shape[0]
    elif limit_n_plots == 0:
        n_plots_x = n_plots_y = sc_matrix_collated.shape[0]
    else:
        n_plots_x = n_plots_y = min(
            [limit_n_plots, sc_matrix_collated.shape[0]]
        )  # sets a max of subplots if you are outputting to screen!

    bins = np.arange(0, max_distance, 0.25)

    size_x, size_y = n_plots_x * 4, n_plots_y * 4

    fig, axs = plt.subplots(
        figsize=(size_x, size_y), ncols=n_plots_x, nrows=n_plots_y, sharex=True
    )

    for i in trange(n_plots_x):
        for j in range(n_plots_y):
            if i != j:
                if mode == "hist":
                    axs[i, j].hist(pixel_size * sc_matrix_collated[i, j, :], bins=bins)
                else:
                    (
                        max_kde,
                        distance_distribution,
                        kde,
                        x_d,
                    ) = distribution_maximum_kernel_density_estimation(
                        sc_matrix_collated,
                        i,
                        j,
                        pixel_size,
                        optimize_kernel_width=optimize_kernel_width,
                        kernel_width=kernel_width,
                        max_distance=max_distance,
                    )
                    axs[i, j].fill_between(x_d, kde, alpha=0.5)
                    axs[i, j].plot(
                        distance_distribution,
                        np.full_like(distance_distribution, -0.01),
                        "|k",
                        markeredgewidth=1,
                    )
                    axs[i, j].vlines(max_kde, 0, kde.max(), colors="r")

            axs[i, j].set_xlim(0, max_distance)
            axs[i, j].set_yticklabels([])

    plt.xlabel("distance, µm")
    plt.ylabel("counts")

    file_extension = output_filename.split(".")[-1]

    if len(file_extension) == 3:
        file_name = f"{output_filename}_PWDhistograms.{file_extension}"
    else:
        file_name = f"{output_filename}_PWDhistograms.png"

    print(f"Output figure: {file_name}\n")
    plt.savefig(file_name)

    if not is_notebook():
        plt.close()


def adjust_colorbar(cbar, pos, c_min, clim):
    """
    Adjusts the colorbar to display correct limits and tick values.

    Parameters:
    - cbar: The colorbar object from plt.colorbar().
    - pos: The image object returned by plt.imshow().
    - c_min: The minimum value to set on the colorbar.
    - clim: The maximum value to set on the colorbar.
    """
    # Set color limits for the colormap
    pos.set_clim(c_min, clim)

    # Get default ticks
    ticks = cbar.get_ticks()

    # Filter ticks to keep them inside [c_min, clim]
    # Safety margin: avoid selecting ticks that are too close to c_min or clim.
    # Prevents generating a huge number of ticks when c_min ≈ clim (especially with Matplotlib <3.7).
    delta = 1e-3
    ticks = [tick for tick in ticks if c_min + delta <= tick <= clim - delta]

    # Ensure c_min and clim are included in ticks
    if c_min not in ticks:
        ticks.insert(0, c_min)
    if clim not in ticks:
        ticks.append(clim)

    # Apply new ticks to colorbar
    cbar.set_ticks(ticks)

    # Format tick labels with two decimals, keeping exact c_min and clim
    yticklabels = [f"{tick:.2f}" for tick in ticks[1:-1]]  # Middle ticks formatted
    yticklabels.insert(0, str(c_min))  # Exact c_min
    yticklabels.append(str(clim))  # Exact clim

    # Apply formatted labels
    cbar.ax.set_yticklabels(yticklabels)


def plot_single_matrix(
    matrix, cmap, matrix_title, fontsize, barcode_names, cm_title, c_min, clim, fig_path
):
    plt.figure(figsize=(15, 15))
    pos = plt.imshow(matrix, cmap=cmap)
    plt.title(matrix_title, fontsize=float(fontsize) * 1.3)
    plt.xlabel("barcode #", fontsize=float(fontsize) * 1.2)
    plt.ylabel("barcode #", fontsize=float(fontsize) * 1.2)
    plt.xticks(np.arange(len(barcode_names)), barcode_names, fontsize=fontsize)
    plt.yticks(np.arange(len(barcode_names)), barcode_names, fontsize=fontsize)
    cbar = plt.colorbar(pos, fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize=float(fontsize) * 0.8)
    cbar.minorticks_on()
    cbar.set_label(cm_title, fontsize=float(fontsize) * 1.0)
    adjust_colorbar(cbar, pos, c_min, clim)
    plt.savefig(fig_path)


def get_matrix_title(
    fig_type,
    n_barcodes,
    n_cells,
    proximity_threshold=None,
    remove_nan=False,
    c_min=-1,
    clim=0,
):
    title = f"{fig_type} | Barcodes: {n_barcodes} | Traces: {str(n_cells)}\n"
    remove_nan_txt = "nonNANs" if remove_nan else "n_cells"
    if "proximity" in fig_type:
        title += f"Threshold: {proximity_threshold}μm | norm: {remove_nan_txt} |"
    c_min_txt = "auto" if c_min == -1 else str(c_min)
    c_max_txt = "auto" if clim == 0 else str(clim)
    title += f" c_min: {c_min_txt} | c_max: {c_max_txt}\n"
    return title


def adjust_cmin_cmax(matrix, cmin=-1, cmax=0):
    if cmax == 0:
        cmax = round(np.nanmax(matrix), 4)
    if cmin == -1:
        arr = matrix.astype("float")
        arr[arr == 0] = np.nan
        cmin = round(np.nanmin(arr), 4)
    return cmin, cmax


def get_plot_path(
    input_filename,
    output_folder,
    file_format,
    c_min,
    c_max,
    mode,
    threshold=0.25,
    remove_nan=False,
):
    input_basename = os.path.basename(input_filename).split(".")[0]
    norm_txt = "_norm" if remove_nan else ""
    threshold_txt = f"_T{threshold}" if threshold != 0.25 else ""
    c_min_max_txt = f"{c_min:.2f}-{c_max:.2f}"
    filename = f"Fig_{input_basename}_{mode}{norm_txt}{threshold_txt}_{c_min_max_txt}.{file_format}"
    out_path = os.path.join(output_folder, filename)
    return out_path


def plot_nan_matrix(
    nan_matrix,
    barcode_names,
    input_filename,
    output_folder,
    file_format="png",
    n_cells=0,
    font_size=22,
    remove_nan=False,
):
    c_min, c_max = adjust_cmin_cmax(nan_matrix)
    m_title = get_matrix_title(
        "NaN percentage", nan_matrix.shape[0], n_cells, c_min=c_min, clim=c_max
    )
    plot_path = get_plot_path(
        input_filename,
        output_folder,
        file_format,
        c_min,
        c_max,
        mode="nan%",
        remove_nan=remove_nan,
    )
    plot_single_matrix(
        nan_matrix,
        "Reds",
        m_title,
        font_size,
        barcode_names,
        "NaN value (%)",
        c_min,
        c_max,
        plot_path,
    )


def plot_him_matrix(
    matrix,
    barcode_names,
    input_filename,
    output_folder,
    file_format="png",
    mode="proximity",
    n_cells=0,
    font_size=22,
    proximity_threshold=0.25,
    remove_nan=False,
    cmtitle="proximity frequency",
    c_min=-1,
    c_max=0,
    c_m="seismic",
):
    m_title = get_matrix_title(
        f"Mode: {mode}",
        str(matrix.shape[0]),
        n_cells,
        proximity_threshold,
        remove_nan,
        c_min,
        c_max,
    )
    c_min, c_max = adjust_cmin_cmax(matrix, c_min, c_max)
    plot_path = get_plot_path(
        input_filename,
        output_folder,
        file_format,
        c_min,
        c_max,
        mode,
        proximity_threshold,
        remove_nan,
    )
    plot_single_matrix(
        matrix,
        c_m,
        m_title,
        font_size,
        barcode_names,
        cmtitle,
        c_min,
        c_max,
        plot_path,
    )
    return plot_path[:-4]


def calculate_nan_matrix(sc_matrices):
    n_barcodes = sc_matrices.shape[0]
    n_cells = sc_matrices.shape[2]
    nan_matrix = np.zeros((n_barcodes, n_barcodes))
    for i in range(n_barcodes):
        for j in range(n_barcodes):
            if i != j:
                distance_distribution = sc_matrices[i, j, :]
                number_nans = len(np.nonzero(np.isnan(distance_distribution))[0])
                nan_percentage = number_nans / n_cells
                nan_matrix[i, j] = nan_percentage
    return nan_matrix


def calculate_contact_probability_matrix(
    i_sc_matrix_collated,
    pixel_size,
    threshold=0.25,
    remove_nan=False,
    min_number_contacts=0,
):
    n_x = n_y = i_sc_matrix_collated.shape[0]
    n_cells = i_sc_matrix_collated.shape[2]
    sc_matrix = np.zeros((n_x, n_y))
    for i in range(n_x):
        for j in range(n_y):
            if i != j:
                distance_distribution = pixel_size * i_sc_matrix_collated[i, j, :]
                number_contacts = distance_distribution.squeeze().shape[0] - len(
                    np.nonzero(np.isnan(distance_distribution))[0]
                )
                number_nans = len(np.nonzero(np.isnan(distance_distribution))[0])
                if number_contacts < min_number_contacts:
                    print(
                        f"$ Rejected {i}-{j} because number contacts: {number_contacts} < {min_number_contacts}"
                    )
                    probability = 0.0
                elif not remove_nan:
                    probability = (
                        len(np.nonzero(distance_distribution < threshold)[0]) / n_cells
                    )
                else:
                    if n_cells == number_nans:
                        probability = np.nan
                    else:
                        n_real_bin = n_cells - number_nans
                        below_threshold_indices = np.nonzero(
                            distance_distribution < threshold
                        )[0]
                        probability = len(below_threshold_indices) / n_real_bin
                sc_matrix[i, j] = probability
    return sc_matrix


# @jit(nopython=True)
def find_optimal_kernel_width(distance_distribution):
    bandwidths = 10 ** np.linspace(-1, 1, 100)
    grid = GridSearchCV(
        KernelDensity(kernel="gaussian"), {"bandwidth": bandwidths}, cv=LeaveOneOut()
    )
    grid.fit(distance_distribution[:, None])
    return grid.best_params_


# @jit(nopython=True)
def retrieve_kernel_density_estimator(
    distance_distribution_0, x_d, optimize_kernel_width=False, kernel_width=0.25
):
    """
    Gets the kernel density function and maximum from a distribution of PWD distances

    Parameters
    ----------
    distance_distribution_0 : nd array
        List of PWD distances.
    x_d : nd array
        x grid.
    optimize_kernel_width : Boolean, optional
        whether to optimize bandwidth. The default is False.

    Returns
    -------
    np array
        kde distribution.
    np array
        Original distribution without NaNs

    """

    nan_array = np.isnan(distance_distribution_0)

    not_nan_array = ~nan_array

    distance_distribution = distance_distribution_0[not_nan_array]

    # instantiate and fit the KDE model
    if optimize_kernel_width:
        kernel_width = find_optimal_kernel_width(distance_distribution)["bandwidth"]
    else:
        kernel_width = kernel_width

    kde = KernelDensity(bandwidth=kernel_width, kernel="gaussian")

    # makes sure the list is not full of NaNs.
    if distance_distribution.shape[0] > 0:
        kde.fit(distance_distribution[:, None])
    else:
        return np.array([0]), np.array([0])

    # score_samples returns the log of the probability density
    logprob = kde.score_samples(x_d[:, None])

    return logprob, distance_distribution


# @jit(nopython=True)
def distribution_maximum_kernel_density_estimation(
    sc_matrix_collated,
    bin1,
    bin2,
    pixel_size,
    optimize_kernel_width=False,
    kernel_width=0.25,
    max_distance=4.0,
):
    """
    calculates the kernel distribution and its maximum from a set of PWD distances

    Parameters
    ----------
    sc_matrix_collated : np array 3 dims
        SC PWD matrix.
    bin1 : int
        first bin.
    bin2 : int
        first bin.
    pixel_size : float
        pixel size in µm
    optimize_kernel_width : Boolean, optional
        does kernel need optimization?. The default is False.

    Returns
    -------
    float
        maximum of kernel.
    np array
        list of PWD distances used.
    np array
        kernel distribution.
    x_d : np array
        x grid.

    """
    distance_distribution_unlimited = (
        pixel_size * sc_matrix_collated[bin1, bin2, :]
    )  # full distribution
    distance_distribution_unlimited = distance_distribution_unlimited[
        ~np.isnan(distance_distribution_unlimited)
    ]  # removes nans

    if bin1 == bin2:
        # protection against bins in the diagonal
        distance_distribution_0 = distance_distribution_unlimited
    else:
        # removes values larger than max_distance
        distance_distribution_0 = distance_distribution_unlimited[
            np.nonzero(distance_distribution_unlimited < max_distance)
        ]
    x_d = np.linspace(0, max_distance, 2000)

    # checks that distribution is not empty
    if distance_distribution_0.shape[0] > 0:
        logprob, distance_distribution = retrieve_kernel_density_estimator(
            distance_distribution_0, x_d, optimize_kernel_width, kernel_width
        )
        if logprob.shape[0] > 1:
            kernel_distribution = 10 * np.exp(logprob)
            maximum_kernel_distribution = x_d[np.argmax(kernel_distribution)]
            return (
                maximum_kernel_distribution,
                distance_distribution,
                kernel_distribution,
                x_d,
            )
        else:
            return np.nan, np.zeros(x_d.shape[0]), np.zeros(x_d.shape[0]), x_d
    else:
        return np.nan, np.zeros(x_d.shape[0]), np.zeros(x_d.shape[0]), x_d


def get_rg_from_pwd(pwd_matrix_0, min_number_pwd=4, threshold=6):
    """
    Calculates the Rg from a 2D pairwise distance matrix
    while taking into account that some of the PWD might be NaN

    PWDmatrix:       numpy array, NxN
    minFracNotNaN:   require a minimal fraction of PWDs to be not NaN, return NaN otherwise

    for the math, see https://en.wikipedia.org/wiki/Radius_of_gyration#Molecular_applications
    """

    pwd_matrix = pwd_matrix_0.copy()

    # check that pwd_matrix is of right shape
    if pwd_matrix.ndim != 2:
        raise SystemExit(
            f"get_rg_from_pwd: Expected 2D input but got {pwd_matrix.ndim}D."
        )
    if pwd_matrix.shape[0] != pwd_matrix.shape[1]:
        raise SystemExit("get_rg_from_pwd: Expected square matrix as input.")

    # make sure the diagonal is NaN
    np.fill_diagonal(pwd_matrix, np.NaN)

    # filters out PWD
    pwd_matrix[pwd_matrix > threshold] = np.nan

    # get the number of PWDs that are not NaN
    num_not_nan = (
        np.sum(~np.isnan(pwd_matrix)) / 2
    )  # default is to compute the sum of the flattened array

    if num_not_nan < min_number_pwd:
        return np.NaN

    # calculate Rg
    sqr = np.square(pwd_matrix)
    sqr = np.nansum(sqr)  # default is to compute the sum of the flattened array

    rg_sq = sqr / (2 * (2 * num_not_nan + pwd_matrix.shape[0]))  # replaces 1/(2*N^2)

    return np.sqrt(rg_sq)


def get_detection_eff_barcodes(sc_matrix_collated):
    """
    Return the detection efficiency of all barcodes.
    Assumes a barcode is detected as soon as one PWD with this barcode is detected.
    """

    # check that pwd_matrix is of right shape
    if sc_matrix_collated.ndim != 3:
        raise SystemExit(
            f"getBarcodeEff: Expected 3D input but got {sc_matrix_collated.ndim}D."
        )
    if sc_matrix_collated.shape[0] != sc_matrix_collated.shape[1]:
        raise SystemExit(
            "getBarcodeEff: Expected axis 0 and 1 to have the same length."
        )

    # make sure the diagonal is NaN
    for i in range(sc_matrix_collated.shape[0]):
        sc_matrix_collated[i, i, :] = np.NaN

    # calculate barcode efficiency
    n_cells = sc_matrix_collated.shape[2]

    eff = np.sum(~np.isnan(sc_matrix_collated), axis=0)
    eff[eff > 1] = 1

    eff0 = eff.copy()
    n_cells_2 = np.nonzero(np.sum(eff0, axis=0) > 2)[0].shape[0]

    eff = np.sum(eff, axis=-1)  # sum over all cells

    eff = eff / n_cells_2

    print(f"\n\n *** n_cells={n_cells} | n_cells_2={n_cells_2}")
    return eff


def get_barcodes_per_cell(sc_matrix_collated):
    """
    Returns the number of barcodes that were detected in each cell of sc_matrix_collated.
    """

    # make sure the diagonal is NaN
    for i in range(sc_matrix_collated.shape[0]):
        sc_matrix_collated[i, i, :] = np.NaN

    num_barcodes = np.sum(~np.isnan(sc_matrix_collated), axis=0)
    num_barcodes[num_barcodes > 1] = 1
    num_barcodes = np.sum(num_barcodes, axis=0)

    return num_barcodes


def get_coordinates_from_pwd_matrix(matrix):
    # multi-dimensional scaling to get coordinates from PWDs
    # make sure mean_sc_matrix is symmetric
    matrix = 0.5 * (matrix + np.transpose(matrix))
    # run metric mds
    verbosity = 0  # default: 0, quite verbose: 2
    mds = manifold.MDS(
        n_components=3,
        metric=True,
        n_init=20,
        max_iter=3000,
        verbose=verbosity,
        eps=1e-9,
        n_jobs=1,
        random_state=1,
        dissimilarity="precomputed",  # euclidean | precomputed
    )

    return mds.fit(matrix).embedding_


def sort_cells_by_number_pwd(him_data):
    # sc_matrix = him_data.data["SCmatrixCollated"]
    sc_matrix = him_data.sc_matrix_selected

    n_cells = sc_matrix.shape[2]

    # finds the number of barcodes detected per cell.
    n_barcode_per_cell = []
    values = []
    dtype = [("cellID", int), ("n_pwd", int)]

    for i_cell in range(n_cells):
        sc_matrix_cell = sc_matrix[:, :, i_cell]
        n_pwd = int(np.count_nonzero(~np.isnan(sc_matrix_cell)) / 2)
        n_barcode_per_cell.append(n_pwd)
        values.append((i_cell, n_pwd))

    values_array = np.array(values, dtype=dtype)  # create a structured array
    sorted_values = np.sort(values_array, order="n_pwd")

    return sc_matrix, sorted_values, n_cells


def kde_fit(x, x_d, bandwidth=0.2, kernel="gaussian"):
    kde = KernelDensity(bandwidth=bandwidth, kernel="gaussian")
    kde.fit(x[:, None])

    logprob = kde.score_samples(x_d[:, None])

    return logprob, kde


def calculate_ensemble_pwd_matrix(sc_matrix, pixel_size, cells_to_plot, mode="median"):
    """
    performs a KDE or median to calculate the max of the PWD distribution

    Parameters
    ----------
    sc_matrix : TYPE
        DESCRIPTION.
    pixel_size : TYPE
        DESCRIPTION.

    Returns
    -------
    matrix = 2D npy array.

    """

    n_barcodes = sc_matrix.shape[0]
    # cells_to_plot = range(sc_matrix.shape[2])

    mean_sc_matrix = np.zeros((n_barcodes, n_barcodes))

    if mode == "median":
        # calculates the median of all values #
        #######################################
        if max(cells_to_plot) > sc_matrix.shape[2]:
            print(
                f"Error with range in cells2plot {max(cells_to_plot)} as it is larger \
                    than the number of available cells {sc_matrix.shape[2]}"
            )

            keep_plotting = False
        else:
            mean_sc_matrix = pixel_size * np.nanmedian(
                sc_matrix[:, :, cells_to_plot], axis=2
            )
            keep_plotting = True

    elif mode == "KDE":
        keep_plotting = True

        if max(cells_to_plot) > sc_matrix.shape[2]:
            print(
                f"Error with range in cells2plot {max(cells_to_plot)} as it is larger \
                    than the number of available cells {sc_matrix.shape[2]}"
            )

            keep_plotting = False
        else:
            for bin1 in trange(n_barcodes):
                for bin2 in range(n_barcodes):
                    if bin1 != bin2:
                        # print(f"cells_to_plot:{cells_to_plot}, ncells:{sc_matrix.shape}")
                        (
                            maximum_kernel_distribution,
                            _,
                            _,
                            _,
                        ) = distribution_maximum_kernel_density_estimation(
                            sc_matrix[:, :, cells_to_plot],
                            bin1,
                            bin2,
                            pixel_size,
                            optimize_kernel_width=False,
                        )
                        mean_sc_matrix[bin1, bin2] = maximum_kernel_distribution

    return mean_sc_matrix, keep_plotting
