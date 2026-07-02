#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This class will contain methods to load, save, plot barcode localizations and statistics
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.table import Table, vstack


class LocalizationTable:
    def __init__(self):
        self.a = 1
        self.original_format = "ecsv"  # Default format
        self.fof_ct_version = "v1.0"
        self.table_namespace = "4dn_FOF-CT_demultiplexing"
        self.xyz_unit = "pixel"
        self.lab_name = "NollmannLab"
        self.experimenter_name = ""
        self.experimenter_contact = ""
        self.description = ""
        self.software_title = ""
        self.software_type = ""
        self.software_authors = ""
        self.software_description = ""
        self.software_repository = ""
        self.software_citation = ""
        self.additional_tables = "4dn_FOF-CT_core"
        self.columns = [
            "Loc_ID",
            "ROI #",
            "CellID #",
            "Barcode #",
            "Spot_ID",
            "Z",
            "X",
            "Y",
            "snr",
            "spot_pixel_percentage",
            "skew",
            "patch_size",
            "object_class",
            "mean_intensity",
            "flux",
            "roundness",
        ]

    def _read_metadata_from_4dn(self, file):
        """
        Reads metadata fields from a .4dn file and stores them as class attributes.
        """
        with open(file, "r") as f:
            for line in f:
                if line.startswith("##FOF-CT_Version="):
                    self.fof_ct_version = line.split("=")[1].strip()
                elif line.startswith("#Table_Namespace="):
                    self.table_namespace = line.split("=")[1].strip()
                elif line.startswith("#XYZ_Unit="):
                    self.xyz_unit = line.split("=")[1].strip()
                elif line.startswith("#Lab_Name:"):
                    self.lab_name = line.split(": ")[1].strip()
                elif line.startswith("#Experimenter_Name:"):
                    self.experimenter_name = line.split(": ")[1].strip()
                elif line.startswith("#Experimenter_Contact:"):
                    self.experimenter_contact = line.split(": ")[1].strip()
                elif line.startswith("#Description:"):
                    self.description = line.split(": ")[1].strip()
                elif line.startswith("#Software_Type:"):
                    self.software_type = line.split(": ")[1].strip()

                elif line.startswith("#Software_Title:"):
                    self.software_title = line.split(": ")[1].strip()
                elif line.startswith("#Software_Authors:"):
                    self.software_authors = line.split(": ")[1].strip()
                elif line.startswith("#Software_Description:"):
                    self.software_description = line.split(": ")[1].strip()
                elif line.startswith("#Software_Repository:"):
                    self.software_repository = line.split(": ")[1].strip()
                elif line.startswith("#Software_PreferredCitationID:"):
                    self.software_citation = line.split(": ")[1].strip()
                elif line.startswith("#Additional_Tables:"):
                    self.additional_tables = line.split(": ")[1].strip()
                elif line.startswith("##columns="):
                    self.columns = line.split("=")[1].split("(")[1].split(")")[0]
                    self.columns = self.columns.split(", ")
                    # self.columns = self._read_column_names_from_4dn(file)
                    print(f"> Columns read: {self.columns}")

    def _convert_4dn_to_astropy(self, fofct_file):
        """
        Converts a .4dn file to an Astropy table with appropriate formatting.
        Also saves a BED file mapping genomic coordinates to barcode numbers.
        """
        column_names = self.columns
        csv_data = pd.read_csv(fofct_file, comment="#", header=None, names=column_names)

        # Rename columns for Astropy compatibility
        csv_data.rename(
            columns={
                "Loc_ID": "Buid",
                "Spot_ID": "id",
                "X": "xcentroid",
                "Y": "ycentroid",
                "Z": "zcentroid",
            },
            inplace=True,
        )

        return Table.from_pandas(csv_data)

    def _convert_astropy_to_4dn(self, table, output_file):
        """
        Converts an Astropy table back to .4dn format with appropriate headers.
        """
        output_file = output_file.strip(".ecsv").strip(".4dn") + ".4dn"

        csv_data = table.to_pandas()
        csv_data.rename(
            columns={
                "Buid": "Loc_ID",
                "id": "Spot_ID",
                "xcentroid": "X",
                "ycentroid": "Y",
                "zcentroid": "Z",
            },
            inplace=True,
        )

        # parses column list for header
        column_list = ", ".join(self.columns)

        header = f"""##FOF-CT_Version={self.fof_ct_version}
##Table_Namespace={self.table_namespace}
##XYZ_Unit={self.xyz_unit}
#Lab_Name: {self.lab_name}
#Experimenter_Name: {self.experimenter_name}
#Experimenter_Contact: {self.experimenter_contact}
#Description: {self.description}
#Software_Title: {self.software_title}
#Software_Type: {self.software_type}
#Software_Authors: {self.software_authors}
#Software_Description: {self.software_description}
#Software_Repository: {self.software_repository}
#Software_PreferredCitationID: {self.software_citation}
#Additional_Tables: {self.additional_tables}
##Columns=({column_list})
"""
        with open(output_file, "w") as f:
            f.write(header)
            csv_data.to_csv(f, index=False, header=False)
        print(f"Saved 4dn spot table with headers: {output_file}")

    def load(self, file):
        """
        Loads barcode_map

        Parameters
        ----------
        filename_barcode_coordinates : string
            filename with barcode_map

        Returns
        -------
        barcode_map : Table()
        unique_barcodes: list
            lis of unique barcodes read from barcode_map

        """

        if not os.path.exists(file):
            print(f"# ERROR: could not find coordinates file: {file}")
            sys.exit()

        file_ext = os.path.splitext(file)[1].lower()
        if file_ext in (".ecsv", ".dat"):
            # print("$ Importing table from pyHiM format")
            barcode_map = read_table_from_ecsv(file)
            self.data = barcode_map
            self.original_format = "ecsv"
        elif file_ext == ".4dn":
            print("$ Importing table from fof-ct format")
            self._read_metadata_from_4dn(file)
            barcode_map = self._convert_4dn_to_astropy(file)
            self.data = barcode_map
            self.original_format = "4dn"
        else:
            raise ValueError("Unsupported file format. Use .ecsv, .dat, or .4dn")

        # print(f"$ Successfully loaded barcode localizations file: {file}")

        unique_barcodes = np.unique(barcode_map["Barcode #"].data)
        # number_unique_barcodes = unique_barcodes.shape[0]

        # print(f"$ Number of barcodes read from barcode_map: {number_unique_barcodes}")
        # print(f"$ Unique Barcodes detected: {unique_barcodes}")

        return barcode_map, unique_barcodes

    def remove_empty_comments(self):
        try:
            if len(self.data.meta["comments"]):
                self.data.meta["comments"] = [
                    com for com in self.data.meta["comments"] if com
                ]
        except KeyError:
            self.data.meta["comments"] = []

    def remove_duplicate_comments(self):
        if len(self.data.meta["comments"]):
            self.data.meta["comments"] = list(dict.fromkeys(self.data.meta["comments"]))

    def save(self, file_name, barcode_map, comments="", format="ecsv"):
        """
        Saves output table

        Parameters
        ----------
        filename_barcode_coordinates : string
            filename of table.
        barcode_map : astropy Table
            Table to be written to file.
        tag : string, optional
            tag to be added to filename. The default is "_".
        ext : string, optional
            file extension. The default is 'ecsv'.
        comments : list of strings, optional
            Will output as comments to the header. The default is [].

        Returns
        -------
        None.

        """
        if format == "4dn":
            self._convert_astropy_to_4dn(barcode_map, file_name)
        else:
            print(f"$ Saving output table as {file_name} ...")
            self.data = barcode_map
            self.remove_empty_comments()
            self.remove_duplicate_comments()
            try:
                barcode_map.meta["comments"].append(comments)
            except KeyError:
                barcode_map.meta["comments"] = [comments]

            barcode_map.write(
                file_name,
                format="ascii.ecsv",
                overwrite=True,
            )

    def append(self, table1, table2):
        """
        appends <table> to self.data

        Parameters
        ----------
        table : astropy table
            table to append to existing self.data table.

        Returns
        -------
        None.

        """

        return vstack([table1, table2])

    def plot_distribution_fluxes(
        self, barcode_map, filename_list=("localization_distribution_fluxes.png",)
    ):
        """
        This function will plot:
        - the number of localizations per barcode
        - the snr distribution per barcode
        - scatterplot of the snr versus z
        - scatterplot of roundness versus skew

        Parameters
        ----------
        barcode_map : TYPE
            DESCRIPTION.
        filename_list: list
            filename

        Returns
        -------
        None.

        """
        from matplotlib.colors import BoundaryNorm

        # initializes figure and font settings explicitly so plots look the same
        # whether this method is called from pyHiM, notebooks, or the CLI.
        figure_size = (30, 15)
        axes_label_size = 24
        tick_label_size = 20
        colorbar_label_size = 24
        save_dpi = 100

        fig, axes = plt.subplots(2, 2, figsize=figure_size)
        ax = axes.ravel()

        # initializes variables
        skew = barcode_map["skew"]
        barcode_id = barcode_map["Barcode #"]
        zcentroid = barcode_map["zcentroid"]
        snr = barcode_map["snr"]
        object_class = barcode_map["object_class"]
        roundness = barcode_map["roundness"]

        # plots data
        # panel 1

        # Sort barcode identities
        unique_barcodes = np.sort(np.unique(barcode_id))

        # Collect SNR values for each barcode
        snr_by_barcode = [snr[barcode_id == bc] for bc in unique_barcodes]

        # Draw violin plot
        parts = ax[0].violinplot(
            snr_by_barcode,
            positions=unique_barcodes,
            widths=0.8,
            showmeans=False,
            showmedians=True,
            showextrema=True,
        )

        ax[0].set_xlabel("barcode_id", fontsize=axes_label_size)
        ax[0].set_ylabel("snr", fontsize=axes_label_size)
        ax[0].set_xticks(unique_barcodes)

        # panel 2
        for body in parts["bodies"]:
            body.set_facecolor("steelblue")
            body.set_edgecolor("black")
            body.set_alpha(0.7)

        parts["cmedians"].set_color("red")
        parts["cbars"].set_color("black")
        parts["cmins"].set_color("black")
        parts["cmaxes"].set_color("black")

        p_2 = ax[1].scatter(snr, zcentroid, c=object_class, cmap="seismic", alpha=0.55)
        ax[1].set_xlabel("snr", fontsize=axes_label_size)
        ax[1].set_ylabel("z_centroid", fontsize=axes_label_size)

        cbar2 = fig.colorbar(
            p_2,
            ax=ax[1],
            fraction=0.046,
            pad=0.04,
        )

        cbar2.set_label("object_class", fontsize=colorbar_label_size)
        cbar2.ax.tick_params(labelsize=tick_label_size)

        # panel 3
        unique_barcodes, counts = np.unique(barcode_id, return_counts=True)

        ax[2].bar(unique_barcodes, counts, width=0.8)
        ax[2].set_xlabel("barcode_id", fontsize=axes_label_size)
        ax[2].set_ylabel("Number of detections", fontsize=axes_label_size)
        ax[2].set_xticks(unique_barcodes)

        # panel 4
        unique_barcodes = np.sort(np.unique(barcode_id))

        cmap = plt.get_cmap("tab20b", len(unique_barcodes))
        norm = BoundaryNorm(np.arange(len(unique_barcodes) + 1) - 0.5, cmap.N)

        # Map barcode IDs to consecutive integers
        barcode_to_idx = {bc: i for i, bc in enumerate(unique_barcodes)}
        color_idx = np.array([barcode_to_idx[bc] for bc in barcode_id])

        p_3 = ax[3].scatter(
            roundness,
            skew,
            c=color_idx,
            cmap=cmap,
            norm=norm,
            alpha=0.5,
        )

        ax[3].set_ylabel("skew", fontsize=axes_label_size)
        ax[3].set_xlabel("roundness", fontsize=axes_label_size)

        cbar = fig.colorbar(
            p_3,
            ax=ax[3],
            ticks=np.arange(len(unique_barcodes)),
            fraction=0.046,
            pad=0.04,
        )

        cbar.set_label("Barcode", fontsize=colorbar_label_size)
        cbar.set_ticklabels(unique_barcodes)
        cbar.ax.tick_params(labelsize=tick_label_size)

        for axis in ax:
            axis.tick_params(axis="both", labelsize=tick_label_size)

        # saves figure
        fig.savefig("".join(filename_list), dpi=save_dpi)

        plt.close(fig)

    def plot_intensity_distribution(
        self, intensities, output_file="intensity_distribution.png"
    ):
        """
        Plots the distribution of intensity values from the localization table
        with a logarithmic y-axis and saves it as a PNG file.
        """

        plt.figure(figsize=(8, 6))
        plt.hist(intensities, bins=50, color="blue", alpha=0.7, edgecolor="black")
        plt.xlabel("Intensity")
        plt.ylabel("Frequency (log scale)")
        plt.title("Localization Intensity Distribution")
        plt.yscale("log")  # Set y-axis to logarithmic scale
        plt.grid(True)
        plt.savefig(output_file)
        plt.close()
        print(f"$ Saved intensity distribution plot as {output_file}")

    def plots_localizations(self, barcode_map_full, filename_list):
        """
        This function plots 3 subplots (xy, xz, yz) with the localizations.
        One figure is produced per ROI.

        Parameters
        ----------
        image : List of numpy ndarray (N-dimensional array)
            3D raw image of format .tif

        label : List of numpy ndarray (N-dimensional array)
            3D labeled image of format .tif

        filename_list: list
            filename
        """

        # indexes table by ROI
        barcode_map_roi, number_rois = decode_rois(barcode_map_full)

        for i_roi in range(number_rois):
            # creates sub Table for this ROI
            barcode_map = barcode_map_roi.groups[i_roi]
            n_roi = barcode_map["ROI #"][0]
            print(f"> Plotting barcode localization map for ROI: {n_roi}")
            color_dict = build_color_dict(barcode_map, key="Barcode #")

            # initializes figure
            fig = plt.figure(constrained_layout=False)
            im_size = 60
            fig.set_size_inches((im_size * 2, im_size))
            gs = fig.add_gridspec(2, 2)
            ax = [
                fig.add_subplot(gs[:, 0]),
                fig.add_subplot(gs[0, 1]),
                fig.add_subplot(gs[1, 1]),
            ]

            # defines variables
            x = barcode_map["xcentroid"]
            y = barcode_map["ycentroid"]
            z = barcode_map["zcentroid"]
            colors = [color_dict[str(x)] for x in barcode_map["Barcode #"]]
            titles = [
                "Z-projection (pixel)",
                "X-projection (pixel)",
                "Y-projection (pixel)",
            ]

            # makes plot
            plots_localization_projection(x, y, ax[0], colors, titles[0])
            plots_localization_projection(x, z, ax[1], colors, titles[1])
            plots_localization_projection(y, z, ax[2], colors, titles[2])

            fig.tight_layout()

            # saves output figure
            filename_list_i = filename_list.copy()
            filename_list_i.insert(-1, f"_ROI{str(n_roi)}")
            fig.savefig("".join(filename_list_i))

    def compares_localizations(
        self, barcode_map_1, barcode_map_2, filename_list, fontsize=20
    ):
        """
        Compares the localizations of two barcode tables

        Parameters
        ----------
        barcode_map_1 : astropy Table
            localization table 1.
        barcode_map_2 : astropy Table
            localization table 2.

        Returns
        -------
        None.

        """

        barcode_map_2.add_index("Buid")
        number_localizations = len(barcode_map_1)

        labels = ["xcentroid", "ycentroid", "zcentroid"]
        diffs = {label: [] for label in labels}
        # iterates over rows in barcode_map_1
        for row in range(number_localizations):
            buid_1 = barcode_map_2[row]["Buid"]
            barcode_found = True

            # finds same Buid in barcode_map_2
            try:
                barcode_map_2.loc[buid_1]
            except KeyError:
                barcode_found = False

            # collects differences in values between same localization in both tables
            if barcode_found:
                for label in labels:
                    a, b = barcode_map_2.loc[buid_1][label], barcode_map_1[row][label]
                    if ~np.isnan(a).any() and ~np.isnan(b).any():
                        # diff = a - b
                        # if np.isnan(diff):
                        #    diff = 0
                        diffs[label].append(a - b)

        # plots figures
        fig, axes = plt.subplots(2, 2)
        ax = axes.ravel()
        fig.set_size_inches((30, 30))

        for label, axis in zip(labels, ax):
            r = np.array(diffs[label])
            axis.hist(r, bins=20)
            axis.set_xlabel(f"{label} correction, px", fontsize=fontsize)
            axis.set_ylabel("counts", fontsize=fontsize)

        ax[3].scatter(
            np.array(diffs["ycentroid"]), np.array(diffs["xcentroid"]), s=3, alpha=0.8
        )
        ax[3].set_xlabel("dx-position, px", fontsize=fontsize)
        ax[3].set_ylabel("dy-position, px", fontsize=fontsize)

        fig.savefig("".join(filename_list))


def decode_rois(data):
    data_indexed = data.group_by("ROI #")

    number_rois = len(data_indexed.groups.keys)

    print(f"\n$ rois detected: {number_rois}")

    return data_indexed, number_rois


def build_color_dict(data, key="Barcode #"):
    unique_barcodes = np.unique(data[key])
    output_array = range(unique_barcodes.shape[0])

    return {
        str(barcode): output for barcode, output in zip(unique_barcodes, output_array)
    }


def plots_localization_projection(coord1, coord2, axis, colors, title="" * 3):
    """
    This function will produce the scatter plot and add title

    Parameters
    ----------
    coord1 : 1D Numpy array, float
        first coordinate (x).
    coord2 : 1D Numpy array, float
        first coordinate (y).
    axis : matplotlib axis
        figure axis handle.
    colors : 1D Numpy array, float
        colorcode used in scatter plot.
    title : string, optional
        title of subpanel. The default is ''*3.

    Returns
    -------
    None.

    """
    axis.scatter(coord1, coord2, s=5, c=colors, alpha=0.9, cmap="hsv")  # nipy_spectral
    axis.set_title(title)


def read_table_from_ecsv(path):
    """Read an astropy Table saved as an ``ecsv`` file."""
    # read ecsv file
    table = Table.read(path, format="ascii.ecsv")

    return table


def create_output_table():
    output = Table(
        names=(
            "Buid",
            "ROI #",
            "CellID #",
            "Barcode #",
            "id",
            "zcentroid",
            "xcentroid",
            "ycentroid",
            "snr",
            "spot_pixel_percentage",
            "skew",
            "patch_size",
            "object_class",
            "mean_intensity",
            "flux",
            "roundness",
        ),
        dtype=(
            "S2",
            "int",
            "int",
            "int",
            "int",
            "f4",
            "f4",
            "f4",
            "f4",
            "f4",
            "f4",
            "int",
            "int",
            "f4",
            "f4",
            "f4",
        ),
    )
    return output


"""
def create_output_table():
    output = Table(
        names=(
            "Buid",
            "ROI #",
            "CellID #",
            "Barcode #",
            "id",
            "zcentroid",
            "xcentroid",
            "ycentroid",
            "sharpness",
            "roundness1",
            "roundness2",
            "npix",
            "sky",
            "peak",
            "flux",
            "mag",
        ),
        dtype=(
            "S2",
            "int",
            "int",
            "int",
            "int",
            "f4",
            "f4",
            "f4",
            "f4",
            "f4",
            "f4",
            "int",
            "f4",
            "f4",
            "f4",
            "f4",
        ),
    )
    return output
"""
