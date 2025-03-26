# -*- coding: utf-8 -*-
"""
trace table management class
"""

import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.table import Table, vstack
from tqdm import tqdm

font = {"weight": "normal", "size": 22}
matplotlib.rc("font", **font)


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


class ChromatinTraceTable:

    def __init__(self, xyz_unit="micron", genome_assembly="mm10"):
        """
        Initializes the ChromatinTraceTable class.
        """
        self.a = 1
        self.xyz_unit = xyz_unit
        self.genome_assembly = genome_assembly
        self.experimenter_name = ""
        self.experimenter_contact = ""
        self.software_title = ""
        self.software_authors = ""
        self.lab_name = ""
        self.software_description = ""
        self.software_repository = ""
        self.software_citation = ""
        self.columns = []
        self.data = None
        self.original_format = "ecsv"  # Default format

    def initialize(self):
        self.data = Table(
            names=(
                "Spot_ID",
                "Trace_ID",
                "x",
                "y",
                "z",
                "Chrom",
                "Chrom_Start",
                "Chrom_End",
                "ROI #",
                "Mask_id",
                "Barcode #",
                "label",
            ),
            dtype=(
                "S2",
                "S2",
                "f4",
                "f4",
                "f4",
                "S2",
                "int",
                "int",
                "int",
                "int",
                "int",
                "S2",
            ),
        )

        self.data.meta["comments"] = [
            f"xyz_unit={self.xyz_unit}",
            f"genome_assembly={self.genome_assembly}",
        ]

    def load(self, file):
        """
        Loads a trace table from a .ecsv or .4dn file.
        """
        if not os.path.exists(file):
            print(f"# ERROR: could not find file: {file}")
            sys.exit()

        file_ext = os.path.splitext(file)[1].lower()
        if file_ext == ".ecsv":
            print("$ Importing table from pyHiM format")
            self.data = read_table_from_ecsv(file)
            self.original_format = "ecsv"
        elif file_ext == ".4dn":
            print("$ Importing table from fof-ct format")
            self._read_metadata_from_4dn(file)
            self.data = self._convert_4dn_to_astropy(file)
            self.original_format = "4dn"
        else:
            raise ValueError("Unsupported file format. Use .ecsv or .4dn")

        print(f"Successfully loaded trace table: {file}")
        return self.data

    def remove_empty_comments(self):
        if len(self.data.meta["comments"]):
            self.data.meta["comments"] = [
                com for com in self.data.meta["comments"] if com
            ]

    def remove_duplicate_comments(self):
        if len(self.data.meta["comments"]):
            self.data.meta["comments"] = list(dict.fromkeys(self.data.meta["comments"]))

    def save(self, file_name, comments=""):
        """
        Saves the trace table in the appropriate format (.ecsv or .4dn).
        """
        if self.original_format == "4dn":
            self._convert_astropy_to_4dn(self.data, file_name)
        else:
            print(f"$ Saving output table as {file_name} ...")
            self.remove_empty_comments()
            self.remove_duplicate_comments()
            try:
                self.data.meta["comments"].append(comments)
            except KeyError:
                self.data.meta["comments"] = [comments]

            save_table_to_ecsv(self.data, file_name)

    def _read_metadata_from_4dn(self, file):
        """
        Reads metadata fields from a .4dn file and stores them as class attributes.
        """
        with open(file, "r") as f:
            for line in f:
                if line.startswith("#experimenter_name:"):
                    self.experimenter_name = line.split(": ")[1].strip()
                elif line.startswith("#experimenter_contact:"):
                    self.experimenter_contact = line.split(": ")[1].strip()
                elif line.startswith("##genome_assembly:"):
                    self.genome_assembly = line.split("=")[1].strip()
                elif line.startswith("#Software_Title:"):
                    self.software_title = line.split(": ")[1].strip()
                elif line.startswith("#Software_Authors:"):
                    self.software_authors = line.split(": ")[1].strip()
                elif line.startswith("#lab_name:"):
                    self.lab_name = line.split(": ")[1].strip()
                elif line.startswith("#Software_Description:"):
                    self.software_description = line.split(": ")[1].strip()
                elif line.startswith("#Software_Repository:"):
                    self.software_repository = line.split(": ")[1].strip()
                elif line.startswith("#Software_PreferredCitationID:"):
                    self.software_citation = line.split(": ")[1].strip()
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
        column_names = self.columns  # self._read_column_names_from_4dn(fofct_file)
        csv_data = pd.read_csv(fofct_file, comment="#", header=None, names=column_names)

        # Rename XYZ columns for Astropy compatibility
        csv_data.rename(columns={"X": "x", "Y": "y", "Z": "z"}, inplace=True)

        # Handle optional Cell_ID column
        if "Cell_ID" in csv_data.columns:
            csv_data.rename(columns={"Cell_ID": "Mask_id"}, inplace=True)
        else:
            print("> No Cell_ID column found, will use -1 for Mask_id")
            csv_data["Mask_id"] = -1  # Placeholder if Cell_ID is missing

        # Handle optional Extra_Cell_ROI_ID column
        if "Extra_Cell_ROI_ID" in csv_data.columns:
            csv_data.rename(columns={"Extra_Cell_ROI_ID": "ROI #"}, inplace=True)
        else:
            print("> No Extra_Cell_ROI_ID column found, will use 0 for Mask_id")
            csv_data["ROI #"] = 0  # Default value if missing

        # Assign Barcode # by ordering and mapping unique genomic positions
        unique_barcodes = (
            csv_data[["Chrom", "Chrom_Start", "Chrom_End"]]
            .drop_duplicates()
            .sort_values(by=["Chrom", "Chrom_Start", "Chrom_End"])
            .reset_index(drop=True)
        )
        unique_barcodes["Barcode #"] = range(1, len(unique_barcodes) + 1)
        barcode_mapping = {
            tuple(row[:3]): row[3]
            for row in unique_barcodes.itertuples(index=False, name=None)
        }

        csv_data["Barcode #"] = csv_data.apply(
            lambda row: barcode_mapping[
                (row["Chrom"], row["Chrom_Start"], row["Chrom_End"])
            ],
            axis=1,
        )
        csv_data["label"] = "None"  # Placeholder for label

        # Save BED file with Barcode # mapping
        bed_file = fofct_file.replace(".4dn", ".bed")
        unique_barcodes.to_csv(bed_file, sep="\t", header=False, index=False)
        print(f"Saved BED file: {bed_file}")

        return Table.from_pandas(csv_data)

    def _convert_astropy_to_4dn(self, table, output_file):
        """
        Converts an Astropy table back to .4dn format with appropriate headers.
        """
        output_file = output_file.strip(".ecsv") + ".4dn"

        csv_data = table.to_pandas()
        csv_data.rename(columns={"x": "X", "y": "Y", "z": "Z"}, inplace=True)

        # Remove extra columns
        csv_data = csv_data.drop(columns=["Barcode #", "label"], errors="ignore")

        # Rename columns
        if "Extra_Cell_ROI_ID" in self.columns:
            csv_data.rename(columns={"ROI #": "Extra_Cell_ROI_ID"}, inplace=True)
        else:
            csv_data = csv_data.drop(columns=["ROI #"], errors="ignore")

        if "Cell_ID" in self.columns:
            csv_data.rename(columns={"Mask_id": "Cell_ID"}, inplace=True)
        else:
            csv_data = csv_data.drop(columns=["Mask_id"], errors="ignore")

        # parses column list for header
        column_list = ", ".join(self.columns)

        header = f"""##FOF-CT_version=v0.1
##Table_namespace=4dn_FOF-CT_core
##genome_assembly={self.genome_assembly}
##XYZ_unit=micron
#Software_Title: {self.software_title}
#Software_Type: SpotLoc+Tracing
#Software_Authors: {self.software_authors}
#Software_Description: {self.software_description}
#Software_Repository: {self.software_repository}
#Software_PreferredCitationID: {self.software_citation}
#lab_name: {self.lab_name}
#experimenter_name: {self.experimenter_name}
#experimenter_contact: {self.experimenter_contact}
#additional_tables:
##columns=({column_list})
"""
        # print(f"> Columnds to export to 4dn table: {csv_data.columns}")
        with open(output_file, "w") as f:
            f.write(header)
            csv_data.to_csv(f, index=False, header=False)
        print(f"Saved 4dn trace table with headers: {output_file}")

    def prevent_roi_conflict(self, table):
        existing_roi = set(self.data["ROI #"])
        conflicting_rois = set(table["ROI #"]) & existing_roi
        roi_mapping = {}
        current_max = max(existing_roi) + 1 if existing_roi else 1
        for roi in conflicting_rois:
            while current_max in existing_roi:
                current_max += 1
            roi_mapping[roi] = current_max
            existing_roi.add(current_max)
        new_roi_numbers = [roi_mapping.get(roi, roi) for roi in table["ROI #"]]
        table["ROI #"] = new_roi_numbers
        return table

    def load_bed_file(self, bed_file):
        """Loads the BED file into a dictionary mapping barcode numbers to genomic coordinates.
        Handles files with inconsistent tab spacing by using regex splitting."""
        bed_dict = {}

        with open(bed_file, "r") as f:
            for line in f:
                # Skip empty lines
                if not line.strip():
                    continue

                # Split on any number of whitespace characters
                # This handles inconsistent tabs/spaces more robustly
                fields = line.strip().split()

                # Ensure we have exactly 4 fields
                if len(fields) != 4:
                    print(f"Warning: Skipping malformed line: {line.strip()}")
                    continue

                try:
                    chrom = fields[0]
                    chrom_start = int(fields[1])
                    chrom_end = int(fields[2])
                    barcode = int(fields[3])

                    bed_dict[barcode] = {
                        "Chrom": chrom,
                        "Chrom_Start": chrom_start,
                        "Chrom_End": chrom_end,
                    }
                except ValueError as e:
                    print(
                        f"Warning: Skipping line with invalid data types: {line.strip()}"
                    )
                    print(f"Error: {e}")

        if not bed_dict:
            raise ValueError("No valid entries found in the BED file.")

        print(f"Successfully loaded {len(bed_dict)} barcode mappings from BED file.")

        return bed_dict

    def impute_genomic_coordinates(self, bed_dict, auto_continue=False):
        """Updates the Chrom, Chrom_Start, and Chrom_End columns in the trace file based on the BED file."""

        unmatched_barcodes = set()
        matched_count = 0
        total_count = 0

        for row in self.data:
            barcode = row["Barcode #"]
            total_count += 1

            if barcode in bed_dict:
                row["Chrom"] = bed_dict[barcode]["Chrom"]
                row["Chrom_Start"] = bed_dict[barcode]["Chrom_Start"]
                row["Chrom_End"] = bed_dict[barcode]["Chrom_End"]
                matched_count += 1
            else:
                unmatched_barcodes.add(barcode)

        if unmatched_barcodes:
            missing_percent = (len(unmatched_barcodes) / total_count) * 100
            print(
                f"Warning: {len(unmatched_barcodes)} unique barcodes ({missing_percent:.1f}% of rows) couldn't be matched:"
            )

            # Only show up to 10 unmatched barcodes to avoid cluttering the output
            if len(unmatched_barcodes) <= 10:
                print(
                    f"  Unmatched barcodes: {', '.join(map(str, sorted(unmatched_barcodes)))}"
                )
            else:
                print(
                    f"  First 10 unmatched barcodes: {', '.join(map(str, sorted(list(unmatched_barcodes))[:10]))}"
                )
                print(f"  ... and {len(unmatched_barcodes) - 10} more")

            # Ask the user if they want to continue if more than 10% of barcodes are unmatched
            if missing_percent > 10 and ~auto_continue:
                response = input(
                    "More than 10% of barcodes couldn't be matched. Continue anyway? (y/n): "
                )
                if response.lower() != "y":
                    print("Operation aborted by user.")
                    return

        print(
            f"Successfully matched {matched_count}/{total_count} rows ({(matched_count/total_count)*100:.1f}%)"
        )

        return matched_count, total_count

    def append(self, table):
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
        table = self.prevent_roi_conflict(table)
        self.data = vstack([self.data, table])

    def filter_traces_by_coordinate(self, coor="z", coor_min=0.0, coor_max=np.inf):
        """
        This function will remove the spots that are outside coordinate limits

        Parameters
        ----------
        coor : string, optional
            which coordinate to process ('x','y' or 'z'). The default is 'z'.
        coor_min : float, optional
            minimum value. The default is 0..
        coor_max : float, optional
            maximum value. The default is np.inf.

        Returns
        -------
        updated trace table is kept in self.data

        """
        trace_table = self.data

        if len(trace_table) > 0:
            # indexes trace file
            trace_table_indexed = trace_table.group_by("Trace_ID")

            # iterates over traces
            print(f"\n$ Will keep localizations with {coor_min} < {coor} < {coor_max}.")
            print(
                f"$ Number of original spots / traces: {len(trace_table)} / {len(trace_table_indexed.groups)}"
            )
            rows_to_remove = []
            for idx, row in enumerate(trace_table):
                coordinate = float(row[coor])

                if coordinate < coor_min or coordinate > coor_max:
                    rows_to_remove.append(idx)
                    # coordinates.append(coordinate)

            print(f"$ Number of spots to remove: {len(rows_to_remove)}")

            trace_table.remove_rows(rows_to_remove)

            if len(trace_table) > 0:
                trace_table_indexed = trace_table.group_by("Trace_ID")
                number_traces_left = len(trace_table_indexed.groups)
            else:
                number_traces_left = 0

            print(
                f"$ Number of spots / traces left: {len(trace_table)} / {number_traces_left}"
            )

        else:
            print("! Error: you are trying to filter an empty trace table!")
        self.data = trace_table

    def filter_by_intensity(self, trace, localizations, intensity_min):
        """
        Filters localizations in the trace file based on intensity from the localization table.
        """
        localizations.add_index("Buid")  # Add an index for fast lookup

        rows_to_remove = []
        number_spots = len(trace.data)
        intensities_kept = list()
        for idx, row in enumerate(trace.data):
            spot_id = row["Spot_ID"]
            try:
                intensity = localizations.loc[spot_id]["peak"]
                if intensity < intensity_min:
                    rows_to_remove.append(idx)
                else:
                    intensities_kept.append(intensity)
            except KeyError:
                continue  # If Spot_ID is not found, keep the entry

        trace.data.remove_rows(rows_to_remove)
        print(
            f"> Removed {len(rows_to_remove)}/{number_spots} localizations below intensity threshold ({intensity_min})."
        )
        print(f"> Number of rows in filtered trace table: {len(trace.data)}")

        return intensities_kept

    def barcode_statistics(self, trace_table):
        """
        calculates the number of times a barcode is repeated in a trace for all traces in a trace table

        Parameters
        ----------
        trace_table : ASTROPY table
            trace table.

        Returns
        -------
        collective_barcode_stats : dict
            dict with barcode identities as keys and a list of the number of times it was present in each trace treated.

        """
        collective_barcode_stats = {}

        trace_table_indexed = trace_table.group_by("Trace_ID")

        # iterates over traces
        print("$ Calculating barcode stats...")

        for trace in tqdm(trace_table_indexed.groups):
            unique_barcodes = list(set(trace["Barcode #"].data))
            barcodes = list(trace["Barcode #"].data)

            # if number_unique_barcodes < number_barcodes:
            barcode_stats = {}
            for barcode in unique_barcodes:
                barcode_rep = barcodes.count(barcode)
                barcode_stats[str(barcode)] = barcode_rep

                if str(barcode) in collective_barcode_stats:
                    collective_barcode_stats[str(barcode)].append(barcode_rep)
                else:
                    collective_barcode_stats[str(barcode)] = [barcode_rep]

        return collective_barcode_stats

    def plots_barcode_statistics(
        self,
        collective_barcode_stats,
        file_name="barcode_stats",
        kind="violin",
        norm=True,
        format="png",
    ):
        """
        plots the collecive_bracode stats (see previous function)

        Parameters
        ----------
        collective_barcode_stats : dict
            dict with barcode identities as keys and a list of the number of times it was present in each trace treated.
        file_name : str, optional
            output filename for saving figure. The default is 'barcode_stats.png'.
        kind : str, optional
            Options for plotting styles: 'violin' or 'matrix'. The default is 'violin'.

        Returns
        -------
        None.

        """
        sorted_barcodes = sorted([int(x) for x in collective_barcode_stats.keys()])
        data = [collective_barcode_stats[str(key)] for key in sorted_barcodes]

        fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(15, 15))

        label, density = ("frequency", True) if norm else ("counts", False)
        ax1.set_title("Relative barcode frequencies", fontsize=30)

        if "violin" in kind:
            self._extracted_from_plots_barcode_statistics_38(ax1, data, sorted_barcodes)
        else:
            bins = range(1, 10)
            matrix = np.zeros((len(sorted_barcodes), len(bins) - 1))
            for idx, barcode_data in enumerate(data):
                matrix[idx, :], _ = np.histogram(
                    barcode_data, bins=bins, density=density
                )
            bin_number = list(bins)
            pos = ax1.imshow(np.transpose(matrix), cmap="Reds")
            ax1.set_xticks(np.arange(matrix.shape[0]), sorted_barcodes, fontsize=15)
            ax1.set_yticks(np.arange(0, len(bins)), bin_number, fontsize=15)
            ax1.set_ylabel("number of barcodes", fontsize=25)
            ax1.set_xlabel("barcode IDs", fontsize=25)
            fig.colorbar(
                pos, ax=ax1, location="bottom", anchor=(0.5, 1), shrink=0.4, label=label
            )
        print(
            f"$ Exporting relative barcode frequencies figure to: {file_name}.{format}"
        )
        fig.savefig(f"{file_name}.{format}")

    # TODO Rename this here and in `plots_barcode_statistics`
    def _extracted_from_plots_barcode_statistics_38(self, ax1, data, sorted_barcodes):
        ax1.set_ylabel("number of barcodes")
        ax1.violinplot(data)

        ax1.set_xticks(np.arange(1, len(sorted_barcodes) + 1), labels=sorted_barcodes)
        ax1.set_xlim(0.25, len(sorted_barcodes) + 0.75)
        ax1.set_ylim(0.0, 10)
        ax1.set_xlabel("barcode id")

    def trace_remove_label(self, label=""):
        """
        This function will remove traces that do not contain the word 'label' in the 'label' column

        Parameters
        ----------
        label : TYPE, string
            the labe to keep. The default is "".

        Returns
        -------
        None.

        """
        trace_table = self.data

        trace_table_new = trace_table.copy()

        rows_to_remove = []

        for idx, row in enumerate(trace_table_new):
            if label in row["label"]:
                rows_to_remove.append(idx)

        trace_table_new.remove_rows(rows_to_remove)

        removed = len(trace_table) - len(trace_table_new)
        print(f"$ Removed {removed} spots that contained the label: {label}")

        self.data = trace_table_new

    def trace_keep_label(self, label=""):
        """
        This function will remove traces that do not contain the word 'label' in the 'label' column

        Parameters
        ----------
        label : TYPE, string
            the labe to keep. The default is "".

        Returns
        -------
        None.

        """
        trace_table = self.data

        trace_table_new = trace_table.copy()

        rows_to_remove = []

        for idx, row in enumerate(trace_table_new):
            if label not in row["label"]:
                rows_to_remove.append(idx)

        trace_table_new.remove_rows(rows_to_remove)

        removed = len(trace_table) - len(trace_table_new)
        print(f"$ Removed {removed} spots that did not contain the label: {label}")

        self.data = trace_table_new

    def filter_repeated_barcodes(self, trace_file="mock"):
        """
        This function will remove the barcodes that are present more than once in a trace.
        All other barcodes are kept.

        Parameters
        ----------


        Returns
        -------
        updated trace table is kept in self.data

        """
        trace_table = self.data
        trace_table_new = trace_table.copy()
        print("\n$ Removing spots with repeated barcodes...")
        if len(trace_table) > 0:
            # indexes trace file
            trace_table_indexed = trace_table.group_by("Trace_ID")

            # iterates over traces
            print(
                f"\n$ Number of original \n spots: {len(trace_table)} \n traces: {len(trace_table_indexed.groups)}"
            )

            # calculates the statistics for the table before processing
            collective_barcode_stats = self.barcode_statistics(trace_table)

            # plots statistics of barcodes and saves in file
            self.plots_barcode_statistics(
                collective_barcode_stats,
                file_name=f"{trace_file.split('.')[0]}_before_filtering",
                kind="matrix",
                norm=True,
            )

            # iterates over traces
            spots_to_remove = []
            for trace in tqdm(trace_table_indexed.groups):
                unique_barcodes = list(set(trace["Barcode #"].data))
                number_unique_barcodes = len(unique_barcodes)
                barcodes = list(trace["Barcode #"].data)
                number_barcodes = len(barcodes)

                if number_unique_barcodes < number_barcodes:
                    trace_indexed_by_barcode = trace.group_by("Barcode #")

                    for row in trace_indexed_by_barcode:
                        barcode = row["Barcode #"].data
                        barcode_rep = barcodes.count(barcode)

                        # if a barcode is more than once I will remove both instances
                        if barcode_rep > 1:
                            spots_to_remove.append(row["Spot_ID"])

            print(f"$ Number of spots to remove: {len(spots_to_remove)}")
            print("$ Removing repeated spots...")

            rows_to_remove = []
            for idx, row in enumerate(trace_table):
                spot_id = row["Spot_ID"]

                if spot_id in spots_to_remove:
                    rows_to_remove.append(idx)

            trace_table_new.remove_rows(rows_to_remove)

            print(f"$ Number of rows to remove: {len(rows_to_remove)}")

            if len(trace_table_new) > 0:
                trace_table_indexed = trace_table_new.group_by("Trace_ID")

            print(
                f"$ After filtering, I see \n spots: {len(trace_table_new)} \n traces: {len(trace_table_indexed.groups)}"
            )

            # calculates the statistics for the table before processing
            collective_barcode_stats_new = self.barcode_statistics(trace_table_new)

            # plots statistics of barcodes and saves in file
            self.plots_barcode_statistics(
                collective_barcode_stats_new,
                file_name=f"{trace_file.split('.')[0]}_filtered",
                kind="matrix",
                norm=False,
            )

        else:
            print("! Error: you are trying to filter an empty trace table!")
        self.data = trace_table_new

    def remove_duplicates_loc(self, localization_table=None):
        """
        Removes duplicated barcodes within each trace.
        If a localization_table is provided, keeps only the spot with the highest intensity ("peak").
        Otherwise, removes all instances of duplicated barcodes.

        Parameters
        ----------
        localization_table : astropy Table, optional
            Localization table with 'Buid' and 'peak' columns. Used to select spot with highest intensity.

        Returns
        -------
        Updates self.data with filtered trace table.
        """
        trace_table = self.data
        trace_table_new = trace_table.copy()
        print("\n$ Removing duplicated barcodes within traces...")

        if len(trace_table) == 0:
            print("! Error: you are trying to filter an empty trace table!")
            return

        trace_table_indexed = trace_table.group_by("Trace_ID")
        trace_table.add_index("Spot_ID")  # Add index for faster lookup
        rows_to_remove = []

        if localization_table is not None:
            print("$ Using intensity to resolve duplicates...")
            localization_table.add_index("Buid")

            for trace in tqdm(trace_table_indexed.groups):
                barcode_groups = trace.group_by("Barcode #").groups
                for group in barcode_groups:
                    if len(group) == 1:
                        continue  # no duplicates

                    peaks = []
                    for row in group:
                        spot_id = row["Spot_ID"]
                        try:
                            peak = localization_table.loc[spot_id]["peak"]
                        except KeyError:
                            peak = -1
                        peaks.append(peak)

                    max_idx = peaks.index(max(peaks))
                    for idx, row in enumerate(group):
                        if idx != max_idx:
                            global_idx = trace_table.loc_indices[row["Spot_ID"]]
                            rows_to_remove.append(global_idx)

        else:
            print(
                "$ No localization table provided. Removing all instances of duplicated barcodes."
            )

            for trace in trace_table_indexed.groups:
                barcode_groups = trace.group_by("Barcode #").groups
                for group in barcode_groups:
                    if len(group) <= 1:
                        continue
                    for row in group:
                        global_idx = trace_table.loc_indices[row["Spot_ID"]]
                        rows_to_remove.append(global_idx)

        trace_table_new.remove_rows(rows_to_remove)

        print(f"$ Number of rows to remove: {len(rows_to_remove)}")

        if len(trace_table_new) > 0:
            number_traces_left = len(trace_table_new.group_by("Trace_ID").groups)
        else:
            number_traces_left = 0

        print(
            f"$ After filtering, I see \n spots: {len(trace_table_new)} \n traces: {number_traces_left}"
        )

        self.data = trace_table_new

    def remove_duplicates(self):
        """
        removes duplicated (identical) spots

        Parameters
        ----------


        Returns
        -------
        trace_table : ASTROPY Table
            output trace table.
        """
        trace_table = self.data
        trace_table_new = trace_table.copy()
        print("\n$ Removing duplicated barcodes within traces...")
        if len(trace_table) > 0:
            # indexes trace file
            trace_table_indexed = trace_table.group_by("Spot_ID")

            # finds barcodes with the same UID and stores UIDs in list
            spots_to_remove = [
                trace["Spot_ID"][0]
                for trace in tqdm(trace_table_indexed.groups)
                if len(trace) > 1
            ]

            # finds row of the first offending barcode
            # this only removes one of the duplicated barcodes --> assumes at most there are two copies
            rows_to_remove = []
            for idx, row in enumerate(trace_table):
                spot_id = row["Spot_ID"]
                if spot_id in spots_to_remove:
                    rows_to_remove.append(idx)

            # removes from table
            trace_table_new.remove_rows(rows_to_remove)

            print(f"$ Number of rows to remove: {len(rows_to_remove)}")

            if len(trace_table_new) > 0:
                trace_table_indexed = trace_table_new.group_by("Trace_ID")
                number_traces_left = len(trace_table_indexed.groups)
            else:
                number_traces_left = 0

            print(
                f"$ After filtering, I see \n spots: {len(trace_table_new)} \n traces: {number_traces_left}"
            )

        else:
            print("! Error: you are trying to filter an empty trace table!")

        self.data = trace_table_new

    def remove_barcode(self, remove_barcode=None):
        """
        Removes a specific barcode from a trace table

        Returns
        -------
        trace_table : ASTROPY Table
            output trace table.
        """

        if remove_barcode is not None:
            print(f"$ Removing barcode <{remove_barcode}>")

            trace_table = self.data
            trace_table_new = trace_table.copy()

            # indexes trace file
            trace_table_indexed = trace_table.group_by("Barcode #")
            number_barcodes_before = len(trace_table_indexed.groups)

            # iterates over traces
            spots_to_remove = []
            for sub_table_barcode in trace_table_indexed.groups:
                barcode_name = list(set(sub_table_barcode["Barcode #"]))
                if int(remove_barcode) in barcode_name:
                    print(f"$ Found barcode: {barcode_name}")
                    spots_to_remove.extend(row["Spot_ID"] for row in sub_table_barcode)
            print(f"$ Number of spots to remove: {len(spots_to_remove)}")

            # builds the list with the rows to remove
            rows_to_remove = []
            for idx, row in enumerate(tqdm(trace_table)):
                spot_id = row["Spot_ID"]

                if spot_id in spots_to_remove:
                    rows_to_remove.append(idx)

            # removes targeted spots
            trace_table_new.remove_rows(rows_to_remove)

            # provides statistics
            trace_table_indexed_new = trace_table_new.group_by("Barcode #")
            number_barcodes_left = len(trace_table_indexed_new.groups)
            print(
                f"\n$ Number of barcodes \n\t original: {number_barcodes_before} \n\t after: {number_barcodes_left}"
            )

        self.data = trace_table_new

    def filter_traces_by_n(self, minimum_number_barcodes=2):
        """
        Removes rows in trace table with less than `minimum_number_barcodes` barcodes

        Parameters
        ----------
        trace_table : ASTROPY Table
            input trace table.
        minimum_number_barcodes : TYPE, optional
            minimum number of barcodes in trace. The default is 1.

        Returns
        -------
        trace_table : ASTROPY Table
            output trace table.

        """

        trace_table = self.data

        # indexes trace file
        trace_table_indexed = trace_table.group_by("Trace_ID")

        # iterates over traces
        print(f"\n$ Removing traces with < {minimum_number_barcodes} spots")
        print(
            f"$ Number of original spots / traces: {len(trace_table)} / {len(trace_table_indexed.groups)}"
        )

        barcodes_to_remove = []
        print("$ Analyzing traces...")

        for trace in tqdm(trace_table_indexed.groups):
            number_unique_barcodes = len(list(set(trace["Barcode #"].data)))

            if number_unique_barcodes < minimum_number_barcodes:
                barcodes_to_remove.append(list(trace["Spot_ID"].data))

        print(f"$ Number of traces to remove: {len(barcodes_to_remove)}")

        list_barcode_to_remove = []
        for barcodes in tqdm(barcodes_to_remove):
            list_barcode_to_remove.extend(iter(barcodes))
        rows_to_remove = []

        print("$ Finding which rows to remove...")
        for idx, row in enumerate(tqdm(trace_table)):
            spot_id = row["Spot_ID"]
            if spot_id in list_barcode_to_remove:
                rows_to_remove.append(idx)

        trace_table.remove_rows(rows_to_remove)
        if len(trace_table) > 0:
            trace_table_indexed = trace_table.group_by("Trace_ID")
            number_traces_left = len(trace_table_indexed.groups)
        else:
            number_traces_left = 0

        print(
            f"$ Number of spots / traces left: {len(trace_table)} / {number_traces_left}"
        )

        self.data = trace_table
