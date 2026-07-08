#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This class:
    - iterates over chromatin traces
        - calculates the pair-wise distances for each single-cell mask
        - outputs are:
            - Table with #cell #PWD #coordinates (e.g. buildsPWDmatrix_3D_order:0_ROI:1.ecsv)
            - NPY array with single cell PWD single cell matrices (e.g. buildsPWDmatrix_3D_HiMscMatrix.npy)
            - NPY array with barcode identities (e.g. buildsPWDmatrix_3D_uniqueBarcodes.ecsv)
            - the files with no "3D" tag contain data analyzed using 2D localizations.

    - Single-cell results are combined together to calculate:
        - Distribution of pairwise distance for each barcode combination
        - Ensemble mean pairwise distance matrix using mean of distribution
        - Ensemble mean pairwise distance matrix using Kernel density estimation
        - Ensemble Hi-M matrix using a predefined threshold
        - For each of these files, there is an image in PNG format saved. Images containing "3D" are for 3D other are for 2D.
"""

import glob
import os
from pathlib import Path

import numpy as np
from astropy.table import unique
from joblib import Parallel, delayed
from sklearn.metrics import pairwise_distances
from tqdm import tqdm

from traceratops.core.chromatin_trace_table import ChromatinTraceTable
from traceratops.core.him_matrix_operations import (
    calculate_contact_probability_matrix,
    plot_distance_histograms,
    plot_matrix,
)


class BuildMatrix:
    def __init__(
        self,
        param,
        acq_params={"zBinning": None, "pixelSizeXY": None, "pixelSizeZ": None},
        colormaps=dict(),
    ):
        self.current_param = param
        self.colormaps = (
            colormaps
            if colormaps
            else {
                "PWD_KDE": "terrain",
                "PWD_median": "terrain",
                "contact": "coolwarm",
                "Nmatrix": "Blues",
            }
        )

        self.initialize_parameters(acq_params)

        # initialize with default values
        self.current_folder = []

    def initialize_parameters(self, acq_params: dict):
        # initializes parameters from current_param
        # TODO: Check this condition
        if not isinstance(self.current_param, dict):
            self.z_binning = acq_params["zBinning"]
            self.pixel_size_xy = acq_params["pixelSizeXY"]
            self.pixel_size_z_0 = acq_params["pixelSizeZ"]
            self.pixel_size_z = self.z_binning * self.pixel_size_z_0
            self.pixel_size = [
                self.pixel_size_xy,
                self.pixel_size_xy,
                self.pixel_size_z,
            ]
            self.log_name_md = self.current_param.param_dict["fileNameMD"]
        else:
            self.log_name_md = "trace_to_matrix.md"

    def calculate_pwd_single_mask(self, x, y, z):
        """
        Calculates PWD between barcodes detected in a given mask. For this:
            - converts xyz pixel coordinates into nm using self.pixel_size dictionary
            - calculates pair-wise distance matrix in nm
            - converts it into pixel units using self.pixel_size['x'] as an isotropic pixelsize.

        Parameters
        ----------
        r1: list of floats with xyz coordinates for spot 1 in microns
        r2: list of floats with xyz coordinates for spot 2 in microns

        Returns
        -------
        Returns pairwise distance matrix between barcodes in microns

        """
        r_mum = np.column_stack((x, y, z))
        return pairwise_distances(r_mum)

    @staticmethod
    def _build_trace_distance_slice(
        trace, unique_barcode_to_index, number_unique_barcodes, mode, distance_threshold
    ):
        """Build one single-trace distance matrix slice.

        The default ``min`` mode is vectorized because it is the runtime-critical
        path used by ``trace_to_matrix``. The ``last`` and legacy iterative
        ``mean`` modes keep their original row-major update order so duplicate
        barcode observations produce byte-for-byte equivalent results.
        """
        barcodes_to_process = np.asarray(trace["Barcode #"].data)
        coordinates = np.column_stack(
            (
                np.asarray(trace["x"].data),
                np.asarray(trace["y"].data),
                np.asarray(trace["z"].data),
            )
        )
        pwd_matrix = pairwise_distances(coordinates)
        barcode_indices = np.fromiter(
            (unique_barcode_to_index[barcode] for barcode in barcodes_to_process),
            dtype=np.intp,
            count=len(barcodes_to_process),
        )

        matrix_slice = np.full(
            (number_unique_barcodes, number_unique_barcodes), np.nan, dtype=float
        )

        if mode == "min":
            different_barcodes = (
                barcodes_to_process[:, None] != barcodes_to_process[None, :]
            )
            valid_distances = different_barcodes & (pwd_matrix < distance_threshold)
            if np.any(valid_distances):
                row_positions, column_positions = np.nonzero(valid_distances)
                rows = barcode_indices[row_positions]
                columns = barcode_indices[column_positions]
                distances = pwd_matrix[row_positions, column_positions]
                flat_matrix = np.full(matrix_slice.size, np.inf, dtype=float)
                np.minimum.at(
                    flat_matrix,
                    np.ravel_multi_index((rows, columns), matrix_slice.shape),
                    distances,
                )
                matrix_slice = flat_matrix.reshape(matrix_slice.shape)
                matrix_slice[np.isinf(matrix_slice)] = np.nan
            return matrix_slice

        for barcode1, ibarcode1 in zip(
            barcodes_to_process, range(len(barcodes_to_process))
        ):
            index_barcode_1 = barcode_indices[ibarcode1]
            for barcode2, ibarcode2 in zip(
                barcodes_to_process, range(len(barcodes_to_process))
            ):
                if barcode1 != barcode2:
                    index_barcode_2 = barcode_indices[ibarcode2]
                    newdistance = pwd_matrix[ibarcode1, ibarcode2]
                    if newdistance < distance_threshold:
                        if mode == "last":
                            matrix_slice[index_barcode_1][index_barcode_2] = newdistance
                        elif mode == "mean":
                            matrix_slice[index_barcode_1][index_barcode_2] = np.nanmean(
                                [
                                    newdistance,
                                    matrix_slice[index_barcode_1][index_barcode_2],
                                ]
                            )
                        else:
                            raise ValueError(
                                f"Unsupported distance matrix mode: {mode}"
                            )

        return matrix_slice

    def build_distance_matrix(self, mode="min", distance_threshold=np.inf, n_jobs=1):
        """
        Builds pairwise distance matrix from a coordinates table

        Parameters
        ----------
        mode : string, optional
            The default is "mean": calculates the mean distance if there are several combinations possible.
            "min": calculates the minimum distance if there are several combinations possible.
            "last": keeps the last distance calculated
        n_jobs : int, optional
            Number of parallel workers used across independent traces. ``1`` keeps
            the serial execution path; ``-1`` uses all available workers.

        Returns
        -------
        self.sc_matrix the single-cell PWD matrix
        self.unique_barcodes list of unique barcodes

        """
        # detects number of unique traces from trace table
        number_matrices = len(unique(self.trace_table.data, keys="Trace_ID"))

        # finds unique barcodes from trace table
        unique_barcodes = unique(self.trace_table.data, keys="Barcode #")[
            "Barcode #"
        ].data
        number_unique_barcodes = unique_barcodes.shape[0]

        print(
            f"$ Found {number_unique_barcodes} barcodes and {number_matrices} traces.",
            "INFO",
        )

        unique_barcode_to_index = {
            barcode: index for index, barcode in enumerate(unique_barcodes)
        }

        # loops over traces
        print("> Processing traces...", "INFO")
        data_traces = self.trace_table.data.group_by("Trace_ID")
        trace_groups = list(data_traces.groups)

        if n_jobs == 1:
            slices = [
                self._build_trace_distance_slice(
                    trace,
                    unique_barcode_to_index,
                    number_unique_barcodes,
                    mode,
                    distance_threshold,
                )
                for trace in tqdm(trace_groups, total=number_matrices)
            ]
        else:
            slices = Parallel(n_jobs=n_jobs)(
                delayed(self._build_trace_distance_slice)(
                    trace,
                    unique_barcode_to_index,
                    number_unique_barcodes,
                    mode,
                    distance_threshold,
                )
                for trace in tqdm(trace_groups, total=number_matrices)
            )

        self.sc_matrix = np.stack(slices, axis=2)
        self.unique_barcodes = unique_barcodes

    def calculate_n_matrix(self):
        number_cells = self.sc_matrix.shape[2]

        if number_cells > 0:
            n_matrix = np.sum(~np.isnan(self.sc_matrix), axis=2)
        else:
            number_barcodes = self.sc_matrix.shape[0]
            n_matrix = np.zeros((number_barcodes, number_barcodes))

        self.n_matrix = n_matrix

    def get_output_prefix(self, file, outputFolder=None):
        """
        Build the output filename prefix for a trace file.

        Outputs are written next to the input trace by default, or inside
        outputFolder when it is provided. In both cases, the input trace
        filename stem is preserved in the generated filenames.
        """
        trace_path = Path(file)
        output_dir = (
            Path(outputFolder) if outputFolder is not None else trace_path.parent
        )
        output_dir.mkdir(parents=True, exist_ok=True)
        return str(output_dir / f"{trace_path.stem}_Matrix")

    def plots_all_matrices(self, output_filename):
        """
        Plots all matrices after analysis

        Parameters
        ----------
        output_filename : str
            Output filename prefix used for generated plots.

        Returns
        -------
        None.

        """
        number_rois = 1  # by default we plot one ROI at a time.

        clim_scale = 1.0  # factor to multiply the clim by. If 1, the clim will be the mean of the PWD distribution of the whole map
        pixel_size = 1  # this is 1 as coordinates are in microns.
        n_cells = self.sc_matrix.shape[2]

        # plots PWD matrix
        # uses KDE
        plot_matrix(
            self.sc_matrix,
            self.unique_barcodes,
            pixel_size,
            number_rois,
            output_filename,
            self.log_name_md,
            figtitle="PWD matrix - KDE",
            mode="KDE",  # median or KDE
            clim=clim_scale * np.nanmean(self.sc_matrix),
            n_cells=n_cells,
            c_m=self.colormaps["PWD_KDE"],
            cmtitle="distance, um",
            filename_ending="_PWDmatrixKDE.png",
        )

        # uses median
        plot_matrix(
            self.sc_matrix,
            self.unique_barcodes,
            pixel_size,
            number_rois,
            output_filename,
            self.log_name_md,
            figtitle="PWD matrix - median",
            mode="median",  # median or KDE
            clim=clim_scale * np.nanmean(self.sc_matrix),
            cmtitle="distance, um",
            n_cells=n_cells,
            c_m=self.colormaps["PWD_median"],
            filename_ending="_PWDmatrixMedian.png",
        )

        # calculates and plots contact probability matrix from merged samples/datasets
        him_matrix = calculate_contact_probability_matrix(
            self.sc_matrix,
            pixel_size,
            remove_nan=True,
        )
        n_cells = self.sc_matrix.shape[2]
        c_scale = him_matrix.max()
        plot_matrix(
            him_matrix,
            self.unique_barcodes,
            pixel_size,
            number_rois,
            output_filename,
            self.log_name_md,
            figtitle="Hi-M matrix",
            mode="counts",
            clim=c_scale,
            n_cells=n_cells,
            c_m=self.colormaps["contact"],
            cmtitle="proximity frequency",
            filename_ending="_HiMmatrix.png",
        )

        # plots n_matrix
        plot_matrix(
            self.n_matrix,
            self.unique_barcodes,
            pixel_size,
            number_rois,
            output_filename,
            self.log_name_md,
            figtitle="N-matrix",
            mode="counts",
            n_cells=n_cells,
            clim=np.max(self.n_matrix),
            c_m=self.colormaps["Nmatrix"],
            cmtitle="number of measurements",
            filename_ending="_Nmatrix.png",
        )

        plot_distance_histograms(
            self.sc_matrix,
            pixel_size,
            output_filename,
            self.log_name_md,
            mode="KDE",
            kernel_width=0.25,
            optimize_kernel_width=False,
        )

    def save_matrices(self, output_filename):
        # saves output
        np.save(f"{output_filename}_PWDscMatrix.npy", self.sc_matrix)
        print(f"$ saved: {output_filename}_PWDscMatrix.npy")

        np.savetxt(
            f"{output_filename}_uniqueBarcodes.ecsv",
            self.unique_barcodes,
            delimiter=" ",
            fmt="%d",
        )

        print(f"$ saved: {output_filename}_uniqueBarcodes.ecsv")

        np.save(f"{output_filename}_Nmatrix.npy", self.n_matrix)
        print(f"$ saved: {output_filename}_Nmatrix.npy")

    def launch_analysis(
        self, file, distance_threshold=np.inf, outputFolder=None, n_jobs=1
    ):
        """
        run analysis for a chromatin trace table.

        Returns
        -------
        None.

        """

        # creates and loads trace table
        self.trace_table = ChromatinTraceTable()
        self.trace_table.load(file)
        output_filename = self.get_output_prefix(file, outputFolder)

        # runs calculation of PWD matrix
        self.build_distance_matrix(
            "min", distance_threshold=distance_threshold, n_jobs=n_jobs
        )  # mean min last

        # calculates N-matrix: number of PWD distances for each barcode combination
        self.calculate_n_matrix()

        # runs plotting operations
        self.plots_all_matrices(output_filename)

        # saves matrix
        self.save_matrices(output_filename)

    def run(self, data_path, matrix_params):
        self.label = "barcode"
        self.current_folder = data_path

        # reads chromatin traces
        files = [
            x
            for x in glob.glob(
                data_path
                + os.sep
                + matrix_params.folder
                + os.sep
                + "data"
                + os.sep
                + "Trace_*.ecsv"
            )
            if "uniqueBarcodes" not in x
        ]

        if not files:
            print("$ No chromatin trace table found !", "WARN")
            return

        print(f"> Will process {len(files)} trace tables with names:")
        for file in files:
            print(f"\t{os.path.basename(file)}")

        for file in files:
            self.launch_analysis(file)

        print(
            f"$ {len(files)} chromatin trace tables processed in {self.current_folder}"
        )
