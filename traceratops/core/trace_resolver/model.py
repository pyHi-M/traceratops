"""Empirical physical-distance models conditioned on genomic separation."""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class DistanceStatistics:
    mean: float
    standard_deviation: float
    n_observations: int


class EmpiricalDistanceModel:
    """Estimate distance statistics from traces with unique barcode identities.

    Sparse separation bins borrow the nearest sufficiently populated bin. If no
    bin is sufficiently populated, the pooled distribution is used. A variance
    floor prevents degenerate training data from creating infinite costs.
    """

    def __init__(self, minimum_observations=3, variance_floor=1e-6):
        if minimum_observations < 1:
            raise ValueError("minimum_observations must be positive")
        if variance_floor <= 0:
            raise ValueError("variance_floor must be positive")
        self.minimum_observations = minimum_observations
        self.variance_floor = variance_floor
        self.statistics = {}
        self.pooled = DistanceStatistics(0.0, 1.0, 0)
        self.genomic_source = "barcode"
        self.fallback_source = "clean_traces"

    @staticmethod
    def _uses_genomic_coordinates(table):
        required = {"Chrom", "Chrom_Start", "Chrom_End"}
        if not required.issubset(table.colnames):
            return False
        midpoints = (
            np.asarray(table["Chrom_Start"], dtype=float)
            + np.asarray(table["Chrom_End"], dtype=float)
        ) / 2
        n_barcodes = len(np.unique(np.asarray(table["Barcode #"])))
        return (
            len(table) > 0
            and len(set(np.asarray(table["Chrom"]).astype(str))) == 1
            and (n_barcodes < 2 or len(np.unique(midpoints)) >= 2)
            and not any(np.ma.getmaskarray(table[name]).any() for name in required)
        )

    def positions(self, table):
        if self._uses_genomic_coordinates(table):
            midpoint = (
                np.asarray(table["Chrom_Start"], dtype=float)
                + np.asarray(table["Chrom_End"], dtype=float)
            ) / 2
            self.genomic_source = "coordinates"
            return midpoint
        self.genomic_source = "barcode"
        return np.asarray(table["Barcode #"], dtype=float)

    def fit(self, table):
        observations = {}
        pooled = []
        if len(table):
            for trace in table.group_by("Trace_ID").groups:
                barcodes = np.asarray(trace["Barcode #"])
                if len(np.unique(barcodes)) != len(barcodes):
                    continue
                positions = self.positions(trace)
                coords = np.column_stack((trace["x"], trace["y"], trace["z"])).astype(
                    float
                )
                for left in range(len(trace)):
                    for right in range(left + 1, len(trace)):
                        separation = float(abs(positions[right] - positions[left]))
                        if separation == 0:
                            continue
                        distance = float(np.linalg.norm(coords[right] - coords[left]))
                        observations.setdefault(separation, []).append(distance)
                        pooled.append(distance)
        if not pooled:
            # Candidate-only files may contain no clean trace. Use the nearest
            # candidate pair between barcode identities as a weak, robust
            # fallback instead of assuming a zero physical scale.
            self.fallback_source = "nearest_candidate_pairs"
            for trace in table.group_by("Trace_ID").groups:
                positions = self.positions(trace)
                coords = np.column_stack((trace["x"], trace["y"], trace["z"])).astype(
                    float
                )
                barcodes = np.asarray(trace["Barcode #"])
                ordered = sorted(
                    set(barcodes),
                    key=lambda value: np.min(positions[barcodes == value]),
                )
                for left_number, left_barcode in enumerate(ordered):
                    for right_barcode in ordered[left_number + 1 :]:
                        left = np.flatnonzero(barcodes == left_barcode)
                        right = np.flatnonzero(barcodes == right_barcode)
                        distances = [
                            float(np.linalg.norm(coords[i] - coords[j]))
                            for i in left
                            for j in right
                        ]
                        if distances:
                            distance = min(distances)
                            separation = float(
                                abs(positions[right[0]] - positions[left[0]])
                            )
                            observations.setdefault(separation, []).append(distance)
                            pooled.append(distance)
        self.statistics = {
            separation: self._summarize(values)
            for separation, values in observations.items()
        }
        self.pooled = (
            self._summarize(pooled) if pooled else DistanceStatistics(0.0, 1.0, 0)
        )
        return self

    def _summarize(self, values):
        values = np.asarray(values, dtype=float)
        standard_deviation = max(float(np.std(values)), np.sqrt(self.variance_floor))
        return DistanceStatistics(
            float(np.mean(values)), standard_deviation, len(values)
        )

    def for_separation(self, separation):
        separation = float(abs(separation))
        populated = {
            delta: stats
            for delta, stats in self.statistics.items()
            if stats.n_observations >= self.minimum_observations
            or self.fallback_source == "nearest_candidate_pairs"
        }
        if separation in populated:
            return populated[separation]
        if populated:
            nearest = min(populated, key=lambda delta: abs(delta - separation))
            return populated[nearest]
        return self.pooled
