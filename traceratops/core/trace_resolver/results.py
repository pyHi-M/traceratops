"""Result and diagnostic objects shared by trace-resolution algorithms."""

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np


@dataclass(frozen=True)
class TraceMultiplicity:
    """Barcode multiplicity measurements for one input trace."""

    n_unique_barcodes: int
    n_repeated_barcodes: int
    fraction_repeated_barcodes: float
    maximum_barcode_multiplicity: int
    n_excess_detections: int

    @classmethod
    def from_barcodes(cls, barcodes):
        _, counts = np.unique(np.asarray(barcodes), return_counts=True)
        n_unique = len(counts)
        n_repeated = int(np.sum(counts > 1))
        return cls(
            n_unique_barcodes=n_unique,
            n_repeated_barcodes=n_repeated,
            fraction_repeated_barcodes=(n_repeated / n_unique if n_unique else 0.0),
            maximum_barcode_multiplicity=int(np.max(counts)) if n_unique else 0,
            n_excess_detections=int(np.sum(counts - 1)),
        )


@dataclass(frozen=True)
class ResolutionResult:
    """Assignments and conservative, non-probabilistic solution diagnostics."""

    assignments: np.ndarray
    status: str
    inferred_n_polymers: int
    best_score: float
    alternative_score: Optional[float]
    raw_score_gap: Optional[float]
    confidence_score: Optional[float]
    n_ambiguous_barcodes: int = 0
    ambiguous_barcodes: Tuple[object, ...] = ()

    @property
    def n_unassigned(self):
        return int(np.sum(self.assignments < 0))

    @property
    def fraction_unassigned(self):
        return (
            self.n_unassigned / len(self.assignments) if len(self.assignments) else 0.0
        )
