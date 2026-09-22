"""Modular barcode-multiplicity classification for trace resolution."""

from dataclasses import dataclass
from enum import Enum

from .results import TraceMultiplicity


class TraceClassification(str, Enum):
    UNCHANGED = "unchanged"
    CLEAN_ONE = "clean_one"
    RESOLVE_TWO = "resolve_two"


@dataclass(frozen=True)
class ClassificationThresholds:
    """Provisional candidate thresholds intended for simulation calibration."""

    min_repeated_barcodes: int = 3
    min_repeated_fraction: float = 0.3
    candidate_rule: str = "both"

    def __post_init__(self):
        if self.min_repeated_barcodes < 1:
            raise ValueError("min_repeated_barcodes must be positive")
        if not 0 <= self.min_repeated_fraction <= 1:
            raise ValueError("min_repeated_fraction must be between 0 and 1")
        if self.candidate_rule not in {"both", "either"}:
            raise ValueError("candidate_rule must be 'both' or 'either'")


def classify_trace(multiplicity, thresholds=None):
    """Classify a trace without embedding a fixed biological interpretation."""
    if not isinstance(multiplicity, TraceMultiplicity):
        multiplicity = TraceMultiplicity.from_barcodes(multiplicity)
    thresholds = thresholds or ClassificationThresholds()
    if multiplicity.n_repeated_barcodes == 0:
        return TraceClassification.UNCHANGED
    tests = (
        multiplicity.n_repeated_barcodes >= thresholds.min_repeated_barcodes,
        multiplicity.fraction_repeated_barcodes >= thresholds.min_repeated_fraction,
    )
    is_two = all(tests) if thresholds.candidate_rule == "both" else any(tests)
    return TraceClassification.RESOLVE_TWO if is_two else TraceClassification.CLEAN_ONE
