"""Barcode-aware reconstruction of chromatin polymers."""

from .classifier import (
    ClassificationThresholds,
    LikelihoodMultiplicityClassifier,
    LikelihoodTraceScores,
    TraceClassification,
    barcode_heterogeneity_pvalue,
    build_multiplicity_matrix,
    classify_trace,
)
from .model import EmpiricalDistanceModel
from .resolver import TraceResolver
from .results import ResolutionResult, TraceMultiplicity

__all__ = [
    "ClassificationThresholds",
    "EmpiricalDistanceModel",
    "LikelihoodMultiplicityClassifier",
    "LikelihoodTraceScores",
    "ResolutionResult",
    "TraceClassification",
    "TraceMultiplicity",
    "TraceResolver",
    "classify_trace",
    "barcode_heterogeneity_pvalue",
    "build_multiplicity_matrix",
]
