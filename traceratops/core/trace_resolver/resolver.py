"""Public interface for barcode-aware polymer assignment."""

import numpy as np

from .beam import beam_search
from .model import EmpiricalDistanceModel
from .results import ResolutionResult
from .scoring import PolymerScorer


class TraceResolver:
    """Resolve one or two polymers while allowing detections to be rejected."""

    def __init__(
        self,
        distance_model=None,
        history_mode="multi",
        distance_score="residual",
        history_length=3,
        beam_width=100,
        rejection_cost=16.0,
        minimum_polymer_size=2,
        minimum_confidence=0.05,
    ):
        self.distance_model = distance_model or EmpiricalDistanceModel()
        self.scorer = PolymerScorer(
            self.distance_model, history_mode, distance_score, history_length
        )
        self.beam_width = beam_width
        self.rejection_cost = rejection_cost
        self.minimum_polymer_size = minimum_polymer_size
        self.minimum_confidence = minimum_confidence

    def resolve(self, trace, n_polymers):
        if n_polymers not in {1, 2}:
            raise ValueError("n_polymers must be 1 or 2")
        positions = self.distance_model.positions(trace)
        states = beam_search(
            trace,
            positions,
            self.scorer,
            n_polymers,
            self.beam_width,
            self.rejection_cost,
            self.minimum_polymer_size,
        )
        best = states[0]
        alternative = states[1] if len(states) > 1 else None
        alternative_score = alternative.score if alternative else None
        gap = alternative.score - best.score if alternative else None
        confidence = (
            gap / (abs(best.score) + abs(alternative.score) + 1e-12)
            if alternative
            else None
        )

        assignments = best.assignments.copy()
        ambiguous_barcodes = ()
        if n_polymers == 2:
            has_two_polymers = all(
                count >= self.minimum_polymer_size for count in best.counts
            )
            is_confident = (
                has_two_polymers
                and confidence is not None
                and confidence >= self.minimum_confidence
            )
            status = "split" if is_confident else "removed_ambiguous"
            inferred = 2 if is_confident else 0
            if not is_confident:
                assignments[:] = -1
        else:
            if confidence is not None and confidence < self.minimum_confidence:
                changed = np.flatnonzero(assignments != alternative.assignments)
                barcodes = np.asarray(trace["Barcode #"])
                unique, counts = np.unique(barcodes, return_counts=True)
                repeated = set(unique[counts > 1].tolist())
                ambiguous_barcodes = tuple(
                    barcode
                    for barcode in np.unique(barcodes[changed]).tolist()
                    if barcode in repeated
                )
                for barcode in ambiguous_barcodes:
                    assignments[barcodes == barcode] = -1
            status = "cleaned"
            inferred = 1

        return ResolutionResult(
            assignments=assignments,
            status=status,
            inferred_n_polymers=inferred,
            best_score=best.score,
            alternative_score=alternative_score,
            raw_score_gap=gap,
            confidence_score=confidence,
            n_ambiguous_barcodes=len(ambiguous_barcodes),
            ambiguous_barcodes=ambiguous_barcodes,
        )
