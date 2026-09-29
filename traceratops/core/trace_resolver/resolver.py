"""Public interface for barcode-aware polymer assignment."""

import numpy as np

from .beam import beam_search
from .model import EmpiricalDistanceModel
from .results import CandidateScore, ResolutionResult
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

    def resolve_one_polymer_candidates(self, trace, minimum_confidence=None):
        """Resolve duplicates using a fresh constrained search per candidate.

        The ordinary one-polymer result supplies the unchanged trace-level score
        diagnostics. Each duplicate decision is instead based on the best full
        trace found while forcing that localization. Other duplicate barcodes
        remain unconstrained and therefore optimize jointly in every run. A
        final search forces all confident winners together and excludes every
        ambiguous barcode from the emitted assignment.
        """
        threshold = (
            self.minimum_confidence
            if minimum_confidence is None
            else minimum_confidence
        )
        if threshold < 0:
            raise ValueError("minimum_confidence must be non-negative")
        global_result = self.resolve(trace, n_polymers=1)
        positions = self.distance_model.positions(trace)
        barcodes = np.asarray(trace["Barcode #"])
        comparisons = []
        ambiguous = []
        confident_winners = {}

        for barcode in sorted(
            set(barcodes), key=lambda value: np.min(positions[barcodes == value])
        ):
            candidates = np.flatnonzero(barcodes == barcode).tolist()
            if len(candidates) == 1:
                continue
            scored = []
            for index in candidates:
                states = beam_search(
                    trace,
                    positions,
                    self.scorer,
                    1,
                    self.beam_width,
                    self.rejection_cost,
                    self.minimum_polymer_size,
                    forced_candidates={barcode: index},
                )
                scored.append((float(states[0].score), index))
            scored.sort(key=lambda item: (item[0], item[1]))
            best_score, best_index = scored[0]
            second_score = scored[1][0]
            gap = second_score - best_score
            confidence = gap / (abs(best_score) + abs(second_score) + 1e-12)
            confident = confidence >= threshold
            if confident:
                confident_winners[barcode] = best_index
            else:
                ambiguous.append(barcode)

            ranks = {index: rank for rank, (_, index) in enumerate(scored, 1)}
            for score, index in scored:
                selected = confident and index == best_index
                decision = (
                    "ambiguous_remove_all"
                    if not confident
                    else "keep" if selected else "reject"
                )
                comparisons.append(
                    CandidateScore(
                        barcode=barcode,
                        candidate_index=index,
                        spot_id=trace["Spot_ID"][index],
                        score=score,
                        rank=ranks[index],
                        best_score=best_score,
                        second_score=second_score,
                        raw_score_gap=gap,
                        confidence=confidence,
                        selected=selected,
                        ambiguity_threshold=threshold,
                        decision=decision,
                    )
                )

        # Candidate comparisons above deliberately leave other duplicates free.
        # Finalize with one joint run so the emitted assignment is a beam-search
        # solution compatible with every confident winner at once. Ambiguous
        # barcodes are skipped entirely rather than being reintroduced by the
        # ordinary one-polymer rule that selects one candidate per barcode.
        joint_states = beam_search(
            trace,
            positions,
            self.scorer,
            1,
            self.beam_width,
            self.rejection_cost,
            self.minimum_polymer_size,
            forced_candidates=confident_winners,
            excluded_barcodes=ambiguous,
        )
        result = ResolutionResult(
            assignments=joint_states[0].assignments.copy(),
            status="cleaned",
            inferred_n_polymers=1,
            best_score=global_result.best_score,
            alternative_score=global_result.alternative_score,
            raw_score_gap=global_result.raw_score_gap,
            confidence_score=global_result.confidence_score,
            n_ambiguous_barcodes=len(ambiguous),
            ambiguous_barcodes=tuple(ambiguous),
        )
        return result, tuple(comparisons)
