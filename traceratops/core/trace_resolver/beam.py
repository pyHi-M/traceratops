"""Beam-search assignment with duplicate and rejection constraints."""

from dataclasses import dataclass
from itertools import product

import numpy as np


@dataclass
class BeamState:
    score: float
    assignments: np.ndarray
    histories: tuple
    counts: tuple


def _canonical_key(state, n_polymers):
    assignments = tuple(state.assignments.tolist())
    if n_polymers == 2:
        swapped = tuple(1 - value if value >= 0 else value for value in assignments)
        assignments = min(assignments, swapped)
    return assignments


def beam_search(
    trace,
    positions,
    scorer,
    n_polymers,
    beam_width,
    rejection_cost,
    minimum_polymer_size,
):
    if n_polymers not in {1, 2}:
        raise ValueError("Only one or two polymers are supported")
    if beam_width < 2:
        raise ValueError("beam_width must be at least 2")
    if rejection_cost < 0:
        raise ValueError("rejection_cost must be non-negative")

    coords = np.column_stack((trace["x"], trace["y"], trace["z"])).astype(float)
    barcodes = np.asarray(trace["Barcode #"])
    ordered_barcodes = sorted(
        set(barcodes), key=lambda barcode: np.min(positions[barcodes == barcode])
    )
    initial = BeamState(
        0.0,
        np.full(len(trace), -1, dtype=int),
        tuple([] for _ in range(n_polymers)),
        (0,) * n_polymers,
    )
    beam = [initial]
    for barcode in ordered_barcodes:
        candidates = np.flatnonzero(barcodes == barcode).tolist()
        choices = [None] + candidates
        barcode_states = []
        for state in beam:
            for selected in product(choices, repeat=n_polymers):
                assigned = [index for index in selected if index is not None]
                if len(set(assigned)) != len(assigned):
                    continue
                assignments = state.assignments.copy()
                histories = tuple(list(history) for history in state.histories)
                counts = list(state.counts)
                cost = state.score + rejection_cost * (len(candidates) - len(assigned))
                for polymer, index in enumerate(selected):
                    if index is None:
                        continue
                    cost += scorer.transition_cost(
                        coords[index], positions[index], histories[polymer]
                    )
                    assignments[index] = polymer
                    histories[polymer].append((positions[index], coords[index]))
                    counts[polymer] += 1
                barcode_states.append(
                    BeamState(cost, assignments, histories, tuple(counts))
                )
        unique = {}
        for state in sorted(barcode_states, key=lambda item: item.score):
            key = _canonical_key(state, n_polymers)
            if key not in unique:
                unique[key] = state
            if len(unique) >= beam_width:
                break
        beam = list(unique.values())

    valid = [
        state
        for state in beam
        if all(count >= minimum_polymer_size for count in state.counts)
    ]
    return sorted(valid or beam, key=lambda state: state.score)
