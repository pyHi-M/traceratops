import numpy as np
import pytest
from astropy.table import Table

from traceratops.trace_3way_coloc import (
    _trace_barcode_arrays,
    _weighted_pair_counts,
    bootstrap_threeway_colocalization,
    compute_threeway_colocalization,
)


def _sample_trace_table():
    return Table(
        {
            "Trace_ID": [1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4],
            "Barcode #": [1, 2, 3, 1, 2, 3, 2, 3, 1, 2, 4],
            "x": [0.0, 0.1, 0.4, 0.0, 0.1, 0.1, 0.0, 0.0, 0.0, 0.1, 0.1],
            "y": [0.0] * 11,
            "z": [0.0] * 11,
        }
    )


def test_compute_threeway_colocalization_counts_partner_presence_and_colocalization():
    frequencies = compute_threeway_colocalization(
        _sample_trace_table(), anchor_barcode=1, distance_cutoff=0.2
    )

    assert frequencies[(2, 3)] == pytest.approx(1 / 2)
    assert frequencies[(2, 4)] == pytest.approx(1)
    assert frequencies[(3, 4)] == 0


def test_weighted_pair_counts_preserves_duplicate_bootstrap_draws():
    barcodes, present, colocated = _trace_barcode_arrays(
        _sample_trace_table(), anchor_barcode=1, distance_cutoff=0.2
    )
    weights = np.array([2, 0, 0, 1])

    total_counts, colocated_counts = _weighted_pair_counts(
        present, colocated, weights=weights
    )

    barcode_to_idx = {barcode: idx for idx, barcode in enumerate(barcodes)}
    idx_2 = barcode_to_idx[2]
    idx_3 = barcode_to_idx[3]
    assert total_counts[idx_2, idx_3] == 3
    assert colocated_counts[idx_2, idx_3] == 1


def test_bootstrap_threeway_colocalization_returns_all_pairs():
    np.random.seed(0)
    means, sems = bootstrap_threeway_colocalization(
        _sample_trace_table(), anchor_barcode=1, distance_cutoff=0.2, n_bootstrap=3
    )

    assert set(means) == {(2, 3), (2, 4), (3, 4)}
    assert set(sems) == set(means)
