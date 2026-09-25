import numpy as np
import pytest
from astropy.table import Table

from traceratops import trace_splitter
from traceratops.core.chromatin_trace_table import ChromatinTraceTable
from traceratops.core.trace_resolver import (
    ClassificationThresholds,
    EmpiricalDistanceModel,
    LikelihoodMultiplicityClassifier,
    TraceClassification,
    TraceMultiplicity,
    TraceResolver,
    assess_likelihood_classifier_reliability,
    classify_trace,
)
from traceratops.core.trace_resolver.model import DistanceStatistics


def _trace(rows, trace_id="merged"):
    return Table(
        rows=[
            (f"spot-{i}", trace_id, barcode, x, y, z)
            for i, (barcode, x, y, z) in enumerate(rows)
        ],
        names=("Spot_ID", "Trace_ID", "Barcode #", "x", "y", "z"),
    )


def _linear_model(maximum_separation=20):
    model = EmpiricalDistanceModel(minimum_observations=1, variance_floor=0.01)
    model.statistics = {
        float(delta): DistanceStatistics(float(delta), 0.1, 10)
        for delta in range(1, maximum_separation + 1)
    }
    model.pooled = DistanceStatistics(1.0, 0.1, 10)
    return model


def test_multiplicity_distinguishes_one_frequent_barcode_from_many_repeats():
    one_frequent = TraceMultiplicity.from_barcodes([1, 1, 1, 1, 1, 1, 2, 3])
    many_repeated = TraceMultiplicity.from_barcodes(np.repeat(np.arange(15), 2))
    thresholds = ClassificationThresholds(3, 0.3, "both")

    assert one_frequent.n_repeated_barcodes == 1
    assert one_frequent.maximum_barcode_multiplicity == 6
    assert classify_trace(one_frequent, thresholds) == TraceClassification.CLEAN_ONE
    assert classify_trace(many_repeated, thresholds) == TraceClassification.RESOLVE_TWO


def test_clean_single_trace_is_unchanged_and_preserves_spot_ids():
    original = _trace([(1, 0, 0, 0), (2, 1, 0, 0), (3, 2, 0, 0)])
    trace_table = ChromatinTraceTable()
    trace_table.data = original.copy()

    diagnostics = trace_splitter.resolve_traces(trace_table)

    assert diagnostics[0]["status"] == "unchanged"
    assert list(trace_table.data["Spot_ID"]) == list(original["Spot_ID"])
    assert list(trace_table.data["Trace_ID"]) == ["merged"] * 3


def test_threshold_classifier_has_no_likelihood_reliability_warning(capsys):
    table = ChromatinTraceTable()
    table.data = _trace([(1, 0, 0, 0), (2, 1, 0, 0), (3, 2, 0, 0)])

    diagnostics = trace_splitter.resolve_traces(table)

    assert "likelihood multiplicity classification" not in capsys.readouterr().out
    assert diagnostics[0]["classifier_reliability_level"] == ""
    assert "classifier_reliability_message" not in diagnostics.meta


def test_likelihood_reliability_warning_once_and_diagnostics(monkeypatch, capsys):
    rows = []
    for trace_id, offset in (("one", 0.0), ("two", 10.0)):
        rows.extend(
            (f"{trace_id}-{barcode}", trace_id, barcode, offset + barcode, 0.0, 0.0)
            for barcode in range(1, 6)
        )
    table = ChromatinTraceTable()
    table.data = Table(
        rows=rows, names=("Spot_ID", "Trace_ID", "Barcode #", "x", "y", "z")
    )
    original_fit = LikelihoodMultiplicityClassifier.fit

    def fit_in_high_risk_regime(classifier, counts, barcode_ids=None):
        fitted = original_fit(classifier, counts, barcode_ids)
        fitted.detection_efficiency = 0.3
        fitted.doublet_prior = 0.2
        fitted.reliability_assessment = assess_likelihood_classifier_reliability(
            len(fitted.barcode_ids), fitted.detection_efficiency, fitted.doublet_prior
        )
        return fitted

    monkeypatch.setattr(
        LikelihoodMultiplicityClassifier, "fit", fit_in_high_risk_regime
    )

    diagnostics = trace_splitter.resolve_traces(
        table, multiplicity_classifier="likelihood"
    )

    warning = "Warning: likelihood multiplicity classification"
    assert capsys.readouterr().out.count(warning) == 1
    assert set(diagnostics["classifier_n_barcodes"]) == {5}
    assert set(diagnostics["classifier_expected_detected_barcodes_per_polymer"]) == {
        1.5
    }
    assert set(diagnostics["classifier_reliability_level"]) == {"high-risk"}
    assert (
        "poor singlet/doublet identifiability"
        in diagnostics.meta["classifier_reliability_message"]
    )


@pytest.mark.parametrize("duplicate_count", [2, 6])
def test_one_polymer_cleanup_keeps_compatible_duplicate(duplicate_count):
    rows = [(1, 0, 0, 0)]
    rows.extend((2, 1 if i == 0 else 10 + i, 0, 0) for i in range(duplicate_count))
    rows.append((3, 2, 0, 0))
    trace = _trace(rows)
    resolver = TraceResolver(
        _linear_model(), rejection_cost=20, minimum_confidence=0.01
    )

    result = resolver.resolve(trace, n_polymers=1)

    assert result.status == "cleaned"
    kept = trace[result.assignments == 0]
    assert list(kept["Barcode #"]) == [1, 2, 3]
    assert kept["x"][1] == 1
    assert result.n_unassigned == duplicate_count - 1


def test_one_polymer_cleanup_forces_unique_barcode_with_poor_geometry():
    trace = _trace(
        [
            (1, 0, 0, 0),
            (2, 1, 0, 0),
            (2, 5, 0, 0),
            # This unique localization has an intentionally extreme residual.
            (3, 100, 0, 0),
        ]
    )
    result = TraceResolver(
        _linear_model(), rejection_cost=0, minimum_confidence=0
    ).resolve(trace, n_polymers=1)

    assert result.status == "cleaned"
    assert result.assignments[3] == 0
    assert list(trace[result.assignments == 0]["Barcode #"]) == [1, 2, 3]


def test_overlapping_traces_resolve_from_barcode_continuity():
    rows = []
    for barcode in range(1, 7):
        rows.append((barcode, float(barcode), 0.0, 0.0))
        rows.append((barcode, float(barcode), 0.2, 0.0))
    trace = _trace(rows)
    resolver = TraceResolver(
        _linear_model(),
        beam_width=300,
        rejection_cost=20,
        minimum_confidence=0.001,
    )

    result = resolver.resolve(trace, n_polymers=2)

    assert result.status == "split"
    assert result.inferred_n_polymers == 2
    for polymer in (0, 1):
        assigned = trace[result.assignments == polymer]
        assert len(assigned) == 6
        assert len(set(assigned["Barcode #"])) == 6


def test_spatially_separated_traces_resolve():
    rows = []
    for barcode in range(1, 6):
        rows.append((barcode, float(barcode), 0.0, 0.0))
        rows.append((barcode, float(barcode), 10.0, 0.0))
    result = TraceResolver(
        _linear_model(),
        beam_width=200,
        rejection_cost=20,
        minimum_confidence=0.001,
    ).resolve(_trace(rows), n_polymers=2)

    assert result.status == "split"
    assert sorted(np.sum(result.assignments == polymer) for polymer in (0, 1)) == [
        5,
        5,
    ]


def test_split_table_preserves_spot_ids_and_writes_complete_diagnostics(
    monkeypatch, capsys
):
    rows = []
    for barcode in range(1, 5):
        rows.append((barcode, float(barcode), 0.0, 0.0))
        rows.append((barcode, float(barcode), 10.0, 0.0))
    original = _trace(rows)
    table = ChromatinTraceTable()
    table.data = original.copy()
    identifiers = iter(("polymer-one", "polymer-two"))
    monkeypatch.setattr(trace_splitter, "generate_unique_id", lambda: next(identifiers))

    diagnostics = trace_splitter.resolve_traces(
        table,
        thresholds=ClassificationThresholds(2, 0.2, "both"),
        rejection_cost=20,
        minimum_confidence=0.001,
        model_minimum_observations=1,
    )

    assert set(table.data["Spot_ID"]) == set(original["Spot_ID"])
    assert set(table.data["Trace_ID"]) == {"polymer-one", "polymer-two"}
    assert diagnostics[0]["status"] == "split"
    assert diagnostics[0]["requested_n_polymers"] == 2
    assert diagnostics[0]["inferred_n_polymers"] == 2
    assert diagnostics.meta["genomic_source"] == "barcode"
    assert diagnostics.meta["fallback_source"] == "nearest_candidate_pairs"
    assert "using nearest candidate pairs" in capsys.readouterr().out
    expected = {
        "input_trace_id",
        "status",
        "n_unique_barcodes",
        "n_repeated_barcodes",
        "fraction_repeated_barcodes",
        "maximum_barcode_multiplicity",
        "n_excess_detections",
        "radius_of_gyration",
        "requested_n_polymers",
        "inferred_n_polymers",
        "best_score",
        "alternative_score",
        "raw_score_gap",
        "confidence_score",
        "n_unassigned",
        "fraction_unassigned",
        "n_ambiguous_barcodes",
    }
    assert expected.issubset(diagnostics.colnames)


def test_low_efficiency_polymers_can_have_unequal_completeness_and_gaps():
    trace = _trace(
        [
            (1, 1.0, 0.0, 0.0),
            (1, 1.0, 0.2, 0.0),
            (4, 4.0, 0.0, 0.0),
            (4, 4.0, 0.2, 0.0),
            (8, 8.0, 0.0, 0.0),
            (12, 12.0, 0.0, 0.0),
            (12, 12.0, 0.2, 0.0),
        ]
    )
    result = TraceResolver(
        _linear_model(),
        beam_width=300,
        rejection_cost=20,
        minimum_confidence=0.001,
    ).resolve(trace, n_polymers=2)

    assert result.status == "split"
    lengths = sorted(np.sum(result.assignments == polymer) for polymer in (0, 1))
    assert lengths == [3, 4]


def test_ambiguous_two_polymer_solution_is_removed():
    rows = []
    for barcode in range(1, 5):
        rows.extend([(barcode, float(barcode), 0, 0)] * 2)
    trace = _trace(rows)

    result = TraceResolver(
        _linear_model(), beam_width=200, rejection_cost=20, minimum_confidence=0.1
    ).resolve(trace, n_polymers=2)

    assert result.status == "removed_ambiguous"
    assert np.all(result.assignments == -1)


def test_ambiguous_one_polymer_barcode_is_rejected_but_trace_is_retained():
    trace = _trace(
        [
            (1, 0, 0, 0),
            (2, 1, 0.1, 0),
            (2, 1, -0.1, 0),
            (3, 2, 0, 0),
        ]
    )
    result = TraceResolver(
        _linear_model(), beam_width=100, rejection_cost=20, minimum_confidence=0.1
    ).resolve(trace, n_polymers=1)

    assert result.status == "cleaned"
    assert result.n_ambiguous_barcodes == 1
    assert np.all(result.assignments[np.asarray(trace["Barcode #"]) == 2] == -1)
    assert np.sum(result.assignments == 0) == 2


@pytest.mark.parametrize("history_mode", ["nearest", "multi"])
@pytest.mark.parametrize("distance_score", ["residual", "likelihood"])
def test_all_history_and_distance_score_combinations(history_mode, distance_score):
    trace = _trace([(1, 0, 0, 0), (2, 1, 0, 0), (3, 2, 0, 0)])
    result = TraceResolver(
        _linear_model(),
        history_mode=history_mode,
        distance_score=distance_score,
        rejection_cost=20,
    ).resolve(trace, n_polymers=1)
    assert result.status == "cleaned"


def test_empirical_model_prefers_genomic_coordinates_and_has_sparse_fallback():
    table = Table(
        rows=[
            ("a", 1, "chr1", 0, 10, 0.0, 0.0, 0.0),
            ("a", 2, "chr1", 100, 110, 1.0, 0.0, 0.0),
            ("b", 1, "chr1", 0, 10, 0.0, 0.0, 0.0),
            ("b", 2, "chr1", 100, 110, 1.2, 0.0, 0.0),
        ],
        names=(
            "Trace_ID",
            "Barcode #",
            "Chrom",
            "Chrom_Start",
            "Chrom_End",
            "x",
            "y",
            "z",
        ),
    )
    model = EmpiricalDistanceModel(minimum_observations=3).fit(table)

    assert model.genomic_source == "coordinates"
    assert model.pooled.n_observations == 2
    assert model.for_separation(1000) == model.pooled
