import gc
import importlib.util
import sys
import weakref
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT=Path(__file__).parents[1]/"scripts"/"analyze_trace_curator_attribution_benchmark.py"
SPEC=importlib.util.spec_from_file_location("attribution_benchmark",SCRIPT)
analysis=importlib.util.module_from_spec(SPEC); sys.modules[SPEC.name]=analysis
SPEC.loader.exec_module(analysis)


def reference_frame(n=30):
    return pd.DataFrame([
        {"simulation_id":"train","Trace_ID":f"t{trace}","Spot_ID":f"s{trace}-{barcode}","Barcode":barcode,"x":float(barcode)+trace*.001,"y":0.,"z":0.,"is_corrupted":False}
        for trace in range(n) for barcode in range(1,6)
    ])


def toy_trace(outlier=50., labels=True):
    return pd.DataFrame({"simulation_id":"eval","Trace_ID":"t","Spot_ID":[f"e{x}" for x in range(5)],"Barcode":[1,2,3,4,5],"x":[1.,2.,outlier,4.,5.],"y":0.,"z":0.,"is_corrupted":[False,False,labels,False,False],"selected_for_corruption":[False,False,labels,False,False],"condition":"c","replicate":1,"seed":99,"detection_efficiency":1.,"displacement":1.,"injected_displacement_um":[0,0,1,0,0]})


def test_reference_and_nested_fold_provenance_do_not_leak():
    reference, calibration, evaluation = analysis.fold_seed_sets([1,2,3,4],4,3)
    assert reference == {1,2}
    assert calibration == {3}
    assert evaluation == {4}
    assert not reference & calibration and not reference & evaluation
    frame=reference_frame(); held=toy_trace()
    model=analysis.fit_reference(frame[frame.simulation_id!="eval"],1)
    before=model.separation[1].copy(); analysis.score_observations(held,model)
    np.testing.assert_array_equal(before,model.separation[1])
    assert set(frame.simulation_id)=={"train"}


def test_obvious_displaced_point_is_attributed_not_its_innocent_neighbor():
    scored=analysis.score_trace(toy_trace(),analysis.fit_reference(reference_frame(),1))
    star=scored[scored.model=="star_top2_mean_all"].set_index("Barcode")
    assert star.loc[3,"score"]>star.loc[2,"score"]
    assert star.loc[3,"score"]>star.loc[4,"score"]
    enriched=scored.merge(toy_trace()[["Spot_ID","is_corrupted"]],on="Spot_ID")
    row=analysis.rank_trace(enriched[enriched.model=="star_top2_mean_all"])
    assert row["top1"] and row["rank"]==1


def test_collateral_distance_binning():
    assert analysis.collateral_distance(9,10)=="1"
    assert analysis.collateral_distance(8,10)=="2"
    assert analysis.collateral_distance(7,10)=="3"
    assert analysis.collateral_distance(6,10)==">=4"


def test_trace_ranking_metrics_include_reciprocal_rank():
    frame=pd.DataFrame({"score":[3.,5.,1.],"is_corrupted":[False,True,False]})
    row=analysis.rank_trace(frame)
    assert row=={"scoreable":True,"rank":1,"top1":True,"top2":True,"reciprocal_rank":1.}


def test_ranking_summary_reports_mean_median_and_unconditional_top_k():
    frame=pd.DataFrame({
        "model":"m","detection_efficiency":.5,"displacement":.4,
        "scoreable":[True,True,False],"rank":[1.,3.,np.nan],
        "top1":[True,False,False],"top2":[True,False,False],
        "reciprocal_rank":[1.,1/3,0.],
    })
    row=analysis.summarize_ranking(frame).iloc[0]
    assert row.top1_attribution_accuracy == 1/3
    assert row.mean_true_target_rank == row.median_true_target_rank == 2
    assert row.mean_reciprocal_rank == (1+1/3)/3


def test_bridge_reports_missing_flank_and_insufficient_reference():
    trace=toy_trace().sort_values("Barcode").reset_index(drop=True)
    model=analysis.fit_reference(reference_frame(),1)
    assert analysis.bridge_score(trace,0,model)[1]=="missing_flank"
    insufficient=analysis.fit_reference(reference_frame(1),20)
    assert analysis.bridge_score(trace,2,insufficient)[1]=="insufficient_reference"


def test_edge_reference_support_is_explicit():
    scored=analysis.score_trace(toy_trace(),analysis.fit_reference(reference_frame(1),20))
    non_bridge=scored[scored.model!="empirical_two_flank_bridge"]
    assert non_bridge.score.isna().all()
    assert set(non_bridge.score_status)=={"insufficient_reference"}


def test_condition_frames_are_releasable_for_bounded_memory_streaming():
    refs=[]; model=analysis.fit_reference(reference_frame(),1)
    for _ in range(3):
        frame=analysis.score_trace(toy_trace(),model); refs.append(weakref.ref(frame)); del frame; gc.collect()
    assert all(ref() is None for ref in refs)


def test_zero_displacement_labels_do_not_change_scores():
    model=analysis.fit_reference(reference_frame(),1)
    first=analysis.score_trace(toy_trace(labels=True),model)
    second=analysis.score_trace(toy_trace(labels=False),model)
    np.testing.assert_allclose(first.score,second.score,equal_nan=True)
