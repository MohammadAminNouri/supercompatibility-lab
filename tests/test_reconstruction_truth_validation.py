from __future__ import annotations

import numpy as np
from sklearn.metrics import adjusted_rand_score

from src.reconstruction import (
    edges_from_dataframe,
    grain_graph_reconstruction,
    operator_groupoid_reconstruction,
    orientations_from_dataframe,
    symmetry_group,
    synthetic_parent_reconstruction_demo_from_or,
    variant_graph_reconstruction,
)
from src.reconstruction_academic import known_truth_validation_metrics
from src.retransformation import NiTiAQLattice, niti_ct_otsuka_ren_orientation_relationship


def _ct_case(noise=0.35):
    lat = NiTiAQLattice(3.01, 2.898, 4.108, 4.646, 97.78)
    r = niti_ct_otsuka_ren_orientation_relationship(lat)
    ps, cs = symmetry_group("cubic"), symmetry_group("monoclinic")
    df, adj, _ = synthetic_parent_reconstruction_demo_from_or(r, ps, cs, n_per_parent=6, noise_deg=noise, seed=11)
    gids, ori = orientations_from_dataframe(df)
    edges = edges_from_dataframe(adj, gids)
    return r, ps, cs, df, gids, ori, edges


def test_ct_synthetic_truth_is_generated_and_recovered_with_same_selected_or():
    r, ps, cs, df, gids, ori, edges = _ct_case()
    methods = [
        variant_graph_reconstruction(gids, ori, edges, r, cs, ps, sigma_deg=2.5, inflation=1.6),
        grain_graph_reconstruction(gids, ori, edges, r, cs, ps, sigma_deg=2.5, inflation=1.6),
        operator_groupoid_reconstruction(gids, ori, edges, r, cs, ps, operator_tol_deg=5.0, parent_consistency_deg=5.0),
    ]
    truth = df["true_parent_id"].to_numpy()
    for res in methods:
        tv = known_truth_validation_metrics(res, truth, edges)
        assert tv["reconstructed_parent_count"] == 2
        assert tv["truth_ARI"] == 1.0
        assert tv["truth_completeness"] == 1.0
        assert tv["truth_homogeneity"] == 1.0
        assert tv["truth_boundary_F1"] == 1.0
        assert tv["validation_status"].startswith("PASS")


def test_singleton_oversegmentation_is_not_misreported_as_accurate():
    r, ps, cs, df, gids, ori, edges = _ct_case(noise=0.0)
    # Reuse a valid result object, but replace clustering labels with one singleton per daughter.
    base = variant_graph_reconstruction(gids, ori, edges, r, cs, ps, sigma_deg=2.5, inflation=1.6)
    from dataclasses import replace
    singleton_ids = np.arange(1, len(gids) + 1, dtype=int)
    bad = replace(base, parent_ids=singleton_ids, diagnostics={**base.diagnostics, "n_reconstructed_parents": len(gids)})
    truth = df["true_parent_id"].to_numpy()
    tv = known_truth_validation_metrics(bad, truth, edges)
    assert adjusted_rand_score(truth, singleton_ids) == 0.0
    assert tv["truth_ARI"] == 0.0
    assert tv["truth_completeness"] < 0.5
    assert tv["truth_boundary_precision"] == 1 / 15
    assert tv["truth_boundary_recall"] == 1.0
    assert tv["truth_boundary_F1"] == 0.125
    assert tv["singleton_parent_fraction"] == 1.0
    assert tv["validation_status"].startswith("FAIL")


def test_validation_adjacency_contains_explicit_truth_boundary_marker():
    r, ps, cs, df, gids, ori, edges = _ct_case()
    # The synthetic topology contains one deliberate inter-parent edge.
    _, adj, _ = synthetic_parent_reconstruction_demo_from_or(r, ps, cs, n_per_parent=6, noise_deg=0.35, seed=11)
    assert "true_parent_boundary" in adj.columns
    assert int(adj["true_parent_boundary"].sum()) == 1
