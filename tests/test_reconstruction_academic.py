from io import BytesIO
from zipfile import ZipFile

import numpy as np
import pandas as pd

from src.reconstruction import (
    edges_from_dataframe,
    grain_graph_reconstruction,
    neighbor_voting_reconstruction,
    operator_groupoid_reconstruction,
    orientation_relationship_presets,
    orientations_from_dataframe,
    symmetry_group,
    synthetic_parent_reconstruction_demo,
)
from src.reconstruction_academic import (
    academic_export_zip,
    boundary_consensus_table,
    build_agreement_summary,
    operator_edge_table,
    operator_frequency_table,
    variant_frequency_table,
)


def _fixture():
    df, adj, _ = synthetic_parent_reconstruction_demo(
        "Kurdjumov–Sachs (FCC parent → BCC child)", n_per_parent=5, noise_deg=0.2, seed=19
    )
    gids, ori = orientations_from_dataframe(df)
    edges = edges_from_dataframe(adj, gids)
    p = orientation_relationship_presets()["Kurdjumov–Sachs (FCC parent → BCC child)"]
    ps = symmetry_group(p.parent_symmetry)
    cs = symmetry_group(p.child_symmetry)
    vote = neighbor_voting_reconstruction(gids, ori, edges, p.matrix_child_to_parent, cs, ps, sigma_deg=2.5, merge_deg=5.0)
    graph = grain_graph_reconstruction(gids, ori, edges, p.matrix_child_to_parent, cs, ps, sigma_deg=2.5, inflation=1.6)
    op = operator_groupoid_reconstruction(gids, ori, edges, p.matrix_child_to_parent, cs, ps, operator_tol_deg=5.0, parent_consistency_deg=5.0)
    return df, adj, gids, ori, edges, p, ps, cs, vote, graph, op


def test_candidate_ambiguity_diagnostics_are_exported_per_daughter():
    *_, vote, _, _ = _fixture()
    cols = set(vote.table.columns)
    assert {"second_best_candidate_fit_deg", "candidate_separation_deg", "candidate_count", "absolute_fit_support", "separation_support"}.issubset(cols)
    assert np.all(vote.table["candidate_count"] >= 1)
    assert np.all(vote.table["candidate_separation_deg"] >= 0)
    assert np.all((vote.table["absolute_fit_support"] >= 0) & (vote.table["absolute_fit_support"] <= 1))


def test_academic_agreement_matrices_and_boundary_consensus():
    df, adj, gids, ori, edges, p, ps, cs, vote, graph, op = _fixture()
    results = {"vote": vote, "graph": graph, "operator": op}
    a = build_agreement_summary(results, edges, ps)
    for mat in (a.ari, a.nmi, a.boundary_jaccard, a.matched_parent_orientation_deg):
        assert mat.shape == (3, 3)
    np.testing.assert_allclose(np.diag(a.ari), 1.0)
    np.testing.assert_allclose(np.diag(a.nmi), 1.0)
    np.testing.assert_allclose(np.diag(a.boundary_jaccard), 1.0)
    consensus = boundary_consensus_table(results, edges, gids)
    assert len(consensus) == len(edges)
    assert consensus["boundary_consensus_fraction"].between(0, 1).all()


def test_operator_edge_and_frequency_tables_are_arpge_style_diagnostics():
    df, adj, gids, ori, edges, p, ps, cs, vote, graph, op = _fixture()
    tbl = operator_edge_table(ori, edges, gids, p.matrix_child_to_parent, ps, cs, 5.0, result=op)
    assert len(tbl) == len(edges)
    assert {"nearest_theoretical_operator_id", "operator_residual_deg", "within_operator_tolerance", "reconstructed_same_parent"}.issubset(tbl.columns)
    assert (tbl["operator_residual_deg"] >= 0).all()
    freq = operator_frequency_table(tbl)
    assert {"operator_id", "edges", "edge_fraction_pct", "mean_operator_residual_deg"}.issubset(freq.columns)


def test_academic_export_bundle_contains_inputs_comparisons_and_method_tables():
    df, adj, gids, ori, edges, p, ps, cs, vote, graph, op = _fixture()
    results = {"vote": vote, "operator": op}
    agreement = build_agreement_summary(results, edges, ps)
    consensus = boundary_consensus_table(results, edges, gids)
    parent_tables = {m: pd.DataFrame({"parent_id": sorted(r.parent_orientations)}) for m, r in results.items()}
    daughter_tables = {m: r.table.copy() for m, r in results.items()}
    variant_tables = {m: variant_frequency_table(r, df) for m, r in results.items()}
    operator_tables = {"operator": operator_edge_table(ori, edges, gids, p.matrix_child_to_parent, ps, cs, 5.0, result=op)}
    comp = pd.DataFrame([{"method": m} for m in results])
    raw = academic_export_zip(
        source_df=df, adjacency_df=adj, results=results, comparison=comp,
        agreement=agreement, parent_tables=parent_tables, daughter_tables=daughter_tables,
        variant_tables=variant_tables, operator_tables=operator_tables,
        metadata={"test": True}, methods_text="academic methods", boundary_consensus=consensus,
    )
    with ZipFile(BytesIO(raw)) as z:
        names = set(z.namelist())
        assert "metadata.json" in names
        assert "input/daughter_grains.csv" in names
        assert "comparison/ARI_clustering_agreement.csv" in names
        assert "comparison/prior_parent_boundary_consensus.csv" in names
        assert any(n.endswith("operator_edges.csv") for n in names)
        assert z.testzip() is None
