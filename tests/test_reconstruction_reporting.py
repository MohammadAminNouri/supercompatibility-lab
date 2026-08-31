import numpy as np

from src.reconstruction import (
    neighbor_voting_reconstruction,
    orientation_relationship_presets,
    orientations_from_dataframe,
    edges_from_dataframe,
    symmetry_group,
    synthetic_parent_reconstruction_demo,
)
from src.reconstruction_reporting import (
    audit_reconstruction_input,
    daughter_assignment_table,
    method_requirements_table,
    parent_summary_table,
    reconstruction_quality_summary,
)


def _result():
    ori, edf, _ = synthetic_parent_reconstruction_demo("Kurdjumov–Sachs (FCC parent → BCC child)", n_per_parent=5, noise_deg=0.2, seed=13)
    ids, gos = orientations_from_dataframe(ori)
    edges = edges_from_dataframe(edf, ids)
    p = orientation_relationship_presets()["Kurdjumov–Sachs (FCC parent → BCC child)"]
    r = neighbor_voting_reconstruction(
        ids, gos, edges, p.matrix_child_to_parent,
        symmetry_group(p.child_symmetry), symmetry_group(p.parent_symmetry),
    )
    return ori, edges, r


def test_method_requirements_are_explicit():
    df = method_requirements_table()
    assert len(df) == 5
    assert {"minimum input", "main controls", "important limitation"}.issubset(df.columns)


def test_input_audit_detects_graph():
    ori, edges, _ = _result()
    a = audit_reconstruction_input(ori, edges, "measured/shared-boundary")
    assert a.n_grains == len(ori)
    assert a.n_edges == len(edges)
    assert a.connected_fraction >= 0.9


def test_parent_and_daughter_tables_are_meaningful():
    ori, _, r = _result()
    d = daughter_assignment_table(r, ori)
    p = parent_summary_table(r, ori)
    assert "OR_fit_misorientation_deg" in d.columns
    assert "support_score_0_to_1" in d.columns
    assert "parent_id" in d.columns
    assert len(p) == r.diagnostics["n_reconstructed_parents"]
    assert np.all(p["supporting_daughter_grains"] >= 1)


def test_quality_summary_ranges():
    ori, _, r = _result()
    q = reconstruction_quality_summary(r)
    assert 0 <= q["fraction_fit_le_2_5deg"] <= 1
    assert 0 <= q["fraction_support_lt_0_45"] <= 1
