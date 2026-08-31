from __future__ import annotations

from io import BytesIO
from zipfile import ZipFile

import numpy as np

from src.reconstruction import (
    edges_from_dataframe,
    neighbor_voting_reconstruction,
    orientation_relationship_presets,
    orientations_from_dataframe,
    symmetry_group,
    synthetic_parent_reconstruction_demo,
)
from src.retransformation import (
    NiTiAQLattice,
    branch_switch_catalog,
    cycle_closure_table,
    cycle_evidence_zip,
    match_new_cycle_daughters,
    niti_aq_orientation_relationship,
    niti_aq_parallelism_residuals,
    niti_ct_otsuka_ren_diagnostics,
    niti_ct_otsuka_ren_orientation_relationship,
    observed_branch_occupancy,
    parent_cycle_summary,
    regenerated_variant_library,
)


def _ks_demo(noise_deg: float = 0.25):
    p = orientation_relationship_presets()["Kurdjumov–Sachs (FCC parent → BCC child)"]
    df, edge_df, _ = synthetic_parent_reconstruction_demo(
        "Kurdjumov–Sachs (FCC parent → BCC child)", n_per_parent=6, noise_deg=noise_deg, seed=11
    )
    gids, ori = orientations_from_dataframe(df)
    edges = edges_from_dataframe(edge_df, gids)
    ps = symmetry_group(p.parent_symmetry)
    cs = symmetry_group(p.child_symmetry)
    r = neighbor_voting_reconstruction(gids, ori, edges, p.matrix_child_to_parent, cs, ps)
    return p, df, gids, ori, ps, cs, r


def test_metric_aware_niti_aq_or_is_proper_and_satisfies_parallelisms():
    lat = NiTiAQLattice()
    r = niti_aq_orientation_relationship(lat)
    d = niti_aq_parallelism_residuals(lat, r)
    assert np.isclose(np.linalg.det(r), 1.0, atol=1e-12)
    np.testing.assert_allclose(r.T @ r, np.eye(3), atol=1e-12)
    assert d["plane_parallelism_residual_deg"] < 1e-8
    assert d["direction_parallelism_residual_deg"] < 1e-8




def test_ct_otsuka_ren_model_derived_or_is_proper_and_exposes_exact_correspondence():
    from src.core import C_A_TO_M, C_M_TO_A
    lat = NiTiAQLattice()
    r = niti_ct_otsuka_ren_orientation_relationship(lat)
    d = niti_ct_otsuka_ren_diagnostics(lat, r)
    assert np.isclose(np.linalg.det(r), 1.0, atol=1e-12)
    np.testing.assert_allclose(r.T @ r, np.eye(3), atol=1e-12)
    np.testing.assert_allclose(d["correspondence_A_to_M"], C_A_TO_M, atol=0)
    np.testing.assert_allclose(d["correspondence_M_to_A"], C_M_TO_A, atol=0)
    assert np.asarray(d["F_child_from_parent"]).shape == (3, 3)
    assert np.all(np.asarray(d["principal_stretches_sorted"]) > 0)
    assert d["misorientation_to_natural_AQ_OR_deg"] > 0


def test_ct_otsuka_ren_model_derived_or_has_twelve_b19p_branches():
    r = niti_ct_otsuka_ren_orientation_relationship(NiTiAQLattice())
    ps = symmetry_group("cubic")
    cs = symmetry_group("monoclinic")
    library, cache = regenerated_variant_library({1: np.eye(3)}, r, ps, cs)
    assert len(cache[1]) == 12
    assert len(library) == 12

def test_niti_aq_has_twelve_symmetry_distinct_b19p_branches():
    r = niti_aq_orientation_relationship(NiTiAQLattice())
    ps = symmetry_group("cubic")
    cs = symmetry_group("monoclinic")
    library, cache = regenerated_variant_library({1: np.eye(3)}, r, ps, cs)
    assert len(cache[1]) == 12
    assert len(library) == 12
    assert {"daughter_Bunge_phi1_deg", "daughter_quaternion_qw", "g11", "g33"}.issubset(library.columns)


def test_round_trip_closure_is_small_for_reconstructed_synthetic_data():
    p, df, gids, ori, ps, cs, r = _ks_demo(0.25)
    library, _ = regenerated_variant_library(r.parent_orientations, p.matrix_child_to_parent, ps, cs)
    closure = cycle_closure_table(gids, ori, r, p.matrix_child_to_parent, ps, cs)
    summary = parent_cycle_summary(closure, library)
    assert len(closure) == len(gids)
    assert closure["cycle_closure_misorientation_deg"].mean() < 0.6
    assert closure["cycle_closure_misorientation_deg"].max() < 1.5
    assert len(summary) == 2
    assert (summary["allowed_regenerated_daughter_branches"] == 24).all()


def test_observed_branch_occupancy_accounts_for_all_grains():
    p, df, gids, ori, ps, cs, r = _ks_demo(0.2)
    closure = cycle_closure_table(gids, ori, r, p.matrix_child_to_parent, ps, cs)
    occ = observed_branch_occupancy(closure, df)
    assert int(occ["daughter_grains"].sum()) == len(gids)
    assert (occ["grain_fraction_within_parent_pct"] > 0).all()


def test_branch_switch_catalog_has_zero_diagonal_and_symmetric_angles():
    rcp = niti_aq_orientation_relationship(NiTiAQLattice())
    ps = symmetry_group("cubic")
    cs = symmetry_group("monoclinic")
    tbl = branch_switch_catalog({1: np.eye(3)}, rcp, ps, cs)
    assert len(tbl) == 12 * 12
    diag = tbl[tbl["from_regenerated_branch_id"] == tbl["to_regenerated_branch_id"]]
    assert np.allclose(diag["daughter_orientation_change_deg"], 0.0, atol=1e-8)
    m = tbl.pivot(index="from_regenerated_branch_id", columns="to_regenerated_branch_id", values="daughter_orientation_change_deg").to_numpy()
    np.testing.assert_allclose(m, m.T, atol=1e-10)


def test_independent_new_cycle_exact_generated_variants_match_parent_and_branch():
    rcp = niti_aq_orientation_relationship(NiTiAQLattice())
    ps = symmetry_group("cubic")
    cs = symmetry_group("monoclinic")
    parent_orientations = {1: np.eye(3), 2: __import__("scipy.spatial.transform", fromlist=["Rotation"]).Rotation.from_euler("ZXZ", [23, 41, 17], degrees=True).as_matrix()}
    _, cache = regenerated_variant_library(parent_orientations, rcp, ps, cs)
    measured = [cache[1][3], cache[2][7]]
    out = match_new_cycle_daughters([101, 202], measured, parent_orientations, rcp, ps, cs)
    assert out["matched_reconstructed_parent_id"].tolist() == [1, 2]
    assert out["matched_regenerated_branch_id"].tolist() == [4, 8]
    assert out["new_cycle_OR_library_fit_deg"].max() < 1e-5


def test_known_parent_id_restricts_new_cycle_matching():
    rcp = niti_aq_orientation_relationship(NiTiAQLattice())
    ps = symmetry_group("cubic")
    cs = symmetry_group("monoclinic")
    parent_orientations = {1: np.eye(3), 2: __import__("scipy.spatial.transform", fromlist=["Rotation"]).Rotation.from_euler("ZXZ", [23, 41, 17], degrees=True).as_matrix()}
    _, cache = regenerated_variant_library(parent_orientations, rcp, ps, cs)
    out = match_new_cycle_daughters([7], [cache[2][1]], parent_orientations, rcp, ps, cs, known_parent_ids=[2])
    assert int(out.iloc[0]["matched_reconstructed_parent_id"]) == 2
    assert out.iloc[0]["parent_assignment_basis"] == "user-supplied parent_id"


def test_cycle_evidence_zip_contains_all_academic_tables():
    p, df, gids, ori, ps, cs, r = _ks_demo(0.2)
    library, _ = regenerated_variant_library(r.parent_orientations, p.matrix_child_to_parent, ps, cs)
    closure = cycle_closure_table(gids, ori, r, p.matrix_child_to_parent, ps, cs)
    summary = parent_cycle_summary(closure, library)
    occ = observed_branch_occupancy(closure, df)
    switch = branch_switch_catalog(r.parent_orientations, p.matrix_child_to_parent, ps, cs)
    raw = cycle_evidence_zip(
        variant_library=library,
        closure=closure,
        parent_summary=summary,
        occupancy=occ,
        switch_catalog=switch,
        metadata={"OR": p.name},
    )
    with ZipFile(BytesIO(raw)) as z:
        names = set(z.namelist())
    assert {
        "regenerated_daughter_variant_library.csv",
        "measured_round_trip_cycle_closure.csv",
        "parent_cycle_summary.csv",
        "observed_branch_occupancy.csv",
        "regenerated_branch_switch_catalog.csv",
        "cycle_metadata.json",
        "README.txt",
    }.issubset(names)


def test_all_five_reconstruction_families_recover_metric_aware_niti_aq_synthetic_parents():
    from scipy.spatial.transform import Rotation
    from src.reconstruction import (
        grain_graph_reconstruction,
        nucleation_growth_reconstruction,
        operator_groupoid_reconstruction,
        reconstruction_accuracy_against_labels,
        unique_child_variants,
        variant_graph_reconstruction,
    )

    rcp = niti_aq_orientation_relationship(NiTiAQLattice())
    ps = symmetry_group("cubic")
    cs = symmetry_group("monoclinic")
    parents = [np.eye(3), Rotation.from_euler("ZXZ", [31, 47, 83], degrees=True).as_matrix()]
    rng = np.random.default_rng(123)
    gids: list[int] = []
    ori: list[np.ndarray] = []
    truth: list[int] = []
    for pid, gp in enumerate(parents, 1):
        variants = __import__("src.reconstruction", fromlist=["unique_child_variants"]).unique_child_variants(gp, rcp, ps, cs)
        for j in range(6):
            g = variants[(2 * j) % len(variants)]
            axis = rng.normal(size=3); axis /= np.linalg.norm(axis)
            noisy = Rotation.from_rotvec(np.deg2rad(rng.normal(0.0, 0.2)) * axis).as_matrix() @ g
            gids.append(len(gids) + 1); ori.append(noisy); truth.append(pid)
    edges = []
    for off in (0, 6):
        for i in range(off, off + 6):
            edges.append((i, off + (i - off + 1) % 6))
    edges.append((5, 6))

    methods = [
        neighbor_voting_reconstruction,
        grain_graph_reconstruction,
        variant_graph_reconstruction,
        nucleation_growth_reconstruction,
        operator_groupoid_reconstruction,
    ]
    for fn in methods:
        out = fn(gids, ori, edges, rcp, cs, ps)
        assert len(set(out.parent_ids.tolist())) == 2, fn.__name__
        assert reconstruction_accuracy_against_labels(out, np.asarray(truth)) >= 0.99, fn.__name__
        assert float(out.diagnostics["mean_fit_deg"]) < 0.5, fn.__name__


def test_all_five_reconstruction_families_recover_ct_otsuka_ren_synthetic_parents():
    from scipy.spatial.transform import Rotation
    from src.reconstruction import (
        grain_graph_reconstruction, nucleation_growth_reconstruction, operator_groupoid_reconstruction,
        reconstruction_accuracy_against_labels, unique_child_variants, variant_graph_reconstruction,
    )
    rcp = niti_ct_otsuka_ren_orientation_relationship(NiTiAQLattice())
    ps = symmetry_group("cubic")
    cs = symmetry_group("monoclinic")
    parents = [np.eye(3), Rotation.from_euler("ZXZ", [31, 47, 83], degrees=True).as_matrix()]
    rng = np.random.default_rng(321)
    gids, ori, truth = [], [], []
    for pid, gp in enumerate(parents, 1):
        variants = unique_child_variants(gp, rcp, ps, cs)
        for j in range(6):
            g = variants[(2 * j) % len(variants)]
            axis = rng.normal(size=3); axis /= np.linalg.norm(axis)
            noisy = Rotation.from_rotvec(np.deg2rad(rng.normal(0.0, 0.2)) * axis).as_matrix() @ g
            gids.append(len(gids) + 1); ori.append(noisy); truth.append(pid)
    edges = []
    for off in (0, 6):
        for i in range(off, off + 6):
            edges.append((i, off + (i - off + 1) % 6))
    edges.append((5, 6))
    methods = [neighbor_voting_reconstruction, grain_graph_reconstruction, variant_graph_reconstruction, nucleation_growth_reconstruction, operator_groupoid_reconstruction]
    for fn in methods:
        out = fn(gids, ori, edges, rcp, cs, ps)
        assert len(set(out.parent_ids.tolist())) == 2, fn.__name__
        assert reconstruction_accuracy_against_labels(out, np.asarray(truth)) >= 0.99, fn.__name__
        assert float(out.diagnostics["mean_fit_deg"]) < 0.5, fn.__name__
