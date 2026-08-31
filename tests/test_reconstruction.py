import numpy as np
from scipy.spatial.transform import Rotation

from src.reconstruction import (
    bunge_euler_to_matrix,
    cubic_rotations,
    edges_from_dataframe,
    grain_graph_reconstruction,
    matrix_to_bunge_euler,
    neighbor_voting_reconstruction,
    nucleation_growth_reconstruction,
    operator_groupoid_reconstruction,
    orientation_relationship_presets,
    orientations_from_dataframe,
    parent_candidates,
    reconstruction_accuracy_against_labels,
    refine_orientation_relationship,
    rotation_angle_deg,
    symmetry_group,
    synthetic_parent_reconstruction_demo,
    unique_child_variants,
    variant_graph_reconstruction,
)


def test_bunge_roundtrip_rotation_matrix():
    g = bunge_euler_to_matrix(37.0, 61.0, 123.0)
    a, b, c = matrix_to_bunge_euler(g)
    g2 = bunge_euler_to_matrix(a, b, c)
    assert rotation_angle_deg(g @ g2.T) < 1e-8


def test_proper_cubic_group_has_24_rotations():
    group = cubic_rotations()
    assert len(group) == 24
    assert all(np.isclose(np.linalg.det(g), 1.0, atol=1e-12) for g in group)


def test_standard_or_variant_counts():
    presets = orientation_relationship_presets()
    expected = {
        "Kurdjumov–Sachs (FCC parent → BCC child)": 24,
        "Nishiyama–Wassermann (FCC parent → BCC child)": 12,
        "Bain (FCC parent → BCC child)": 3,
        "Pitsch (FCC parent → BCC child)": 12,
        "Burgers (BCC parent → HCP child)": 12,
    }
    for name, n_expected in expected.items():
        p = presets[name]
        ps = symmetry_group(p.parent_symmetry)
        cs = symmetry_group(p.child_symmetry)
        variants = unique_child_variants(np.eye(3), p.matrix_child_to_parent, ps, cs)
        assert len(variants) == n_expected


def test_ks_child_has_24_parent_candidates():
    p = orientation_relationship_presets()["Kurdjumov–Sachs (FCC parent → BCC child)"]
    ps = symmetry_group(p.parent_symmetry)
    cs = symmetry_group(p.child_symmetry)
    child = unique_child_variants(np.eye(3), p.matrix_child_to_parent, ps, cs)[0]
    candidates = parent_candidates(child, p.matrix_child_to_parent, cs, ps)
    assert len(candidates) == 24
    assert min(__import__("src.reconstruction", fromlist=["misorientation_deg"]).misorientation_deg(g, np.eye(3), ps) for g, _ in candidates) < 1e-7


def _demo():
    p = orientation_relationship_presets()["Kurdjumov–Sachs (FCC parent → BCC child)"]
    df, edge_df, _ = synthetic_parent_reconstruction_demo(n_per_parent=6, noise_deg=0.35, seed=11)
    gids, ori = orientations_from_dataframe(df)
    edges = edges_from_dataframe(edge_df, gids)
    ps = symmetry_group(p.parent_symmetry)
    cs = symmetry_group(p.child_symmetry)
    return p, df, gids, ori, edges, ps, cs


def test_all_five_reconstruction_methods_recover_synthetic_parents():
    p, df, gids, ori, edges, ps, cs = _demo()
    funcs = [
        neighbor_voting_reconstruction,
        grain_graph_reconstruction,
        variant_graph_reconstruction,
        nucleation_growth_reconstruction,
        operator_groupoid_reconstruction,
    ]
    for f in funcs:
        r = f(gids, ori, edges, p.matrix_child_to_parent, cs, ps)
        acc = reconstruction_accuracy_against_labels(r, df.true_parent_id.to_numpy())
        assert acc >= 0.99, (f.__name__, acc, r.diagnostics)
        assert len(set(r.parent_ids.tolist())) == 2
        assert r.diagnostics["mean_fit_deg"] < 1.0


def test_or_refinement_recovers_known_two_degree_perturbation():
    p = orientation_relationship_presets()["Kurdjumov–Sachs (FCC parent → BCC child)"]
    df, edge_df, _ = synthetic_parent_reconstruction_demo(n_per_parent=5, noise_deg=0.25, seed=11)
    gids, ori = orientations_from_dataframe(df)
    edges = edges_from_dataframe(edge_df, gids)
    ps = symmetry_group(p.parent_symmetry); cs = symmetry_group(p.child_symmetry)
    axis = np.array([0.3, 0.4, 0.8660254038])
    axis /= np.linalg.norm(axis)
    pert = Rotation.from_rotvec(np.deg2rad(2.0) * axis).as_matrix() @ p.matrix_child_to_parent
    refined, d = refine_orientation_relationship(ori, edges, pert, cs, ps, max_rotation_deg=4.0)
    initial_err = rotation_angle_deg(pert @ p.matrix_child_to_parent.T)
    final_err = rotation_angle_deg(refined @ p.matrix_child_to_parent.T)
    assert initial_err > 1.9
    assert final_err < 0.4
    assert d["objective_after"] < d["objective_before"]


def test_standard_or_matrices_satisfy_declared_parallelisms():
    presets = orientation_relationship_presets()
    cases = {
        "Kurdjumov–Sachs (FCC parent → BCC child)": ([1,1,1],[1,-1,0],[1,1,0],[1,-1,1]),
        "Nishiyama–Wassermann (FCC parent → BCC child)": ([1,1,1],[0,1,-1],[1,1,0],[0,0,1]),
        "Bain (FCC parent → BCC child)": ([0,1,0],[0,0,1],[0,1,0],[1,0,1]),
        "Pitsch (FCC parent → BCC child)": ([1,0,0],[0,1,1],[1,1,0],[1,-1,1]),
        "Burgers (BCC parent → HCP child)": ([1,1,0],[1,-1,1],[0,0,1],[1,0,0]),
    }
    for name, (np_, dp, nc, dc) in cases.items():
        r = presets[name].matrix_child_to_parent
        np_ = np.asarray(np_, float); np_ /= np.linalg.norm(np_)
        dp = np.asarray(dp, float); dp /= np.linalg.norm(dp)
        nc = np.asarray(nc, float); nc /= np.linalg.norm(nc)
        dc = np.asarray(dc, float); dc /= np.linalg.norm(dc)
        assert abs(abs(float((r @ nc) @ np_)) - 1.0) < 1e-10
        assert abs(abs(float((r @ dc) @ dp)) - 1.0) < 1e-10


def test_neighbor_voting_reconstructs_all_builtin_or_families_on_synthetic_data():
    presets = orientation_relationship_presets()
    for name, p in presets.items():
        df, edge_df, _ = synthetic_parent_reconstruction_demo(name, n_per_parent=6, noise_deg=0.25, seed=13)
        gids, ori = orientations_from_dataframe(df)
        edges = edges_from_dataframe(edge_df, gids)
        ps = symmetry_group(p.parent_symmetry); cs = symmetry_group(p.child_symmetry)
        r = neighbor_voting_reconstruction(gids, ori, edges, p.matrix_child_to_parent, cs, ps)
        assert len(set(r.parent_ids.tolist())) == 2, name
        assert reconstruction_accuracy_against_labels(r, df.true_parent_id.to_numpy()) >= 0.99, name
        assert r.diagnostics["mean_fit_deg"] < 0.5, name


def test_orientation_import_can_explicitly_transpose_specimen_to_crystal_convention():
    import pandas as pd
    g = bunge_euler_to_matrix(17.0, 43.0, 91.0)
    # Build equivalent quaternion for the transposed orientation, then ask importer to transpose it back.
    from src.reconstruction import quaternion_wxyz
    qw, qx, qy, qz = quaternion_wxyz(g.T)
    df = pd.DataFrame({"grain_id": [1], "qw": [qw], "qx": [qx], "qy": [qy], "qz": [qz]})
    _, ori = orientations_from_dataframe(df, convention="specimen_to_crystal")
    assert rotation_angle_deg(ori[0] @ g.T) < 1e-8
