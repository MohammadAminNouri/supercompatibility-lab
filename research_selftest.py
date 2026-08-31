from __future__ import annotations

"""Self-contained scientific audit for the flat Streamlit deployment.

This script extracts the embedded research engine and data directly from
``streamlit_app.py`` into a temporary directory, imports that extracted engine,
and runs deterministic numerical/reconstruction benchmarks.  It therefore does
not require ``src/``, ``tests/`` or ``data/`` directories to exist in GitHub.
"""

import ast
import base64
from pathlib import Path
import sys
import tempfile
import zlib

import numpy as np

ROOT = Path(__file__).resolve().parent
APP = ROOT / "streamlit_app.py"


def _literal_assignment(tree: ast.AST, name: str):
    for node in getattr(tree, "body", []):
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == name:
                    return ast.literal_eval(node.value)
    raise RuntimeError(f"Could not find literal assignment {name!r} in {APP.name}")


def extract_engine() -> Path:
    source = APP.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(APP))
    modules = _literal_assignment(tree, "_EMBEDDED_MODULES")
    data = _literal_assignment(tree, "_EMBEDDED_DATA")

    temp_root = Path(tempfile.mkdtemp(prefix="supercompatibility_final_audit_"))
    src_dir = temp_root / "src"
    data_dir = temp_root / "data"
    src_dir.mkdir(parents=True)
    data_dir.mkdir(parents=True)

    for fullname, payload in modules.items():
        raw = zlib.decompress(base64.b64decode(payload))
        if fullname == "src":
            target = src_dir / "__init__.py"
        else:
            target = src_dir / f"{fullname.rsplit('.', 1)[-1]}.py"
        target.write_bytes(raw)

    for name, payload in data.items():
        (data_dir / name).write_bytes(zlib.decompress(base64.b64decode(payload)))

    sys.path.insert(0, str(temp_root))
    return temp_root


def assert_close(actual, expected, atol, label: str):
    if not np.allclose(actual, expected, atol=atol, rtol=0.0):
        raise AssertionError(f"{label}: got {actual!r}, expected {expected!r} ± {atol}")


def main() -> None:
    extracted = extract_engine()

    from src.compatibility_methods import laminate_volume_fraction_scan, pairwise_twin_compatibility
    from src.core import LatticeInput, cmc_degeneracy, cmc_matrix, normalized_metrics, smc_matrix
    from src.distances import all_cofactor_systems, compatibility_dashboard
    from src.presets import PRESETS
    from src.ptmc import stretch_from_lattice
    from src.reconstruction import (
        edges_from_dataframe,
        grain_graph_reconstruction,
        neighbor_voting_reconstruction,
        nucleation_growth_reconstruction,
        operator_groupoid_reconstruction,
        orientation_relationship_presets,
        orientations_from_dataframe,
        symmetry_group,
        synthetic_parent_reconstruction_demo,
        unique_child_variants,
        variant_graph_reconstruction,
    )
    from src.symmetry import stretch_variants
    from src.reconstruction_academic import known_truth_validation_metrics

    # 1) Published binary NiTi metric/PTMC benchmark.
    binary = LatticeInput(3.01, 2.898, 4.108, 4.646, 97.78)
    assert_close(binary.ratios(), (0.9627906977, 1.3647840532, 1.5435215947), 1e-10, "binary ratios")
    dash = compatibility_dashboard(binary)
    assert_close(dash.lambda2, 0.9650480588, 2e-10, "binary lambda2")

    variants = stretch_variants(stretch_from_lattice(binary).U)
    if len(variants) != 12:
        raise AssertionError(f"stretch variants: got {len(variants)}, expected 12")
    pairs = pairwise_twin_compatibility(variants)
    compatible_pairs = int(pairs.rank_one_compatible.sum())
    if len(pairs) != 66 or compatible_pairs != 42:
        raise AssertionError(f"rank-one pairs: got {compatible_pairs}/{len(pairs)}, expected 42/66")

    # 2) C1-compatible benchmark: CMC planes + SMC shear vectors.
    a_b2 = 3.01
    c1 = LatticeInput(a_b2, 0.9628 * a_b2, np.sqrt(2.0) * a_b2, 1.5435 * a_b2, 97.78)
    c1_dash = compatibility_dashboard(c1)
    assert_close(c1_dash.lambda2, 1.0, 2e-12, "C1 lambda2")
    if c1_dash.best_epsilon is None:
        raise AssertionError("C1 metric intercompatibility epsilon was not computed")
    assert_close(c1_dash.best_epsilon, 0.215, 0.003, "C1 best epsilon")

    ma, mm = normalized_metrics(c1)
    deg = cmc_degeneracy(cmc_matrix(ma, mm))
    if deg.order != 1 or len(deg.habit_planes) != 2:
        raise AssertionError(f"CMC degeneracy: order={deg.order}, planes={len(deg.habit_planes)}; expected 1 and 2")
    plane_targets = [np.array([1.0, -1.0, 2.41966]), np.array([1.0, -1.0, -0.31568])]
    for target in plane_targets:
        if not any(np.allclose(p, target, atol=3e-4, rtol=0.0) for p in deg.habit_planes):
            raise AssertionError(f"Missing CMC habit plane near {target}")

    smc = smc_matrix(ma, mm)
    shear_targets = {
        2.41966: np.array([0.36938, -0.36938, -0.05378]),
        -0.31568: np.array([0.115568, -0.115568, 0.216812]),
    }
    for p in deg.habit_planes:
        key = min(shear_targets, key=lambda k: abs(p[2] - k))
        assert_close(smc @ p, shear_targets[key], 8e-5, f"SMC shear for plane z={key}")

    # 3) Full classical cofactor cross-check + direct all-volume-fraction scan.
    passing = [x for x in all_cofactor_systems(c1) if x.all_pass]
    if not passing:
        raise AssertionError("C1 benchmark has no cofactor-passing domain system")
    scan = laminate_volume_fraction_scan(stretch_from_lattice(c1).U, passing[0].a, passing[0].n, 101)
    if scan.max_middle_stretch_residual >= 1e-12:
        raise AssertionError(
            f"all-volume-fraction middle-stretch residual too large: {scan.max_middle_stretch_residual:.3e}"
        )

    # 4) Standard OR family variant-count invariants.
    expected_counts = {
        "Kurdjumov–Sachs (FCC parent → BCC child)": 24,
        "Nishiyama–Wassermann (FCC parent → BCC child)": 12,
        "Bain (FCC parent → BCC child)": 3,
        "Pitsch (FCC parent → BCC child)": 12,
        "Burgers (BCC parent → HCP child)": 12,
    }
    presets = orientation_relationship_presets()
    for name, expected in expected_counts.items():
        p = presets[name]
        ps = symmetry_group(p.parent_symmetry)
        cs = symmetry_group(p.child_symmetry)
        count = len(unique_child_variants(np.eye(3), p.matrix_child_to_parent, ps, cs))
        if count != expected:
            raise AssertionError(f"{name}: got {count} variants, expected {expected}")

    # 5) Five independent parent-reconstruction method families on deterministic synthetic data.
    ks = presets["Kurdjumov–Sachs (FCC parent → BCC child)"]
    ps = symmetry_group(ks.parent_symmetry)
    cs = symmetry_group(ks.child_symmetry)
    df, edf, _ = synthetic_parent_reconstruction_demo(n_per_parent=6, noise_deg=0.35, seed=11)
    gids, ori = orientations_from_dataframe(df)
    edges = edges_from_dataframe(edf, gids)
    methods = [
        neighbor_voting_reconstruction,
        grain_graph_reconstruction,
        variant_graph_reconstruction,
        nucleation_growth_reconstruction,
        operator_groupoid_reconstruction,
    ]
    reconstruction_rows = []
    for fn in methods:
        r = fn(gids, ori, edges, ks.matrix_child_to_parent, cs, ps)
        tv = known_truth_validation_metrics(r, df.true_parent_id.to_numpy(), edges)
        nparents = int(r.diagnostics["n_reconstructed_parents"])
        mean_fit = float(r.diagnostics["mean_fit_deg"])
        if tv["truth_ARI"] < 0.99 or tv["truth_completeness"] < 0.99 or tv["truth_boundary_F1"] < 0.99 or nparents != 2:
            raise AssertionError(
                f"{r.method}: parents={nparents}, truth_ARI={tv['truth_ARI']:.3f}, completeness={tv['truth_completeness']:.3f}, boundary_F1={tv['truth_boundary_F1']:.3f}; expected exact two-parent recovery"
            )
        reconstruction_rows.append((r.method, nparents, float(tv["truth_ARI"]), float(tv["truth_completeness"]), float(tv["truth_boundary_F1"]), mean_fit))

    print("SUPERCOMPATIBILITY LAB FINAL SCIENTIFIC SELF-TEST: PASS")
    print(f"extracted engine: {extracted}")
    print(f"binary ratios: {tuple(float(x) for x in binary.ratios())}")
    print(f"binary lambda2: {dash.lambda2:.10f}")
    print(f"stretch variants: {len(variants)}")
    print(f"rank-one compatible pairs: {compatible_pairs}/{len(pairs)}")
    print(f"C1 lambda2: {c1_dash.lambda2:.16f}")
    print(f"C1 best epsilon: {c1_dash.best_epsilon:.9f}")
    print(f"CMC habit planes: {[np.asarray(p).round(6).tolist() for p in deg.habit_planes]}")
    print(f"max all-volume-fraction |sigma2-1|: {scan.max_middle_stretch_residual:.3e}")
    for method, nparents, ari, completeness, boundary_f1, mean_fit in reconstruction_rows:
        print(f"reconstruction {method}: parents={nparents}, truth_ARI={ari:.3f}, completeness={completeness:.3f}, boundary_F1={boundary_f1:.3f}, mean_fit={mean_fit:.3f} deg")


if __name__ == "__main__":
    main()
