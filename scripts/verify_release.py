"""Deterministic release audit for numerical and reconstruction benchmarks."""
from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np

from src.compatibility_methods import laminate_volume_fraction_scan, pairwise_twin_compatibility
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


def main() -> None:
    inp = PRESETS["Published binary NiTi example"]
    c1 = PRESETS["C1-compatible teaching example"]
    assert inp is not None and c1 is not None

    dash = compatibility_dashboard(inp)
    variants = stretch_variants(stretch_from_lattice(inp).U)
    pair = pairwise_twin_compatibility(variants)
    print("Binary benchmark")
    print(f"  ratios = {tuple(round(x, 7) for x in inp.ratios())}")
    print(f"  lambda2 = {dash.lambda2:.10f}")
    print(f"  stretch variants = {len(variants)}")
    print(f"  rank-one-compatible variant pairs = {int(pair.rank_one_compatible.sum())}/{len(pair)}")

    cc = [x for x in all_cofactor_systems(c1) if x.all_pass]
    assert cc, "C1 benchmark should contain at least one cofactor-passing domain system."
    scan = laminate_volume_fraction_scan(stretch_from_lattice(c1).U, cc[0].a, cc[0].n, 101)
    print("C1 all-volume-fraction audit")
    print(f"  selected domain = {cc[0].domain_type}")
    print(f"  max |sigma2-1| = {scan.max_middle_stretch_residual:.3e}")
    assert scan.max_middle_stretch_residual < 1e-12

    p = orientation_relationship_presets()["Kurdjumov–Sachs (FCC parent → BCC child)"]
    ps = symmetry_group(p.parent_symmetry); cs = symmetry_group(p.child_symmetry)
    assert len(unique_child_variants(np.eye(3), p.matrix_child_to_parent, ps, cs)) == 24
    df, edf, _ = synthetic_parent_reconstruction_demo(n_per_parent=6, noise_deg=0.35, seed=11)
    gids, ori = orientations_from_dataframe(df); edges = edges_from_dataframe(edf, gids)
    funcs = [
        neighbor_voting_reconstruction,
        grain_graph_reconstruction,
        variant_graph_reconstruction,
        nucleation_growth_reconstruction,
        operator_groupoid_reconstruction,
    ]
    print("Synthetic parent-reconstruction audit")
    for fn in funcs:
        r = fn(gids, ori, edges, p.matrix_child_to_parent, cs, ps)
        tv = known_truth_validation_metrics(r, df.true_parent_id.to_numpy(), edges)
        print(f"  {r.method}: parents={r.diagnostics['n_reconstructed_parents']}, truth_ARI={tv['truth_ARI']:.3f}, completeness={tv['truth_completeness']:.3f}, boundary_F1={tv['truth_boundary_F1']:.3f}, mean_fit={r.diagnostics['mean_fit_deg']:.3f} deg")
        assert tv["truth_ARI"] >= 0.99
        assert tv["truth_completeness"] >= 0.99
        assert tv["truth_boundary_F1"] >= 0.99
        assert int(r.diagnostics["n_reconstructed_parents"]) == 2

    print("Release audit: PASS")


if __name__ == "__main__":
    main()
