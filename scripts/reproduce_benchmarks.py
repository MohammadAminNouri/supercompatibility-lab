"""Print the numerical benchmark quantities used to audit the research engine."""
from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import json
from src.core import cmc_degeneracy, cmc_matrix, normalized_metrics, rank_matches, smc_matrix
from src.distances import compatibility_dashboard
from src.presets import PRESETS
from src.ptmc import stretch_from_lattice


def summarize(name: str) -> dict[str, object]:
    inp = PRESETS[name]
    if inp is None:
        raise RuntimeError(f"Preset {name!r} is not a lattice input")
    ma, mm = normalized_metrics(inp)
    cmc = cmc_matrix(ma, mm)
    deg = cmc_degeneracy(cmc)
    smc = smc_matrix(ma, mm)
    stretch = stretch_from_lattice(inp)
    dash = compatibility_dashboard(inp, cc1_tol=1e-5, cc2_tol=1e-5)
    matches = rank_matches(inp)
    return {
        "preset": name,
        "input": {
            "a_B2_A": inp.a_b2,
            "a_B19p_A": inp.a_b19p,
            "b_B19p_A": inp.b_b19p,
            "c_B19p_A": inp.c_b19p,
            "beta_deg": inp.beta_deg,
            "ratios": list(inp.ratios()),
        },
        "stretch_eigenvalues": stretch.eigenvalues.tolist(),
        "cmc_eigenvalues": deg.eigenvalues.tolist(),
        "cmc_order": deg.order,
        "habit_planes": [p.tolist() for p in deg.habit_planes],
        "smc_shears_same_plane_scaling": [(smc @ p).tolist() for p in deg.habit_planes],
        "best_metric_match": None if not matches else {
            "epsilon": matches[0].epsilon,
            "angle_deg": matches[0].angle_deg,
            "twin": matches[0].twin.label,
        },
        "dashboard": {
            "lambda2": dash.lambda2,
            "abs_lambda2_minus_1": dash.lambda2_abs_residual,
            "best_cc2_normalized": dash.best_cc2_normalized,
            "best_cc3_margin": dash.best_cc3_margin,
            "ptmc_all_pass": dash.ptmc_all_pass,
            "best_metric_epsilon": dash.best_epsilon,
        },
    }


def main() -> None:
    names = [
        "Published binary NiTi example",
        "C1-compatible teaching example",
        "Rounded supercompatible target",
    ]
    print(json.dumps([summarize(n) for n in names], indent=2))


if __name__ == "__main__":
    main()
