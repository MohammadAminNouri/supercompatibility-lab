"""Batch-screen B2/B19' lattice rows from a CSV without launching Streamlit."""
from __future__ import annotations

import sys

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import argparse

import pandas as pd

from src.core import LatticeInput
from src.distances import compatibility_dashboard

REQUIRED = ["a_B2_A", "a_B19p_A", "b_B19p_A", "c_B19p_A", "beta_deg"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--cc1-tol", type=float, default=1e-5)
    parser.add_argument("--cc2-tol", type=float, default=1e-5)
    args = parser.parse_args()

    df = pd.read_csv(args.input_csv)
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        raise SystemExit("Missing columns: " + ", ".join(missing))

    rows = []
    for _, r in df.iterrows():
        inp = LatticeInput(*(float(r[c]) for c in REQUIRED))
        d = compatibility_dashboard(inp, cc1_tol=args.cc1_tol, cc2_tol=args.cc2_tol)
        rows.append({
            **dict(r),
            "lambda1": d.lambda1,
            "lambda2": d.lambda2,
            "lambda3": d.lambda3,
            "abs_lambda2_minus_1": d.lambda2_abs_residual,
            "det_U": d.det_u,
            "cmc_relative_zero": d.cmc_relative_zero,
            "best_cc2_normalized": d.best_cc2_normalized,
            "best_cc3_margin": d.best_cc3_margin,
            "ptmc_cofactor_pass": d.ptmc_all_pass,
            "best_metric_epsilon": d.best_epsilon,
            "best_metric_twin": d.best_ct_twin,
        })
    pd.DataFrame(rows).to_csv(args.output_csv, index=False)


if __name__ == "__main__":
    main()
