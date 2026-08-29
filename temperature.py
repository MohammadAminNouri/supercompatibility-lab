from __future__ import annotations

import numpy as np
import pandas as pd

from .core import LatticeInput
from .distances import compatibility_dashboard


REQUIRED_COLUMNS = [
    "temperature_K",
    "a_B2_A",
    "a_B19p_A",
    "b_B19p_A",
    "c_B19p_A",
    "beta_deg",
]


def _axis_text(axis: np.ndarray) -> str:
    v = np.asarray(axis, float)
    nz = np.abs(v[np.abs(v) > 1e-8])
    if len(nz) == 0:
        return "[0 0 0]"
    w = np.rint(v / np.min(nz)).astype(int)
    for x in w:
        if x != 0:
            if x < 0:
                w *= -1
            break
    return "[" + " ".join(str(int(x)) for x in w) + "]"


def validate_temperature_frame(df: pd.DataFrame) -> pd.DataFrame:
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError("Missing required columns: " + ", ".join(missing))
    out = df[REQUIRED_COLUMNS].copy()
    for c in REQUIRED_COLUMNS:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    if out.isna().any().any():
        bad = out.index[out.isna().any(axis=1)].tolist()[:5]
        raise ValueError(f"Non-numeric or missing values in temperature data near row(s): {bad}")
    if (out["temperature_K"] <= 0).any():
        raise ValueError("Temperature must be > 0 K.")
    for c in ["a_B2_A", "a_B19p_A", "b_B19p_A", "c_B19p_A"]:
        if (out[c] <= 0).any():
            raise ValueError(f"{c} must be positive.")
    if ((out["beta_deg"] <= 0) | (out["beta_deg"] >= 180)).any():
        raise ValueError("beta_deg must lie between 0 and 180 degrees.")
    return out.sort_values("temperature_K").reset_index(drop=True)


def temperature_sweep(df: pd.DataFrame, cc1_tol: float = 1e-5, cc2_tol: float = 1e-5) -> pd.DataFrame:
    clean = validate_temperature_frame(df)
    rows: list[dict[str, object]] = []
    for r in clean.itertuples(index=False):
        inp = LatticeInput(r.a_B2_A, r.a_B19p_A, r.b_B19p_A, r.c_B19p_A, r.beta_deg)
        d = compatibility_dashboard(inp, cc1_tol=cc1_tol, cc2_tol=cc2_tol)
        rows.append({
            "temperature_K": r.temperature_K,
            "a_B2_A": r.a_B2_A,
            "a_B19p_A": r.a_B19p_A,
            "b_B19p_A": r.b_B19p_A,
            "c_B19p_A": r.c_B19p_A,
            "beta_deg": r.beta_deg,
            "lambda1": d.lambda1,
            "lambda2": d.lambda2,
            "lambda3": d.lambda3,
            "abs_lambda2_minus_1": d.lambda2_abs_residual,
            "det_U": d.det_u,
            "cmc_relative_zero": d.cmc_relative_zero,
            "cofactor_cc2_normalized": d.best_cc2_normalized,
            "cofactor_cc3_margin": d.best_cc3_margin,
            "cofactor_domain": d.best_cofactor_type,
            "cofactor_axis": _axis_text(d.best_cofactor_axis),
            "ptmc_cofactor_pass": d.ptmc_all_pass,
            "epsilon": np.nan if d.best_epsilon is None else d.best_epsilon,
            "metric_shear_match_available": d.best_epsilon is not None,
        })
    return pd.DataFrame(rows)
