from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.multioutput import MultiOutputRegressor
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .core import LatticeInput
from .distances import compatibility_dashboard


TARGET_COLUMNS = ["a_B2_A", "a_B19p_A", "b_B19p_A", "c_B19p_A", "beta_deg"]


@dataclass
class MLResult:
    model: Any
    model_name: str
    feature_columns: list[str]
    target_columns: list[str]
    cv_metrics: pd.DataFrame
    n_samples: int


def _model(name: str, seed: int):
    if name == "Linear regression":
        return Pipeline([
            ("scale", StandardScaler()),
            ("reg", LinearRegression()),
        ])
    if name == "Random forest":
        return RandomForestRegressor(n_estimators=350, min_samples_leaf=2, random_state=seed, n_jobs=-1)
    if name == "Extra trees":
        return ExtraTreesRegressor(n_estimators=350, min_samples_leaf=2, random_state=seed, n_jobs=-1)
    raise ValueError(f"Unknown model: {name}")


def validate_ml_frame(df: pd.DataFrame, feature_columns: list[str]) -> pd.DataFrame:
    required = feature_columns + TARGET_COLUMNS
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError("Missing required columns: " + ", ".join(missing))
    out = df[required].copy()
    for c in required:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    out = out.dropna()
    if len(out) < 8:
        raise ValueError("At least 8 complete rows are required; 20+ is strongly recommended for any meaningful validation.")
    if not feature_columns:
        raise ValueError("Choose at least one numeric composition/processing feature.")
    return out


def train_lattice_model(
    df: pd.DataFrame,
    feature_columns: list[str],
    model_name: str = "Extra trees",
    folds: int = 5,
    seed: int = 42,
) -> MLResult:
    clean = validate_ml_frame(df, feature_columns)
    X = clean[feature_columns].to_numpy(float)
    y = clean[TARGET_COLUMNS].to_numpy(float)
    n_splits = int(np.clip(folds, 2, min(10, len(clean))))
    cv = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
    estimator = _model(model_name, seed)
    pred = cross_val_predict(estimator, X, y, cv=cv, n_jobs=1)
    rows = []
    for j, target in enumerate(TARGET_COLUMNS):
        rows.append({
            "target": target,
            "MAE": mean_absolute_error(y[:, j], pred[:, j]),
            "RMSE": float(np.sqrt(mean_squared_error(y[:, j], pred[:, j]))),
            "R2": r2_score(y[:, j], pred[:, j]),
        })
    estimator.fit(X, y)
    return MLResult(estimator, model_name, feature_columns, TARGET_COLUMNS.copy(), pd.DataFrame(rows), len(clean))


def predict_lattice(result: MLResult, candidates: pd.DataFrame) -> pd.DataFrame:
    missing = [c for c in result.feature_columns if c not in candidates.columns]
    if missing:
        raise ValueError("Candidate file is missing feature columns: " + ", ".join(missing))
    x = candidates[result.feature_columns].apply(pd.to_numeric, errors="coerce")
    if x.isna().any().any():
        raise ValueError("Candidate feature columns contain missing/non-numeric values.")
    pred = result.model.predict(x.to_numpy(float))
    out = candidates.copy().reset_index(drop=True)
    for j, target in enumerate(result.target_columns):
        out[f"pred_{target}"] = pred[:, j]
    return out


def physics_screen_predictions(predicted: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for idx, r in predicted.iterrows():
        vals = [r[f"pred_{c}"] for c in TARGET_COLUMNS]
        if not all(np.isfinite(vals)) or min(vals[:4]) <= 0 or not (0 < vals[4] < 180):
            row = dict(r)
            row.update({"physics_valid": False})
            rows.append(row)
            continue
        inp = LatticeInput(*[float(v) for v in vals])
        d = compatibility_dashboard(inp, cc1_tol=1e-3, cc2_tol=1e-3)
        row = dict(r)
        row.update({
            "physics_valid": True,
            "abs_lambda2_minus_1": d.lambda2_abs_residual,
            "cmc_relative_zero": d.cmc_relative_zero,
            "best_cc2_normalized": d.best_cc2_normalized,
            "best_cc3_margin": d.best_cc3_margin,
            "ptmc_cofactor_within_1e-3": d.ptmc_all_pass,
            "epsilon_if_exact_cmc": np.nan if d.best_epsilon is None else d.best_epsilon,
        })
        rows.append(row)
    out = pd.DataFrame(rows)
    if "physics_valid" in out:
        out = out.sort_values(
            ["physics_valid", "abs_lambda2_minus_1", "best_cc2_normalized", "cmc_relative_zero"],
            ascending=[False, True, True, True],
            na_position="last",
        ).reset_index(drop=True)
    return out
