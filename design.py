from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution
from scipy.stats import qmc

from .core import LatticeInput, cmc_degeneracy, cmc_matrix, normalized_metrics
from .ptmc import CofactorDomainSpec, evaluate_domain_spec, stretch_from_metrics


@dataclass(frozen=True)
class DesignWeights:
    lambda2: float = 1.0
    cc2: float = 1.0
    cc3_penalty: float = 1.0
    cmc: float = 1.0
    proximity: float = 0.05


@dataclass(frozen=True)
class InverseDesignResult:
    input: LatticeInput
    objective: float
    lambda2_residual: float
    cc2_normalized: float
    cc3_margin: float
    cmc_relative_zero: float
    axis: np.ndarray
    domain_type: str
    partner_axis: np.ndarray | None


def _score(
    x: np.ndarray,
    a_b2: float,
    domain_spec: CofactorDomainSpec,
    weights: DesignWeights,
    reference: np.ndarray | None,
) -> tuple[float, dict[str, float]]:
    a, b, c, beta = [float(z) for z in x]
    if min(a, b, c) <= 0 or not (1.0 < beta < 179.0):
        return 1e6, {}
    inp = LatticeInput(a_b2, a * a_b2, b * a_b2, c * a_b2, beta)
    ma, mm = normalized_metrics(inp)
    stretch = stretch_from_metrics(ma, mm)
    try:
        cc = evaluate_domain_spec(stretch.U, domain_spec)
    except ValueError:
        return 100.0, {"lambda2_residual": abs(float(stretch.eigenvalues[1]) - 1.0), "cc2_normalized": 1.0, "cc3_margin": -1.0, "cmc_relative_zero": 1.0}
    deg = cmc_degeneracy(cmc_matrix(ma, mm), rtol=1e-8)
    l2 = abs(cc.cc1_residual)
    cc2 = cc.cc2_normalized
    cc3p = max(0.0, -cc.cc3_margin)
    cmc = deg.relative_zero
    proximity = 0.0
    if reference is not None:
        # Dimensionless scaled distance; beta scaled by 10 degrees.
        delta = np.array([a-reference[0], b-reference[1], c-reference[2], (beta-reference[3])/10.0])
        proximity = float(np.linalg.norm(delta))
    value = (
        weights.lambda2 * l2
        + weights.cc2 * cc2
        + weights.cc3_penalty * cc3p
        + weights.cmc * cmc
        + weights.proximity * proximity
    )
    return float(value), {
        "lambda2_residual": l2,
        "cc2_normalized": cc2,
        "cc3_margin": cc.cc3_margin,
        "cmc_relative_zero": cmc,
    }


def inverse_lattice_design(
    a_b2: float,
    domain_spec: CofactorDomainSpec,
    bounds: dict[str, tuple[float, float]],
    weights: DesignWeights = DesignWeights(),
    reference: tuple[float, float, float, float] | None = None,
    seed: int = 42,
    maxiter: int = 80,
) -> InverseDesignResult:
    ordered_bounds = [bounds["a"], bounds["b"], bounds["c"], bounds["beta_deg"]]
    ref = None if reference is None else np.asarray(reference, float)

    def objective(x: np.ndarray) -> float:
        return _score(x, a_b2, domain_spec, weights, ref)[0]

    opt = differential_evolution(
        objective,
        ordered_bounds,
        seed=seed,
        maxiter=int(np.clip(maxiter, 10, 300)),
        popsize=12,
        tol=1e-8,
        polish=True,
        workers=1,
        updating="immediate",
    )
    a, b, c, beta = [float(v) for v in opt.x]
    inp = LatticeInput(a_b2, a*a_b2, b*a_b2, c*a_b2, beta)
    score, metrics = _score(opt.x, a_b2, domain_spec, weights, ref)
    return InverseDesignResult(
        input=inp,
        objective=score,
        lambda2_residual=metrics["lambda2_residual"],
        cc2_normalized=metrics["cc2_normalized"],
        cc3_margin=metrics["cc3_margin"],
        cmc_relative_zero=metrics["cmc_relative_zero"],
        axis=np.asarray(domain_spec.axis, float),
        domain_type=domain_spec.domain_type,
        partner_axis=None if domain_spec.partner_axis is None else np.asarray(domain_spec.partner_axis, float),
    )


def _pareto_mask(values: np.ndarray) -> np.ndarray:
    n = len(values)
    keep = np.ones(n, dtype=bool)
    for i in range(n):
        if not keep[i]:
            continue
        # A point is dominated if another point is <= in every objective and < in at least one.
        dominated = np.all(values <= values[i], axis=1) & np.any(values < values[i], axis=1)
        dominated[i] = False
        if np.any(dominated):
            keep[i] = False
    return keep


def pareto_lattice_scan(
    a_b2: float,
    domain_spec: CofactorDomainSpec,
    bounds: dict[str, tuple[float, float]],
    n_samples: int = 800,
    seed: int = 42,
) -> pd.DataFrame:
    """Latin-hypercube scan with a Pareto flag for multi-objective exploration."""
    n_samples = int(np.clip(n_samples, 100, 5000))
    sampler = qmc.LatinHypercube(d=4, seed=seed)
    sample = sampler.random(n_samples)
    lo = np.array([bounds[k][0] for k in ("a", "b", "c", "beta_deg")], float)
    hi = np.array([bounds[k][1] for k in ("a", "b", "c", "beta_deg")], float)
    x = qmc.scale(sample, lo, hi)
    rows: list[dict[str, float]] = []
    for a, b, c, beta in x:
        inp = LatticeInput(a_b2, a*a_b2, b*a_b2, c*a_b2, beta)
        ma, mm = normalized_metrics(inp)
        stretch = stretch_from_metrics(ma, mm)
        try:
            cc = evaluate_domain_spec(stretch.U, domain_spec)
        except ValueError:
            continue
        deg = cmc_degeneracy(cmc_matrix(ma, mm), rtol=1e-8)
        rows.append({
            "a": a,
            "b": b,
            "c": c,
            "beta_deg": beta,
            "abs_lambda2_minus_1": abs(cc.cc1_residual),
            "cc2_normalized": cc.cc2_normalized,
            "cc3_penalty": max(0.0, -cc.cc3_margin),
            "cmc_relative_zero": deg.relative_zero,
        })
    df = pd.DataFrame(rows)
    vals = df[["abs_lambda2_minus_1", "cc2_normalized", "cc3_penalty", "cmc_relative_zero"]].to_numpy()
    # Work in log-scaled residuals for dominance stability while preserving order.
    vals = np.log10(vals + 1e-14)
    df["pareto"] = _pareto_mask(vals)
    return df.sort_values(["pareto", "abs_lambda2_minus_1", "cc2_normalized"], ascending=[False, True, True]).reset_index(drop=True)
