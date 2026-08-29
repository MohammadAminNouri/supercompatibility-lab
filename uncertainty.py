from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .core import LatticeInput, cmc_degeneracy, cmc_matrix, normalized_metrics
from .ptmc import CofactorDomainSpec, evaluate_domain_spec, stretch_from_metrics


@dataclass(frozen=True)
class LatticeUncertainty:
    sigma_a_b2: float = 0.0
    sigma_a_b19p: float = 0.0
    sigma_b_b19p: float = 0.0
    sigma_c_b19p: float = 0.0
    sigma_beta_deg: float = 0.0


@dataclass(frozen=True)
class UncertaintySummary:
    samples: pd.DataFrame
    n_valid: int
    within_tolerance_fraction: float
    cc1_fraction: float
    cc2_fraction: float
    cc3_fraction: float
    cmc_fraction: float


def monte_carlo_uncertainty(
    inp: LatticeInput,
    sigma: LatticeUncertainty,
    domain_spec: CofactorDomainSpec,
    n: int = 2000,
    seed: int = 42,
    cc1_tol: float = 1e-3,
    cc2_tol: float = 1e-3,
    cmc_tol: float = 1e-3,
) -> UncertaintySummary:
    """Propagate independent Gaussian lattice-parameter uncertainty.

    Reported fractions are *within user-selected numerical/experimental
    tolerances*, not probabilities of satisfying exact equalities (which have
    zero probability under a continuous uncertainty model).
    """
    n = int(np.clip(n, 100, 20000))
    rng = np.random.default_rng(seed)
    means = np.array([inp.a_b2, inp.a_b19p, inp.b_b19p, inp.c_b19p, inp.beta_deg], float)
    scales = np.array([
        sigma.sigma_a_b2,
        sigma.sigma_a_b19p,
        sigma.sigma_b_b19p,
        sigma.sigma_c_b19p,
        sigma.sigma_beta_deg,
    ], float)
    if np.any(scales < 0) or not np.all(np.isfinite(scales)):
        raise ValueError("All uncertainty standard deviations must be finite and non-negative.")
    draws = rng.normal(means, scales, size=(n, 5))
    valid = (
        (draws[:, :4] > 0).all(axis=1)
        & (draws[:, 4] > 0)
        & (draws[:, 4] < 180)
    )
    draws = draws[valid]
    rows: list[dict[str, float | bool]] = []
    for a_b2, a_m, b_m, c_m, beta in draws:
        sample = LatticeInput(float(a_b2), float(a_m), float(b_m), float(c_m), float(beta))
        ma, mm = normalized_metrics(sample)
        stretch = stretch_from_metrics(ma, mm)
        try:
            cc = evaluate_domain_spec(
                stretch.U,
                domain_spec,
                cc1_tol=cc1_tol,
                cc2_tol=cc2_tol,
            )
        except ValueError:
            # If a sampled lattice leaves the symmetry/domain relation assumed
            # by the selected specification, count it as outside tolerance
            # rather than silently reclassifying the domain system.
            rows.append({
                "lambda2": float(stretch.eigenvalues[1]),
                "abs_lambda2_minus_1": abs(float(stretch.eigenvalues[1]) - 1.0),
                "cc2_normalized": np.nan,
                "cc3_margin": np.nan,
                "cmc_relative_zero": np.nan,
                "cc1_pass": False,
                "cc2_pass": False,
                "cc3_pass": False,
                "cmc_pass": False,
                "all_within_tolerance": False,
                "domain_relation_valid": False,
            })
            continue
        deg = cmc_degeneracy(cmc_matrix(ma, mm), rtol=cmc_tol)
        cmc_pass = deg.relative_zero <= cmc_tol
        all_pass = cc.cc1_pass and cc.cc2_pass and cc.cc3_pass and cmc_pass
        rows.append({
            "lambda2": cc.lambda2,
            "abs_lambda2_minus_1": abs(cc.cc1_residual),
            "cc2_normalized": cc.cc2_normalized,
            "cc3_margin": cc.cc3_margin,
            "cmc_relative_zero": deg.relative_zero,
            "cc1_pass": cc.cc1_pass,
            "cc2_pass": cc.cc2_pass,
            "cc3_pass": cc.cc3_pass,
            "cmc_pass": cmc_pass,
            "all_within_tolerance": all_pass,
            "domain_relation_valid": True,
        })
    frame = pd.DataFrame(rows)
    if frame.empty:
        raise ValueError("No physically valid samples were generated; check uncertainty magnitudes.")
    return UncertaintySummary(
        samples=frame,
        n_valid=len(frame),
        within_tolerance_fraction=float(frame["all_within_tolerance"].mean()),
        cc1_fraction=float(frame["cc1_pass"].mean()),
        cc2_fraction=float(frame["cc2_pass"].mean()),
        cc3_fraction=float(frame["cc3_pass"].mean()),
        cmc_fraction=float(frame["cmc_pass"].mean()),
    )
