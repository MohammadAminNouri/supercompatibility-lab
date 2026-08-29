from __future__ import annotations

"""Independent compatibility diagnostics beyond the scalar cofactor checklist."""

from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd

Array = np.ndarray


@dataclass(frozen=True)
class RankOneCertificate:
    singular_values: Array
    normalized_rank_one_residual: float
    exact_within_tol: bool


@dataclass(frozen=True)
class TwinCompatibilityCertificate:
    chi: Array
    middle_residual: float
    bracket_margin_low: float
    bracket_margin_high: float
    compatible_within_tol: bool


@dataclass(frozen=True)
class LaminateScanResult:
    table: pd.DataFrame
    max_middle_stretch_residual: float
    rms_middle_stretch_residual: float


@dataclass(frozen=True)
class TripletOrthorhombicResult:
    alpha: float
    beta: float
    gamma: float
    tc_i_raw: float
    tc_ii_raw: float
    tc_i_normalized: float
    tc_ii_normalized: float


def hadamard_rank_one_certificate(a: Array, b: Array, tol: float = 1e-8) -> RankOneCertificate:
    """Check whether A-B is rank one using its singular spectrum.

    For a non-zero exact rank-one jump, sigma2=sigma3=0. The normalized
    residual uses sigma2/sigma1 because sigma2 is the first forbidden singular
    value. If A=B, the jump is rank zero and is accepted as trivially compatible.
    """
    d = np.asarray(a, float) - np.asarray(b, float)
    s = np.linalg.svd(d, compute_uv=False)
    if s[0] <= 1e-15:
        r = 0.0
    else:
        r = float(s[1] / s[0])
    return RankOneCertificate(s, r, bool(r <= tol))


def single_variant_ips_residual(U: Array) -> dict[str, float | bool]:
    vals = np.linalg.eigvalsh(0.5 * (np.asarray(U, float) + np.asarray(U, float).T))
    vals = np.sort(vals)
    return {
        "lambda1": float(vals[0]),
        "lambda2": float(vals[1]),
        "lambda3": float(vals[2]),
        "abs_lambda2_minus_1": float(abs(vals[1] - 1.0)),
        "bracket_lambda1_le_1": bool(vals[0] <= 1.0 + 1e-10),
        "bracket_lambda3_ge_1": bool(vals[2] >= 1.0 - 1e-10),
    }


def twin_rank_one_certificate(Ui: Array, Uj: Array, tol: float = 1e-7) -> TwinCompatibilityCertificate:
    """Ball-James/Bhattacharya rank-one compatibility test for two stretches.

    The eigenvalues chi_i of Uj^{-T} Ui^T Ui Uj^{-1} must satisfy
    chi1 <= 1, chi2 = 1, chi3 >= 1.
    """
    ui = np.asarray(Ui, float); uj = np.asarray(Uj, float)
    invj = np.linalg.inv(uj)
    b = invj.T @ ui.T @ ui @ invj
    b = 0.5 * (b + b.T)
    chi = np.sort(np.linalg.eigvalsh(b))
    low = float(1.0 - chi[0])
    high = float(chi[2] - 1.0)
    mid = float(abs(chi[1] - 1.0))
    ok = low >= -tol and high >= -tol and mid <= tol
    return TwinCompatibilityCertificate(chi, mid, low, high, bool(ok))


def pairwise_twin_compatibility(variants: Iterable[Array], tol: float = 1e-7) -> pd.DataFrame:
    vs = list(variants)
    rows = []
    for i in range(len(vs)):
        for j in range(i + 1, len(vs)):
            c = twin_rank_one_certificate(vs[i], vs[j], tol=tol)
            rows.append({
                "variant_i": i + 1,
                "variant_j": j + 1,
                "chi1": float(c.chi[0]),
                "chi2": float(c.chi[1]),
                "chi3": float(c.chi[2]),
                "abs_chi2_minus_1": c.middle_residual,
                "chi1_condition_margin": c.bracket_margin_low,
                "chi3_condition_margin": c.bracket_margin_high,
                "rank_one_compatible": c.compatible_within_tol,
            })
    return pd.DataFrame(rows)


def laminate_volume_fraction_scan(U: Array, a: Array, n: Array, points: int = 101) -> LaminateScanResult:
    """Directly test the all-volume-fraction IPS implication of a twin laminate.

    If R V = U + a⊗n is the twin equation, the average laminate gradient is
        F(f) = U + f a⊗n.
    A planar austenite/laminate interface exists exactly when the middle singular
    value of F(f) equals one. Under the full cofactor conditions this is true for
    every f in [0,1].
    """
    U = np.asarray(U, float); a = np.asarray(a, float); n = np.asarray(n, float)
    fs = np.linspace(0.0, 1.0, max(3, int(points)))
    rows = []
    residuals = []
    for f in fs:
        F = U + f * np.outer(a, n)
        sv = np.sort(np.linalg.svd(F, compute_uv=False))
        r = float(abs(sv[1] - 1.0))
        residuals.append(r)
        rows.append({
            "volume_fraction_f": float(f),
            "sigma1": float(sv[0]),
            "sigma2": float(sv[1]),
            "sigma3": float(sv[2]),
            "abs_sigma2_minus_1": r,
            "det_F": float(np.linalg.det(F)),
        })
    arr = np.asarray(residuals)
    return LaminateScanResult(pd.DataFrame(rows), float(np.max(arr)), float(np.sqrt(np.mean(arr**2))))


def orthorhombic_principal_stretches(a_parent: float, a_o: float, b_o: float, c_o: float) -> tuple[float, float, float]:
    if min(a_parent, a_o, b_o, c_o) <= 0:
        raise ValueError("All lattice lengths must be positive.")
    alpha = a_o / a_parent
    beta = b_o / (np.sqrt(2.0) * a_parent)
    gamma = c_o / (np.sqrt(2.0) * a_parent)
    return float(alpha), float(beta), float(gamma)


def triplet_condition_orthorhombic(alpha: float, beta: float, gamma: float) -> TripletOrthorhombicResult:
    """Cubic->orthorhombic triplet-condition residuals.

    This specialized representation uses positive principal stretches alpha,beta,gamma.
    TC-I:  alpha^2 gamma^2 + 2 gamma^2 beta^2 - 3 alpha^2 beta^2 = 0
    TC-II: 2 alpha^2 + beta^2 - 3 gamma^2 = 0

    The formulas are transformation-class specific and are not a B2->B19' monoclinic
    certification criterion.
    """
    if min(alpha, beta, gamma) <= 0:
        raise ValueError("Principal stretches must be positive.")
    t1 = alpha**2 * gamma**2 + 2.0 * gamma**2 * beta**2 - 3.0 * alpha**2 * beta**2
    t2 = 2.0 * alpha**2 + beta**2 - 3.0 * gamma**2
    s1 = alpha**2 * gamma**2 + 2.0 * gamma**2 * beta**2 + 3.0 * alpha**2 * beta**2
    s2 = 2.0 * alpha**2 + beta**2 + 3.0 * gamma**2
    return TripletOrthorhombicResult(
        float(alpha), float(beta), float(gamma), float(t1), float(t2),
        float(abs(t1) / max(s1, 1e-15)), float(abs(t2) / max(s2, 1e-15)),
    )
