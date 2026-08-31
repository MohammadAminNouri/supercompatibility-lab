from __future__ import annotations

"""Equation-level analytical engine for the B2 -> B19' compatibility workflow.

The functions in this module are intentionally explicit.  Each source-derived
calculation has a stable equation key that maps to the equation registry in
``src.provenance``.  This makes numerical outputs auditable in exported
research records and prevents a UI label from becoming the only documentation
of what was computed.
"""

from dataclasses import dataclass
from math import cos, sin, sqrt
from typing import Iterable, Literal

import numpy as np

from .core import C_A_TO_M, C_M_TO_A, Array, LatticeInput, cmc_matrix, normalized_metrics
from .symmetry import correspondence_intersection_subgroup, cubic_point_group


@dataclass(frozen=True)
class AnalyticCMCResult:
    matrix: Array
    eigenvalues_labeled: dict[str, float]
    eigenvalues_sorted: Array
    K: float
    Delta: float
    quadratic_coefficients: dict[str, float]


@dataclass(frozen=True)
class CompatibilityFamilyResult:
    name: str
    equality_residual: float
    inequality_margins: tuple[float, ...]
    equality_met: bool
    inequalities_met: bool
    met: bool
    equation_key: str
    note: str = ""


@dataclass(frozen=True)
class DegeneracyFamilyResult:
    name: str
    met: bool
    order: int
    residuals: dict[str, float]
    habit_plane: Array | None
    equation_key: str
    note: str = ""


@dataclass(frozen=True)
class IPSGeometryResult:
    U: Array
    lambdas: Array
    eigenvectors: Array
    branch: int
    v: Array
    v_prime: Array
    rotation: Array
    F: Array
    d: Array
    m: Array
    singular_values_F_minus_I: Array
    rank_one_residual: float
    plane_invariance_residual: float
    tau: float
    delta: float
    trace_identity_residual: float
    determinant_identity_residual: float


@dataclass(frozen=True)
class CorrespondenceCoset:
    index: int
    representative: Array
    elements: tuple[Array, ...]
    correspondence_matrices: tuple[Array, ...]


@dataclass(frozen=True)
class CMCSymmetryEntry:
    symmetry: Array
    proportionality: float
    zero_set_residual: float
    exact_form_residual: float


@dataclass(frozen=True)
class FTCResult:
    C_m_to_a: Array
    T_m_to_a: Array
    F_a: Array
    reconstruction_residual: float


def _beta_rad(beta_deg: float) -> float:
    if not np.isfinite(beta_deg) or not 0.0 < beta_deg < 180.0:
        raise ValueError("beta must be finite and strictly between 0° and 180°.")
    return float(np.deg2rad(beta_deg))


def niti_cmc_analytic(a: float, b: float, c: float, beta_deg: float) -> AnalyticCMCResult:
    """Closed-form normalized CMC for B2 -> B19' in the built-in correspondence.

    ``a``, ``b`` and ``c`` are the dimensionless ratios a_B19'/a_B2,
    b_B19'/a_B2 and c_B19'/a_B2.
    """
    if min(a, b, c) <= 0 or not all(np.isfinite([a, b, c])):
        raise ValueError("Normalized lattice ratios must be positive finite numbers.")
    beta = _beta_rad(beta_deg)
    cb = cos(beta)
    q = np.array(
        [
            [(b*b + c*c)/4.0 - 1.0, (b*b - c*c)/4.0, -a*c*cb/2.0],
            [(b*b - c*c)/4.0, (b*b + c*c)/4.0 - 1.0,  a*c*cb/2.0],
            [-a*c*cb/2.0,             a*c*cb/2.0,             a*a - 1.0],
        ],
        dtype=float,
    )
    K = 2.0*a*a + c*c - 4.0
    Delta = 4.0*a**4 + c**4 + 4.0*a*a*c*c*cos(2.0*beta)
    # Round-off can make a theoretically zero discriminant slightly negative.
    if Delta < -1e-12:
        raise ValueError("Analytical CMC discriminant became negative beyond round-off.")
    root = sqrt(max(0.0, Delta))
    labeled = {
        "q1": 0.5*(b*b - 2.0),
        "q2": 0.25*(K - root),
        "q3": 0.25*(K + root),
    }
    coeff = {
        # For 4 q_CMC = 0 as printed in the normalized explicit equation:
        # c²(x-y)² + b²(x+y)² + 4a²z² - 4ac(x-y)z cosβ -4(x²+y²+z²)=0.
        "x2": c*c + b*b - 4.0,
        "y2": c*c + b*b - 4.0,
        "z2": 4.0*a*a - 4.0,
        "xy": 2.0*(b*b - c*c),
        "xz": -4.0*a*c*cb,
        "yz": 4.0*a*c*cb,
    }
    return AnalyticCMCResult(q, labeled, np.sort(np.array(list(labeled.values()))), K, Delta, coeff)


def niti_cmc_from_input(inp: LatticeInput) -> AnalyticCMCResult:
    return niti_cmc_analytic(*inp.ratios(), inp.beta_deg)


def cmc_quadratic_explicit(direction: Iterable[float], a: float, b: float, c: float, beta_deg: float) -> float:
    """Return the explicit equation-44 polynomial (four times u^T CMC u)."""
    x, y, z = [float(v) for v in direction]
    beta = _beta_rad(beta_deg)
    return float(
        c*c*(x-y)**2 + b*b*(x+y)**2 + 4.0*a*a*z*z
        - 4.0*a*c*(x-y)*z*cos(beta) - 4.0*(x*x+y*y+z*z)
    )


def first_order_families(
    a: float,
    b: float,
    c: float,
    beta_deg: float,
    *,
    tol: float = 2e-6,
    c2b_interpretation: Literal["equation", "table"] = "equation",
) -> tuple[CompatibilityFamilyResult, ...]:
    """Evaluate C1, C2a, C2b and C3 with source-transparent residuals.

    The source contains an internal C2b inequality discrepancy.  The analytical
    equation/derivation and plotted C2b example require b <= sqrt(2), while a
    later distance table prints b >= sqrt(2).  Both are available explicitly;
    the default follows the analytical equation rather than silently resolving
    the discrepancy.
    """
    beta = _beta_rad(beta_deg)
    s = sin(beta)
    s2 = s*s
    rt2 = sqrt(2.0)
    out: list[CompatibilityFamilyResult] = []

    c1eq = b*b - 2.0
    c1m = 2.0*a*a + c*c - a*a*c*c*s2 - 2.0
    out.append(CompatibilityFamilyResult(
        "C1", c1eq, (c1m,), abs(c1eq) <= tol, c1m >= -tol,
        abs(c1eq) <= tol and c1m >= -tol, "CT-C1"
    ))

    denom = 1.0 - a*a*s2
    if abs(denom) < 1e-14:
        c_relation = float("nan")
    else:
        c_relation = c*c - 2.0*(1.0-a*a)/denom

    c2a_margins = (a*s - 1.0, rt2 - b)
    c2a_eq_ok = np.isfinite(c_relation) and abs(c_relation) <= tol
    c2a_ineq_ok = min(c2a_margins) >= -tol
    out.append(CompatibilityFamilyResult(
        "C2a", float(c_relation), tuple(float(x) for x in c2a_margins), bool(c2a_eq_ok), bool(c2a_ineq_ok),
        bool(c2a_eq_ok and c2a_ineq_ok), "CT-C2a"
    ))

    # Two equalities are required; a squared-form sum makes the residual
    # dimensionless and exactly comparable with the paper's distance table.
    c2b_eq = abs(a*a - 1.0) + abs(s2 - 1.0)
    if c2b_interpretation == "equation":
        margin = rt2 - b
        note = "Uses analytical equation/derivation: b ≤ √2. A later source table prints the opposite sign; see source-discrepancy record."
        key = "CT-C2b-EQ"
    elif c2b_interpretation == "table":
        margin = b - rt2
        note = "Uses the inequality printed in the later source distance table: b ≥ √2; this conflicts with the analytical equation and plotted C2b example."
        key = "CT-C2b-TABLE"
    else:
        raise ValueError("c2b_interpretation must be 'equation' or 'table'.")
    out.append(CompatibilityFamilyResult(
        "C2b", float(c2b_eq), (float(margin),), c2b_eq <= tol, margin >= -tol,
        c2b_eq <= tol and margin >= -tol, key, note
    ))

    c3_margins = (1.0-a, b-rt2)
    c3_eq_ok = np.isfinite(c_relation) and abs(c_relation) <= tol
    c3_ineq_ok = min(c3_margins) >= -tol
    out.append(CompatibilityFamilyResult(
        "C3", float(c_relation), tuple(float(x) for x in c3_margins), bool(c3_eq_ok), bool(c3_ineq_ok),
        bool(c3_eq_ok and c3_ineq_ok), "CT-C3"
    ))
    return tuple(out)


def higher_order_degeneracy_families(a: float, b: float, c: float, beta_deg: float, tol: float = 2e-6) -> tuple[DegeneracyFamilyResult, ...]:
    beta = _beta_rad(beta_deg)
    s = sin(beta)
    s2 = s*s
    rt2 = sqrt(2.0)
    denom = 1.0 - a*a*s2
    d1_relation = float("nan") if abs(denom) < 1e-14 else c*c - 2.0*(1.0-a*a)/denom
    d1_b = b*b - 2.0
    d1_branch = (a < 1.0 + tol) or (a*s >= 1.0 - tol)
    d1_met = np.isfinite(d1_relation) and abs(d1_relation) <= tol and abs(d1_b) <= tol and d1_branch
    d1_plane: Array | None = None
    if d1_met and abs(2.0-c*c) > 1e-12:
        rad = (1.0-a*a)/(2.0-c*c)
        if rad >= -1e-12:
            val = 2.0*sqrt(max(0.0, rad))
            # Appendix B gives x-y- val*z=0 for a<1 and the opposite sign
            # for the high-a branch.
            zcoef = -val if a < 1.0 else val
            d1_plane = np.array([1.0, -1.0, zcoef])

    d2_res = {"a2_minus_1": a*a-1.0, "c2_minus_2": c*c-2.0, "sin2beta_minus_1": s2-1.0}
    d2_met = max(abs(v) for v in d2_res.values()) <= tol
    d2_plane = np.array([1.0, 1.0, 0.0]) if d2_met else None

    e_res = {**d2_res, "b2_minus_2": b*b-2.0}
    e_met = max(abs(v) for v in e_res.values()) <= tol

    return (
        DegeneracyFamilyResult("D1", bool(d1_met), 2, {"c_relation": float(d1_relation), "b2_minus_2": float(d1_b)}, d1_plane, "CT-D1"),
        DegeneracyFamilyResult("D2", bool(d2_met), 2, {k: float(v) for k, v in d2_res.items()}, d2_plane, "CT-D2"),
        DegeneracyFamilyResult("E", bool(e_met), 3, {k: float(v) for k, v in e_res.items()}, None, "CT-E"),
    )


def paper_defined_distances(a: float, b: float, c: float, beta_deg: float) -> dict[str, dict[str, float | bool | str]]:
    """Source Table-2 distances, including an explicit internal-consistency audit.

    C1/C2b/frontier values reproduce the table directly from the raw Kudoh
    lattice parameters.  For the shared C2a/C3 equality, the printed expression
    is a squared-polynomial distance, but the printed numerical value 0.269706
    is reproduced by the *linear lattice-ratio* distance |c-c_target|.  Both are
    therefore returned, rather than silently choosing one.
    """
    beta = _beta_rad(beta_deg)
    s2 = sin(beta)**2
    denom = 1.0 - a*a*s2
    if abs(denom) < 1e-14:
        c_target = float("nan")
        shared_printed = float("inf")
        shared_value_reproducing = float("inf")
    else:
        target_sq = 2.0*(1.0-a*a)/denom
        c_target = sqrt(max(0.0, target_sq)) if target_sq >= 0 else float("nan")
        shared_printed = abs(c*c-target_sq)
        shared_value_reproducing = abs(c-c_target) if np.isfinite(c_target) else float("inf")
    c1_front = abs(2.0*a*a + c*c - a*a*c*c*s2 - 2.0)
    c1_ineq = 2.0*a*a + c*c - a*a*c*c*s2 >= 2.0
    c2a_ineq = a*sin(beta) >= 1.0 and b <= sqrt(2.0)
    c2a_front = abs(a*a*s2 - 1.0) + abs(b*b - 2.0)
    c2b_eq = abs(a*a - 1.0) + abs(s2 - 1.0)
    c2b_table_ineq = b >= sqrt(2.0)
    c2b_front = abs(b*b - 2.0)
    c3_ineq = a <= 1.0 and b >= sqrt(2.0)
    c3_front = abs(a*a - 1.0) + abs(b*b - 2.0)
    shared = {
        "distance_equality_printed_expression": shared_printed,
        "distance_equality_table_value_reproduction": shared_value_reproducing,
        "c_target": c_target,
        "source_consistency": "The table's printed numerical value is reproduced by |c-c_target|, not by the visually/parsed squared expression. Both are retained.",
    }
    return {
        "C1": {"distance_equality": abs(b*b-2.0), "inequality_met": c1_ineq, "distance_inequality_frontier": c1_front, "source_rule": "Table distance"},
        "C2a": {**shared, "inequality_met": c2a_ineq, "distance_inequality_frontier": c2a_front, "source_rule": "Table distance with equality-expression/value discrepancy retained"},
        "C2b": {"distance_equality": c2b_eq, "inequality_met": c2b_table_ineq, "distance_inequality_frontier": c2b_front, "source_rule": "Table distance; inequality sign conflicts with analytical C2b equation"},
        "C3": {**shared, "inequality_met": c3_ineq, "distance_inequality_frontier": c3_front, "source_rule": "Table distance with equality-expression/value discrepancy retained"},
    }


def appendix_c_o4(beta_deg: float) -> tuple[float, float]:
    """Closed-form O4 Type-I supercompatible normalized (a,c) as functions of beta."""
    beta = _beta_rad(beta_deg)
    sb, cb = sin(beta), cos(beta)
    if abs(sb) < 1e-14:
        raise ValueError("Appendix-C solution is undefined for sin(beta)=0.")
    Rrad = 1.0 - 12.0*sqrt(2.0)*cb/sb
    if Rrad < -1e-12:
        raise ValueError("No real Appendix-C O4 solution at this beta (inner radicand < 0).")
    R = sqrt(max(0.0, Rrad))
    a_rad = 6.0 - 4.0*sqrt(2.0)*cb/sb - 2.0*R
    if a_rad < -1e-12:
        raise ValueError("No real Appendix-C O4 solution at this beta (a radicand < 0).")
    root_a = sqrt(max(0.0, a_rad))
    a = 0.5*root_a
    c_rad = 9.0 + 4.0*sqrt(2.0)*cb*sb + R + cos(2.0*beta)*(3.0-R)
    if c_rad < -1e-12:
        raise ValueError("No real Appendix-C O4 solution at this beta (c radicand < 0).")
    c = 0.5*(-cb*root_a + sqrt(max(0.0, c_rad)))
    return float(a), float(c)


def o2_type_i_family(beta_deg: float) -> tuple[float, float, float]:
    beta = _beta_rad(beta_deg)
    b = sqrt(2.0)
    c = sqrt(2.0)
    a = c/sin(beta)
    return float(a), float(b), float(c)


def o2_type_ii_family(beta_deg: float) -> tuple[float, float, float]:
    beta = _beta_rad(beta_deg)
    a = 1.0
    b = sqrt(2.0)
    c = sqrt(2.0)/sin(beta)
    return float(a), float(b), float(c)


def stretch_shear_from_lambdas(lambda1: float, lambda3: float, tol: float = 1e-12) -> tuple[float, float]:
    """Exact inverse of source equations (7)-(8) for λ2=1.

    Returns non-negative shear amplitude tau and normal dilatation delta.
    """
    if lambda1 <= 0 or lambda3 <= 0:
        raise ValueError("Principal stretches must be positive.")
    delta = lambda1*lambda3 - 1.0
    tau2 = lambda1*lambda1 + lambda3*lambda3 - 1.0 - (1.0+delta)**2
    if tau2 < -tol:
        raise ValueError("The supplied stretches do not correspond to a real IPS shear/dilatation pair with λ2=1.")
    return float(sqrt(max(0.0, tau2))), float(delta)


def lambdas_from_stretch_shear(tau: float, delta: float) -> tuple[float, float]:
    if tau < 0 or not np.isfinite(tau) or not np.isfinite(delta) or 1.0+delta <= 0:
        raise ValueError("Require tau >= 0 and 1+delta > 0.")
    B = 1.0 + (1.0+delta)**2 + tau*tau
    C = (1.0+delta)**2
    disc = B*B - 4.0*C
    if disc < -1e-12:
        raise ValueError("Invalid IPS invariants.")
    root = sqrt(max(0.0, disc))
    mu1, mu3 = (B-root)/2.0, (B+root)/2.0
    return float(sqrt(max(0.0, mu1))), float(sqrt(max(0.0, mu3)))


def _rodrigues(axis: Array, angle: float) -> Array:
    k = np.asarray(axis, float)
    k /= np.linalg.norm(k)
    K = np.array([[0.0, -k[2], k[1]], [k[2], 0.0, -k[0]], [-k[1], k[0], 0.0]])
    return np.eye(3) + sin(angle)*K + (1.0-cos(angle))*(K@K)


def ips_geometry_from_stretch(U: Array, branch: Literal[-1, 1] = 1, lambda2_tol: float = 2e-6) -> IPSGeometryResult:
    """Construct source Eq. (4)-(8) IPS geometry from a symmetric stretch tensor."""
    u = 0.5*(np.asarray(U, float) + np.asarray(U, float).T)
    vals, vecs = np.linalg.eigh(u)
    l1, l2, l3 = [float(v) for v in vals]
    if abs(l2-1.0) > lambda2_tol:
        raise ValueError(f"IPS construction requires λ2≈1; got λ2={l2:.12g}.")
    if not l1 < 1.0 + lambda2_tol or not l3 > 1.0 - lambda2_tol:
        raise ValueError("IPS construction expects λ1≤1≤λ3.")
    den1 = max(1e-30, 1.0-l1*l1)
    den3 = max(1e-30, l3*l3-1.0)
    e1, e2, e3 = vecs[:, 0], vecs[:, 1], vecs[:, 2]
    v = e1/sqrt(den1) + float(branch)*e3/sqrt(den3)
    vp = u @ v
    vn = v/np.linalg.norm(v)
    vpn = vp/np.linalg.norm(vp)
    # Signed rotation around e2 mapping v' to v.
    angle = float(np.arctan2(e2 @ np.cross(vpn, vn), vpn @ vn))
    R = _rodrigues(e2, angle)
    F = R @ u
    A = F - np.eye(3)
    uu, ss, vh = np.linalg.svd(A)
    d = ss[0]*uu[:, 0]
    m = vh[0, :]
    # Fix sign deterministically.
    nz = np.flatnonzero(np.abs(m) > 1e-12)
    if len(nz) and m[nz[0]] < 0:
        m, d = -m, -d
    rank_res = float(np.linalg.norm(A - np.outer(d, m), ord="fro"))
    # Any vector in m^⊥ must be invariant; project F-I onto that plane.
    P = np.eye(3)-np.outer(m, m)
    plane_res = float(np.linalg.norm(A @ P, ord="fro"))
    delta = float(m @ d)  # det(I+d m^T)=1+m^T d
    tau_vec = d - delta*m
    tau = float(np.linalg.norm(tau_vec))
    trace_res = float((l1*l1 + l3*l3 + 1.0) - (2.0 + (1.0+delta)**2 + tau*tau))
    det_res = float((l1*l1*l3*l3) - (1.0+delta)**2)
    return IPSGeometryResult(u, vals, vecs, int(branch), v, vp, R, F, d, m, ss, rank_res, plane_res, tau, delta, trace_res, det_res)


def shear_shear_geometry(tau_perp: float, delta: float, d_norm: float | None = None) -> dict[str, float]:
    """Source Eq. (12)-(13) scalar geometry: tan(phi)=delta/tau_perp and |a|=2 cos(phi)|d|."""
    if tau_perp < 0:
        raise ValueError("tau_perp must be non-negative.")
    phi = float(np.arctan2(delta, tau_perp))
    norm = float(sqrt(tau_perp*tau_perp + delta*delta) if d_norm is None else d_norm)
    return {"phi_deg": float(np.rad2deg(phi)), "cos_phi": float(cos(phi)), "d_norm": norm, "twin_shear_norm": float(2.0*cos(phi)*norm)}


def correspondence_left_cosets() -> tuple[CorrespondenceCoset, ...]:
    """Explicit left-coset partition g_i H_C^A and corresponding C matrices."""
    G = cubic_point_group()
    H = correspondence_intersection_subgroup()

    def key(m: Array) -> tuple[int, ...]:
        return tuple(np.rint(np.asarray(m).ravel()).astype(int).tolist())

    gmap = {key(g): g for g in G}
    remaining = set(gmap)
    out: list[CorrespondenceCoset] = []
    idx = 1
    while remaining:
        k0 = min(remaining)
        rep = gmap[k0]
        elems = {key(rep @ h): rep @ h for h in H}
        keys = set(elems)
        remaining -= keys
        ordered = tuple(elems[k] for k in sorted(elems))
        cms = tuple(g @ C_A_TO_M for g in ordered)
        out.append(CorrespondenceCoset(idx, rep, ordered, cms))
        idx += 1
    return tuple(out)


def cmc_symmetry_group(cmc: Array, proportional_tol: float = 1e-10, exact_tol: float = 1e-10) -> tuple[CMCSymmetryEntry, ...]:
    """Enumerate parent symmetries preserving the CMC zero set.

    The paper defines preservation through the quadratic zero set.  Two non-zero
    quadratic forms have the same zero set when proportional in this finite
    symmetry search, so the implementation reports both a proportional zero-set
    residual and the stricter exact-form residual.
    """
    q = 0.5*(np.asarray(cmc, float)+np.asarray(cmc, float).T)
    qq = float(np.sum(q*q))
    if qq <= 1e-30:
        return tuple(CMCSymmetryEntry(g, 1.0, 0.0, 0.0) for g in cubic_point_group())
    out: list[CMCSymmetryEntry] = []
    qnorm = np.linalg.norm(q, ord="fro")
    for g in cubic_point_group():
        transformed = g.T @ q @ g
        alpha = float(np.sum(q*transformed)/qq)
        zero_res = float(np.linalg.norm(transformed-alpha*q, ord="fro")/max(qnorm, 1e-30))
        exact_res = float(np.linalg.norm(transformed-q, ord="fro")/max(qnorm, 1e-30))
        if zero_res <= proportional_tol or exact_res <= exact_tol:
            out.append(CMCSymmetryEntry(g, alpha, zero_res, exact_res))
    return tuple(out)


def ftc_from_C_and_F(C_m_to_a: Array, F_a: Array) -> FTCResult:
    """Given correspondence C^{M->A} and active distortion F_A, infer T^{M->A} from C=T F."""
    C = np.asarray(C_m_to_a, float)
    F = np.asarray(F_a, float)
    T = C @ np.linalg.inv(F)
    res = float(np.linalg.norm(C - T @ F, ord="fro"))
    return FTCResult(C, T, F, res)


def ftc_from_T_and_F(T_m_to_a: Array, F_a: Array) -> FTCResult:
    T = np.asarray(T_m_to_a, float)
    F = np.asarray(F_a, float)
    C = T @ F
    return FTCResult(C, T, F, float(np.linalg.norm(C-T@F, ord="fro")))


def verify_analytic_cmc_against_general(inp: LatticeInput) -> dict[str, float]:
    ma, mm = normalized_metrics(inp)
    general = cmc_matrix(ma, mm)
    analytic = niti_cmc_from_input(inp)
    return {
        "matrix_frobenius_residual": float(np.linalg.norm(general-analytic.matrix, ord="fro")),
        "eigenvalue_max_abs_residual": float(np.max(np.abs(np.linalg.eigvalsh(general)-analytic.eigenvalues_sorted))),
    }
