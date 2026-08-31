from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Literal

import numpy as np


Array = np.ndarray


# B2 -> B19' correspondence used by the model.
C_A_TO_M = np.array(
    [
        [0.0, 1.0, -1.0],
        [0.0, 1.0, 1.0],
        [1.0, 0.0, 0.0],
    ],
    dtype=float,
)
C_M_TO_A = np.linalg.inv(C_A_TO_M)


@dataclass(frozen=True)
class LatticeInput:
    a_b2: float
    a_b19p: float
    b_b19p: float
    c_b19p: float
    beta_deg: float

    def ratios(self) -> tuple[float, float, float]:
        if self.a_b2 <= 0:
            raise ValueError("B2 lattice parameter must be positive.")
        return (
            self.a_b19p / self.a_b2,
            self.b_b19p / self.a_b2,
            self.c_b19p / self.a_b2,
        )


@dataclass(frozen=True)
class DegeneracyResult:
    eigenvalues: Array
    eigenvectors: Array
    habit_planes: tuple[Array, ...]
    order: int
    exact: bool
    relative_zero: float


@dataclass(frozen=True)
class TwinResult:
    label: str
    twin_type: Literal["Type I", "Type II"]
    parent_element: str
    shear_amplitude: float
    twin_shear_vector: Array
    twin_plane_normal: Array
    martensite_element: Array


@dataclass(frozen=True)
class MatchResult:
    twin: TwinResult
    habit_plane: Array
    shear_direction: Array
    epsilon: float
    angle_deg: float


def validate_lattice(inp: LatticeInput) -> None:
    vals = [inp.a_b2, inp.a_b19p, inp.b_b19p, inp.c_b19p]
    if not all(np.isfinite(vals)) or not np.isfinite(inp.beta_deg):
        raise ValueError("All lattice parameters must be finite numbers.")
    if any(v <= 0 for v in vals):
        raise ValueError("All lattice lengths must be positive.")
    if not 0 < inp.beta_deg < 180:
        raise ValueError("Monoclinic angle beta must be between 0 and 180 degrees.")


def normalized_metrics(inp: LatticeInput) -> tuple[Array, Array]:
    """Return dimensionless A and M metrics after scaling lengths by a_B2."""
    validate_lattice(inp)
    a, b, c = inp.ratios()
    beta = np.deg2rad(inp.beta_deg)
    m_a = np.eye(3)
    m_m = np.array(
        [
            [a * a, 0.0, a * c * np.cos(beta)],
            [0.0, b * b, 0.0],
            [a * c * np.cos(beta), 0.0, c * c],
        ],
        dtype=float,
    )
    return m_a, m_m


def physical_metrics(inp: LatticeInput) -> tuple[Array, Array]:
    validate_lattice(inp)
    beta = np.deg2rad(inp.beta_deg)
    m_a = (inp.a_b2**2) * np.eye(3)
    m_m = np.array(
        [
            [inp.a_b19p**2, 0.0, inp.a_b19p * inp.c_b19p * np.cos(beta)],
            [0.0, inp.b_b19p**2, 0.0],
            [inp.a_b19p * inp.c_b19p * np.cos(beta), 0.0, inp.c_b19p**2],
        ],
        dtype=float,
    )
    return m_a, m_m


def cmc_matrix(m_a: Array, m_m: Array) -> Array:
    """Compatibility by metric correspondence matrix."""
    q = C_M_TO_A.T @ m_m @ C_M_TO_A - m_a
    return 0.5 * (q + q.T)


def smc_matrix(m_a: Array, m_m: Array) -> Array:
    """Shear by metric correspondence matrix."""
    return np.linalg.inv(m_a) - C_A_TO_M @ np.linalg.inv(m_m) @ C_A_TO_M.T


def cmc_value(direction: Iterable[float], cmc: Array) -> float:
    u = np.asarray(direction, dtype=float)
    return float(u.T @ cmc @ u)


def reciprocal_norm(plane: Iterable[float], metric: Array) -> float:
    p = np.asarray(plane, dtype=float)
    return float(np.sqrt(max(0.0, p.T @ np.linalg.inv(metric) @ p)))


def direct_norm(direction: Iterable[float], metric: Array) -> float:
    u = np.asarray(direction, dtype=float)
    return float(np.sqrt(max(0.0, u.T @ metric @ u)))


def normalize_plane(plane: Iterable[float], metric: Array) -> Array:
    p = np.asarray(plane, dtype=float)
    n = reciprocal_norm(p, metric)
    if n <= 1e-15:
        raise ValueError("Plane vector cannot be zero.")
    return p / n


def normalize_direction(direction: Iterable[float], metric: Array) -> Array:
    u = np.asarray(direction, dtype=float)
    n = direct_norm(u, metric)
    if n <= 1e-15:
        raise ValueError("Direction vector cannot be zero.")
    return u / n


def canonical_plane(plane: Iterable[float], atol: float = 1e-12) -> Array:
    """Scale a plane vector so its first non-zero component is +1."""
    p = np.asarray(plane, dtype=float).copy()
    nz = np.flatnonzero(np.abs(p) > atol)
    if len(nz) == 0:
        return p
    p /= p[nz[0]]
    if p[nz[0]] < 0:
        p *= -1.0
    p[np.abs(p) < atol] = 0.0
    return p


def _first_order_planes(evals: Array, evecs: Array, zero_idx: int) -> tuple[Array, Array]:
    other = [i for i in range(3) if i != zero_idx]
    i, j = other
    if evals[i] * evals[j] >= 0:
        raise ValueError("First-order plane factorization requires opposite-sign non-zero eigenvalues.")
    pos = i if evals[i] > 0 else j
    neg = j if evals[i] > 0 else i
    planes: list[Array] = []
    for sign in (1.0, -1.0):
        coeff = np.zeros(3)
        coeff[neg] = np.sqrt(-evals[neg])
        coeff[pos] = sign * np.sqrt(evals[pos])
        p = evecs @ coeff
        planes.append(canonical_plane(p))
    return planes[0], planes[1]


def cmc_degeneracy(cmc: Array, rtol: float = 2e-6, atol: float = 1e-10) -> DegeneracyResult:
    evals, evecs = np.linalg.eigh(cmc)
    # Scale all zero tests to the spectral magnitude so the degeneracy diagnosis
    # is invariant to a consistent change of length units in a metric matrix.
    scale = max(float(np.max(np.abs(evals))), np.finfo(float).tiny)
    zero_mask = np.abs(evals) <= max(atol, rtol) * scale
    n_zero = int(np.sum(zero_mask))
    relative_zero = float(np.min(np.abs(evals)) / scale)

    if n_zero == 3:
        return DegeneracyResult(evals, evecs, tuple(), 3, True, relative_zero)

    if n_zero == 2:
        nonzero = int(np.flatnonzero(~zero_mask)[0])
        # q = lambda*(v.u)^2 = 0 -> one plane normal to that eigenvector.
        plane = canonical_plane(evecs[:, nonzero])
        return DegeneracyResult(evals, evecs, (plane,), 2, True, relative_zero)

    if n_zero == 1:
        zero_idx = int(np.flatnonzero(zero_mask)[0])
        other = [i for i in range(3) if i != zero_idx]
        if evals[other[0]] * evals[other[1]] <= 0:
            planes = _first_order_planes(evals, evecs, zero_idx)
            # Sort for stable display: higher third component first for the common B2/B19' case.
            planes = tuple(sorted(planes, key=lambda p: tuple(-np.round(p, 12))))
            return DegeneracyResult(evals, evecs, planes, 1, True, relative_zero)

    return DegeneracyResult(evals, evecs, tuple(), 0, False, relative_zero)


def shear_direction_from_plane(plane: Iterable[float], smc: Array, normalize: bool = False, metric: Array | None = None) -> Array:
    p = np.asarray(plane, dtype=float)
    if normalize:
        if metric is None:
            raise ValueError("Metric is required when normalize=True.")
        p = normalize_plane(p, metric)
    return smc @ p


def c1_status(inp: LatticeInput) -> dict[str, float | bool]:
    a, b, c = inp.ratios()
    beta = np.deg2rad(inp.beta_deg)
    equality_residual = b * b - 2.0
    inequality_value = 2.0 * a * a + c * c - (a * a * c * c * np.sin(beta) ** 2) - 2.0
    return {
        "b_ratio": b,
        "target_b_ratio": float(np.sqrt(2.0)),
        "equality_residual": float(equality_residual),
        "inequality_margin": float(inequality_value),
        "equality_met": bool(abs(equality_residual) <= 2e-6),
        "inequality_met": bool(inequality_value >= -2e-6),
    }


def compatibility_conditions(inp: LatticeInput, tol: float = 2e-5) -> list[dict[str, object]]:
    """Evaluate C1, C2a, C2b and C3 in normalized B2/B19' notation."""
    a, b, c = inp.ratios()
    beta = np.deg2rad(inp.beta_deg)
    s = np.sin(beta)
    sqrt2 = np.sqrt(2.0)
    out: list[dict[str, object]] = []

    c1_eq = b - sqrt2
    c1_ineq = 2 * a * a + c * c - a * a * c * c * s * s - 2
    out.append({"name": "C1", "equality_residual": float(c1_eq), "inequality_margin": float(c1_ineq), "met": abs(c1_eq) <= tol and c1_ineq >= -tol})

    denom = 1 - a * a * s * s
    rad = (1 - a * a) / denom if abs(denom) > 1e-14 else np.nan
    c_target = sqrt2 * np.sqrt(rad) if np.isfinite(rad) and rad >= 0 else np.nan
    c2a_eq = c - c_target if np.isfinite(c_target) else np.nan
    c2a_margin = min(a * s - 1.0, sqrt2 - b)
    out.append({"name": "C2a", "equality_residual": float(c2a_eq) if np.isfinite(c2a_eq) else np.nan, "inequality_margin": float(c2a_margin), "met": np.isfinite(c2a_eq) and abs(c2a_eq) <= tol and c2a_margin >= -tol})

    c2b_eq = max(abs(a - 1.0), abs(inp.beta_deg - 90.0) / 90.0)
    c2b_margin = sqrt2 - b
    out.append({"name": "C2b", "equality_residual": float(c2b_eq), "inequality_margin": float(c2b_margin), "met": abs(a - 1.0) <= tol and abs(inp.beta_deg - 90.0) <= 90 * tol and c2b_margin >= -tol})

    c3_eq = c - c_target if np.isfinite(c_target) else np.nan
    c3_margin = min(1.0 - a, b - sqrt2)
    out.append({"name": "C3", "equality_residual": float(c3_eq) if np.isfinite(c3_eq) else np.nan, "inequality_margin": float(c3_margin), "met": np.isfinite(c3_eq) and abs(c3_eq) <= tol and c3_margin >= -tol})
    return out


def _reflection_matrix(plane: Iterable[float]) -> Array:
    n = np.asarray(plane, dtype=float)
    n /= np.linalg.norm(n)
    return np.eye(3) - 2.0 * np.outer(n, n)


def _rotation_pi_matrix(axis: Iterable[float]) -> Array:
    u = np.asarray(axis, dtype=float)
    u /= np.linalg.norm(u)
    return 2.0 * np.outer(u, u) - np.eye(3)


def type_i_twin(
    m_a: Array,
    m_m: Array,
    parent_plane: Iterable[float],
    label: str,
    parent_element: str,
) -> TwinResult:
    p_a = np.asarray(parent_plane, dtype=float)
    g = _reflection_matrix(p_a)
    c_int = C_M_TO_A @ g @ C_A_TO_M
    p_m = np.linalg.inv(C_M_TO_A).T @ p_a
    n_m = np.linalg.inv(m_m) @ p_m
    s2 = float(np.trace(c_int.T @ m_m @ c_int @ np.linalg.inv(m_m)) - 3.0)
    if s2 < -1e-8:
        raise ValueError("Negative squared twin shear encountered.")
    shear = float(np.sqrt(max(0.0, s2)))
    a_m = -(c_int + np.eye(3)) @ n_m
    a_a_dir = C_A_TO_M @ a_m
    if direct_norm(a_a_dir, m_a) <= 1e-14 or shear <= 1e-14:
        twin_vec = np.zeros(3)
    else:
        twin_vec = shear * normalize_direction(a_a_dir, m_a)
    n_a = np.linalg.inv(m_a) @ p_a
    n_a = normalize_direction(n_a, m_a)
    return TwinResult(label, "Type I", parent_element, shear, twin_vec, n_a, p_m)


def type_ii_twin(
    m_a: Array,
    m_m: Array,
    parent_axis: Iterable[float],
    label: str,
    parent_element: str,
) -> TwinResult:
    a_a = np.asarray(parent_axis, dtype=float)
    g = _rotation_pi_matrix(a_a)
    c_int = C_M_TO_A @ g @ C_A_TO_M
    a_m = C_M_TO_A @ a_a
    p_m = m_m @ a_m
    c_star = np.linalg.inv(c_int).T
    jp_m = -(c_star - np.eye(3)) @ p_m
    # Reciprocal-space coordinate mapping back to the austenite basis.
    jp_a = C_M_TO_A.T @ jp_m
    s2 = float(np.trace(c_int @ np.linalg.inv(m_m) @ c_int.T @ m_m) - 3.0)
    if s2 < -1e-8:
        raise ValueError("Negative squared twin shear encountered.")
    shear = float(np.sqrt(max(0.0, s2)))
    twin_vec = np.zeros(3) if shear <= 1e-14 else shear * normalize_direction(a_a, m_a)
    n_a = np.linalg.inv(m_a) @ jp_a
    if direct_norm(n_a, m_a) <= 1e-14:
        n_a = np.zeros(3)
    else:
        n_a = normalize_direction(n_a, m_a)
    return TwinResult(label, "Type II", parent_element, shear, twin_vec, n_a, jp_m)


def validated_twin_candidates(m_a: Array, m_m: Array) -> list[TwinResult]:
    """A compact set of B2 two-fold symmetry candidates validated against benchmark values."""
    candidates = [
        type_i_twin(m_a, m_m, (0, 0, 1), "Compound / Type-I equivalent", "parent mirror (001)"),
        type_ii_twin(m_a, m_m, (0, 0, 1), "Compound / Type-II equivalent", "parent 180° axis [001]"),
        type_i_twin(m_a, m_m, (1, 0, 1), "Type I family A", "parent mirror (101)"),
        type_ii_twin(m_a, m_m, (1, 0, 1), "Type II family A", "parent 180° axis [101]"),
        type_i_twin(m_a, m_m, (1, 0, 0), "Type I family B", "parent mirror (100)"),
        type_ii_twin(m_a, m_m, (1, 0, 0), "Type II family B", "parent 180° axis [100]"),
    ]
    return [c for c in candidates if c.shear_amplitude > 1e-12 and direct_norm(c.twin_plane_normal, m_a) > 1e-12]


def shear_match(
    habit_plane: Iterable[float],
    smc: Array,
    twin: TwinResult,
    m_a: Array,
) -> MatchResult:
    m_unit = normalize_plane(habit_plane, m_a)
    d = smc @ m_unit
    if twin.shear_amplitude <= 1e-14:
        return MatchResult(twin, np.asarray(habit_plane, float), d, float("inf"), float("nan"))
    mismatch = 2.0 * float(m_unit.T @ twin.twin_plane_normal) * d - twin.twin_shear_vector
    epsilon = direct_norm(mismatch, m_a) / twin.shear_amplitude
    dn = direct_norm(d, m_a)
    an = direct_norm(twin.twin_shear_vector, m_a)
    if dn <= 1e-14 or an <= 1e-14:
        angle = float("nan")
    else:
        cosang = float((d.T @ m_a @ twin.twin_shear_vector) / (dn * an))
        angle = float(np.rad2deg(np.arccos(np.clip(cosang, -1.0, 1.0))))
        angle = min(angle, 180.0 - angle)
    return MatchResult(twin, np.asarray(habit_plane, float), d, float(epsilon), angle)


def rank_matches(inp: LatticeInput) -> list[MatchResult]:
    m_a, m_m = normalized_metrics(inp)
    cmc = cmc_matrix(m_a, m_m)
    deg = cmc_degeneracy(cmc)
    if not deg.habit_planes:
        return []
    smc = smc_matrix(m_a, m_m)
    results: list[MatchResult] = []
    for twin in validated_twin_candidates(m_a, m_m):
        for plane in deg.habit_planes:
            results.append(shear_match(plane, smc, twin, m_a))
    return sorted(results, key=lambda r: r.epsilon)


def format_vector(v: Iterable[float], digits: int = 5) -> str:
    a = np.asarray(v, dtype=float)
    return "[" + ", ".join(f"{x:.{digits}f}" for x in a) + "]"
