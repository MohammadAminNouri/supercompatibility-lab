from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Literal

import numpy as np

from .core import C_M_TO_A, Array, LatticeInput, normalized_metrics


DomainType = Literal["Type I", "Type II"]
DomainLabel = Literal["Type I", "Type II", "Compound 1", "Compound 2"]


@dataclass(frozen=True)
class StretchResult:
    U: Array
    U2: Array
    eigenvalues: Array
    eigenvectors: Array
    determinant: float


@dataclass(frozen=True)
class CofactorResult:
    domain_type: DomainLabel
    axis: Array
    a: Array
    n: Array
    lambda1: float
    lambda2: float
    lambda3: float
    cc1_residual: float
    cc2_raw: float
    cc2_normalized: float
    cc3_margin: float
    cc2_simplified_residual: float
    cc1_pass: bool
    cc2_pass: bool
    cc3_pass: bool
    partner_axis: Array | None = None

    @property
    def all_pass(self) -> bool:
        return self.cc1_pass and self.cc2_pass and self.cc3_pass


@dataclass(frozen=True)
class CofactorDomainSpec:
    domain_type: DomainLabel
    axis: Array
    partner_axis: Array | None = None
    compound_solution: int | None = None


def symmetric_sqrt(matrix: Array) -> Array:
    """Principal square-root of a real symmetric positive-definite matrix."""
    a = 0.5 * (np.asarray(matrix, float) + np.asarray(matrix, float).T)
    vals, vecs = np.linalg.eigh(a)
    if np.min(vals) <= 0:
        raise ValueError("Stretch construction requires a positive-definite metric matrix.")
    return (vecs * np.sqrt(vals)) @ vecs.T


def cofactor_matrix(matrix: Array) -> Array:
    """Matrix of signed 2x2 minors; remains valid for singular 3x3 matrices."""
    m = np.asarray(matrix, float)
    if m.shape != (3, 3):
        raise ValueError("cofactor_matrix expects a 3x3 matrix.")
    out = np.empty((3, 3), dtype=float)
    for i in range(3):
        for j in range(3):
            minor = np.delete(np.delete(m, i, axis=0), j, axis=1)
            out[i, j] = ((-1) ** (i + j)) * np.linalg.det(minor)
    return out


def stretch_from_metrics(m_a: Array, m_m: Array) -> StretchResult:
    """Build the transformation stretch tensor from metrics and correspondence.

    It solves U^T M_A U = C^T M_M C by mapping through the symmetric square-root
    of M_A. For the built-in cubic B2 parent this reduces exactly to
    U^2 = C^T M_M C after normalization.
    """
    ma = 0.5 * (np.asarray(m_a, float) + np.asarray(m_a, float).T)
    mm = 0.5 * (np.asarray(m_m, float) + np.asarray(m_m, float).T)
    a_half = symmetric_sqrt(ma)
    a_inv_half = np.linalg.inv(a_half)
    target = C_M_TO_A.T @ mm @ C_M_TO_A
    # In orthonormalized parent coordinates, V^T V = A^{-1/2} target A^{-1/2}.
    gram = a_inv_half @ target @ a_inv_half
    v = symmetric_sqrt(gram)
    # U below is represented in the parent orthonormalized frame. This is the
    # correct frame for the classical PTMC/cofactor equations.
    U = 0.5 * (v + v.T)
    vals, vecs = np.linalg.eigh(U)
    return StretchResult(
        U=U,
        U2=U @ U,
        eigenvalues=vals,
        eigenvectors=vecs,
        determinant=float(np.linalg.det(U)),
    )


def stretch_from_lattice(inp: LatticeInput) -> StretchResult:
    ma, mm = normalized_metrics(inp)
    return stretch_from_metrics(ma, mm)


def _unit(v: Iterable[float]) -> Array:
    x = np.asarray(v, dtype=float)
    n = float(np.linalg.norm(x))
    if n <= 1e-15:
        raise ValueError("Axis must be non-zero.")
    return x / n


def twofold_related_stretch(U: Array, axis: Iterable[float]) -> Array:
    """Return U-hat = Q U Q^T with Q the 180° rotation around `axis`."""
    e = _unit(axis)
    q = -np.eye(3) + 2.0 * np.outer(e, e)
    return q @ U @ q.T


def domain_rank_one(U: Array, axis: Iterable[float], domain_type: DomainType) -> tuple[Array, Array]:
    """Classical Type-I/Type-II rank-one domain solution.

    For U-hat = Q U Q^T, Q = -I + 2 e⊗e, returns vectors (a,n) in
    R U-hat = U + a⊗n. These formulae are exact in geometrically nonlinear
    theory; the scaling of a and n is the conventional one used in the
    cofactor-condition literature.
    """
    e = _unit(axis)
    u = np.asarray(U, dtype=float)
    if domain_type == "Type I":
        uinv_e = np.linalg.solve(u, e)
        denom = float(uinv_e @ uinv_e)
        n = e
        a = 2.0 * (uinv_e / denom - u @ e)
    elif domain_type == "Type II":
        ue = u @ e
        denom = float(ue @ ue)
        a = ue
        n = 2.0 * (e - (u @ u @ e) / denom)
    else:
        raise ValueError(f"Unknown domain type: {domain_type}")
    return a, n


def _evaluate_cofactor_rank_one(
    U: Array,
    axis: Iterable[float],
    domain_label: DomainLabel,
    a: Array,
    n: Array,
    cc1_tol: float,
    cc2_tol: float,
    cc3_tol: float,
    simplified_residual: float = float("nan"),
    partner_axis: Array | None = None,
) -> CofactorResult:
    u = np.asarray(U, float)
    vals, _ = np.linalg.eigh(u)
    l1, l2, l3 = [float(v) for v in vals]
    e = _unit(axis)
    b = u @ u - np.eye(3)
    cof_b = cofactor_matrix(b)
    cc2 = float(a @ u @ cof_b @ n)
    scale = float(np.linalg.norm(a) * np.linalg.norm(u, 2) * np.linalg.norm(cof_b, 2) * np.linalg.norm(n))
    cc2n = abs(cc2) / max(scale, 1e-15)
    cc3 = float(np.trace(u @ u) - np.linalg.det(u @ u) - (a @ a) * (n @ n) / 4.0 - 2.0)
    return CofactorResult(
        domain_type=domain_label,
        axis=e,
        a=np.asarray(a, float),
        n=np.asarray(n, float),
        lambda1=l1,
        lambda2=l2,
        lambda3=l3,
        cc1_residual=l2 - 1.0,
        cc2_raw=cc2,
        cc2_normalized=cc2n,
        cc3_margin=cc3,
        cc2_simplified_residual=float(simplified_residual),
        cc1_pass=abs(l2 - 1.0) <= cc1_tol,
        cc2_pass=cc2n <= cc2_tol,
        cc3_pass=cc3 >= -cc3_tol,
        partner_axis=None if partner_axis is None else _unit(partner_axis),
    )


def compound_domain_rank_one(
    U: Array,
    axis1: Iterable[float],
    axis2: Iterable[float],
    which: Literal[1, 2],
) -> tuple[Array, Array]:
    """Rank-one vectors for a compound domain generated by two perpendicular axes.

    The two axes must generate the same symmetry-related stretch. The formulas
    are the compound-domain solutions from the nonlinear cofactor framework.
    """
    u = np.asarray(U, float)
    e1 = _unit(axis1)
    e2 = _unit(axis2)
    if abs(float(e1 @ e2)) > 1e-7:
        raise ValueError("Compound-domain generating axes must be perpendicular.")
    v1 = twofold_related_stretch(u, e1)
    v2 = twofold_related_stretch(u, e2)
    if not np.allclose(v1, v2, atol=1e-7, rtol=0.0) or np.allclose(v1, u, atol=1e-9, rtol=0.0):
        raise ValueError("Axes do not generate the same non-trivial compound-related stretch.")
    u2 = u @ u
    uinv2 = np.linalg.inv(u2)
    if which == 1:
        denom = float(e1 @ uinv2 @ e1)
        xi = 2.0 * float(e2 @ uinv2 @ e1) / denom
        return xi * (u @ e2), e1
    if which == 2:
        denom = float(e1 @ u2 @ e1)
        eta = -2.0 * float(e2 @ u2 @ e1) / denom
        return eta * (u @ e1), e2
    raise ValueError("which must be 1 or 2")


def compound_cofactor_conditions(
    U: Array,
    axis1: Iterable[float],
    axis2: Iterable[float],
    which: Literal[1, 2],
    cc1_tol: float = 1e-5,
    cc2_tol: float = 1e-5,
    cc3_tol: float = 1e-10,
) -> CofactorResult:
    a, n = compound_domain_rank_one(U, axis1, axis2, which)
    label: DomainLabel = "Compound 1" if which == 1 else "Compound 2"
    axis = axis1 if which == 1 else axis2
    partner = axis2 if which == 1 else axis1
    return _evaluate_cofactor_rank_one(U, axis, label, a, n, cc1_tol, cc2_tol, cc3_tol, partner_axis=partner)


def cofactor_conditions(
    U: Array,
    axis: Iterable[float],
    domain_type: DomainType,
    cc1_tol: float = 1e-5,
    cc2_tol: float = 1e-5,
    cc3_tol: float = 1e-10,
) -> CofactorResult:
    """Evaluate CC1–CC3 for a non-compound Type-I or Type-II domain description.

    When one symmetry-related stretch has two non-parallel two-fold generators,
    it is a compound-domain case and should be evaluated with
    ``compound_cofactor_conditions``. The full CC2 equation remains valid in
    either case, but the simplified Type-I/II residual is only asserted here.
    """
    u = np.asarray(U, float)
    e = _unit(axis)
    a, n = domain_rank_one(u, e, domain_type)
    if domain_type == "Type I":
        simplified = float(np.linalg.norm(np.linalg.solve(u, e)) - 1.0)
    else:
        simplified = float(np.linalg.norm(u @ e) - 1.0)
    return _evaluate_cofactor_rank_one(
        u, e, domain_type, a, n, cc1_tol, cc2_tol, cc3_tol, simplified
    )


def enumerate_domain_specs(U: Array, axes: Iterable[Iterable[float]], tol: float = 1e-8) -> list[CofactorDomainSpec]:
    """Classify non-trivial two-fold generators into Type-I/II or Compound systems.

    If two perpendicular generators produce the same non-trivial U-hat, they
    define one compound-domain pair and are represented by its two rank-one
    solutions rather than four misleading Type-I/Type-II labels.
    """
    u = np.asarray(U, float)
    axis_list = [_unit(a) for a in axes]
    used: set[int] = set()
    specs: list[CofactorDomainSpec] = []
    for i, e in enumerate(axis_list):
        if i in used:
            continue
        v = twofold_related_stretch(u, e)
        if np.allclose(v, u, atol=tol, rtol=0.0):
            used.add(i)
            continue
        same = [
            j for j, f in enumerate(axis_list)
            if j not in used
            and not np.allclose(twofold_related_stretch(u, f), u, atol=tol, rtol=0.0)
            and np.allclose(twofold_related_stretch(u, f), v, atol=tol, rtol=0.0)
        ]
        # A conventional compound pair has exactly two perpendicular generators.
        if len(same) == 2:
            e1, e2 = axis_list[same[0]], axis_list[same[1]]
            if abs(float(e1 @ e2)) <= 1e-7:
                specs.extend([
                    CofactorDomainSpec("Compound 1", e1, e2, 1),
                    CofactorDomainSpec("Compound 2", e2, e1, 2),
                ])
                used.update(same)
                continue
        # Otherwise retain the standard two rank-one descriptions for each
        # distinct generator. This branch also makes unusual symmetry cases
        # visible instead of silently discarding them.
        for j in same or [i]:
            f = axis_list[j]
            specs.append(CofactorDomainSpec("Type I", f))
            specs.append(CofactorDomainSpec("Type II", f))
            used.add(j)
    return specs


def evaluate_domain_spec(
    U: Array,
    spec: CofactorDomainSpec,
    cc1_tol: float = 1e-5,
    cc2_tol: float = 1e-5,
    cc3_tol: float = 1e-10,
) -> CofactorResult:
    if spec.domain_type in ("Type I", "Type II"):
        return cofactor_conditions(U, spec.axis, spec.domain_type, cc1_tol, cc2_tol, cc3_tol)
    if spec.partner_axis is None or spec.compound_solution not in (1, 2):
        raise ValueError("Compound domain specification is incomplete.")
    # Store specs as (normal-axis for that solution, other generator). The
    # compound formula itself expects a consistent e1,e2 ordering.
    if spec.compound_solution == 1:
        e1, e2 = spec.axis, spec.partner_axis
        which = 1
    else:
        # Compound 2 spec.axis is e2 and partner_axis is e1.
        e1, e2 = spec.partner_axis, spec.axis
        which = 2
    return compound_cofactor_conditions(U, e1, e2, which, cc1_tol, cc2_tol, cc3_tol)


def ptmc_summary(inp: LatticeInput) -> dict[str, float]:
    stretch = stretch_from_lattice(inp)
    l1, l2, l3 = stretch.eigenvalues
    return {
        "lambda1": float(l1),
        "lambda2": float(l2),
        "lambda3": float(l3),
        "lambda2_residual": float(l2 - 1.0),
        "det_U": float(stretch.determinant),
        "volume_change": float(stretch.determinant - 1.0),
    }
