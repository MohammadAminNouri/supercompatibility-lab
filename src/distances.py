from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .core import (
    LatticeInput,
    cmc_degeneracy,
    cmc_matrix,
    compatibility_conditions,
    normalized_metrics,
    shear_match,
    smc_matrix,
)
from .ptmc import CofactorResult, enumerate_domain_specs, evaluate_domain_spec, stretch_from_metrics
from .symmetry import full_twin_explorer, parent_twofold_axes


@dataclass(frozen=True)
class CompatibilityDashboard:
    lambda1: float
    lambda2: float
    lambda3: float
    lambda2_abs_residual: float
    det_u: float
    cmc_relative_zero: float
    cmc_order: int
    cmc_exact: bool
    best_cc2_normalized: float
    best_cc3_margin: float
    best_cofactor_type: str
    best_cofactor_axis: np.ndarray
    ptmc_all_pass: bool
    best_epsilon: float | None
    best_ct_twin: str | None
    c1_equality_distance: float
    c1_inequality_margin: float


def best_cofactor_system(inp: LatticeInput, cc1_tol: float = 1e-5, cc2_tol: float = 1e-5) -> CofactorResult:
    candidates = all_cofactor_systems(inp, cc1_tol=cc1_tol, cc2_tol=cc2_tol)
    if not candidates:
        raise ValueError("No non-trivial cofactor domain systems were generated.")
    # CC1 is common to all systems for a fixed U. Prefer CC3-valid systems, then
    # the smallest full normalized CC2 residual. NaN simplified residuals belong
    # to compound domains and are not used for ranking.
    candidates.sort(key=lambda r: (0 if r.cc3_pass else 1, r.cc2_normalized))
    return candidates[0]


def all_cofactor_systems(inp: LatticeInput, cc1_tol: float = 1e-5, cc2_tol: float = 1e-5) -> list[CofactorResult]:
    ma, mm = normalized_metrics(inp)
    stretch = stretch_from_metrics(ma, mm)
    specs = enumerate_domain_specs(stretch.U, parent_twofold_axes())
    out = [evaluate_domain_spec(stretch.U, spec, cc1_tol=cc1_tol, cc2_tol=cc2_tol) for spec in specs]
    return sorted(out, key=lambda r: (r.cc2_normalized, -r.cc3_margin, r.domain_type))


def compatibility_dashboard(inp: LatticeInput, cc1_tol: float = 1e-5, cc2_tol: float = 1e-5) -> CompatibilityDashboard:
    ma, mm = normalized_metrics(inp)
    cmc = cmc_matrix(ma, mm)
    deg = cmc_degeneracy(cmc)
    stretch = stretch_from_metrics(ma, mm)
    best_cc = best_cofactor_system(inp, cc1_tol=cc1_tol, cc2_tol=cc2_tol)

    best_epsilon: float | None = None
    best_twin: str | None = None
    if deg.habit_planes:
        smc = smc_matrix(ma, mm)
        matches = []
        for entry in full_twin_explorer(ma, mm):
            for plane in deg.habit_planes:
                m = shear_match(plane, smc, entry.twin, ma)
                if np.isfinite(m.epsilon):
                    matches.append((m.epsilon, f"{entry.double_coset} · {entry.parent_element}"))
        if matches:
            matches.sort(key=lambda x: x[0])
            best_epsilon, best_twin = matches[0]

    conditions = compatibility_conditions(inp)
    c1 = next(x for x in conditions if x["name"] == "C1")
    l1, l2, l3 = [float(v) for v in stretch.eigenvalues]
    return CompatibilityDashboard(
        lambda1=l1,
        lambda2=l2,
        lambda3=l3,
        lambda2_abs_residual=abs(l2 - 1.0),
        det_u=float(stretch.determinant),
        cmc_relative_zero=float(deg.relative_zero),
        cmc_order=int(deg.order),
        cmc_exact=bool(deg.exact),
        best_cc2_normalized=float(best_cc.cc2_normalized),
        best_cc3_margin=float(best_cc.cc3_margin),
        best_cofactor_type=best_cc.domain_type,
        best_cofactor_axis=best_cc.axis,
        ptmc_all_pass=bool(best_cc.all_pass),
        best_epsilon=best_epsilon,
        best_ct_twin=best_twin,
        c1_equality_distance=abs(float(c1["equality_residual"])),
        c1_inequality_margin=float(c1["inequality_margin"]),
    )


def dashboard_rows(inp: LatticeInput, cc1_tol: float = 1e-5, cc2_tol: float = 1e-5) -> list[dict[str, object]]:
    d = compatibility_dashboard(inp, cc1_tol=cc1_tol, cc2_tol=cc2_tol)
    return [
        {
            "criterion": "PTMC middle stretch",
            "quantity": "|λ₂ − 1|",
            "value": d.lambda2_abs_residual,
            "target": "0",
            "interpretation": "Single-variant A/M invariant-plane compatibility condition.",
        },
        {
            "criterion": "CMC degeneracy",
            "quantity": "min|qᵢ| / max|qᵢ|",
            "value": d.cmc_relative_zero,
            "target": "0 plus opposite-sign remaining eigenvalues",
            "interpretation": "Distance of the CMC quadratic form from an exact plane degeneracy.",
        },
        {
            "criterion": "Cofactor CC2",
            "quantity": "normalized residual",
            "value": d.best_cc2_normalized,
            "target": "0",
            "interpretation": "Best full cofactor residual across the classified Type-I, Type-II and Compound systems.",
        },
        {
            "criterion": "Cofactor CC3",
            "quantity": "margin",
            "value": d.best_cc3_margin,
            "target": "≥ 0",
            "interpretation": "Sufficiency inequality for the same best cofactor system.",
        },
        {
            "criterion": "Metric shear/shear",
            "quantity": "ε",
            "value": np.nan if d.best_epsilon is None else d.best_epsilon,
            "target": "0",
            "interpretation": "Defined once an exact CMC habit plane exists; compares A/M IPS shear with M/M twin shear.",
        },
    ]
