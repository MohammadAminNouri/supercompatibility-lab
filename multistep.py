from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .core import Array, cmc_degeneracy
from .ptmc import symmetric_sqrt


@dataclass(frozen=True)
class GeneralCell:
    a: float
    b: float
    c: float
    alpha_deg: float = 90.0
    beta_deg: float = 90.0
    gamma_deg: float = 90.0


@dataclass(frozen=True)
class StageResult:
    parent_metric: Array
    product_metric: Array
    correspondence: Array
    cmc: Array
    cmc_relative_zero: float
    cmc_order: int
    lambda1: float
    lambda2: float
    lambda3: float
    det_u: float


def metric_from_cell(cell: GeneralCell) -> Array:
    if min(cell.a, cell.b, cell.c) <= 0:
        raise ValueError("Cell lengths must be positive.")
    for angle in (cell.alpha_deg, cell.beta_deg, cell.gamma_deg):
        if not 0 < angle < 180:
            raise ValueError("Cell angles must lie between 0 and 180 degrees.")
    alpha, beta, gamma = np.deg2rad([cell.alpha_deg, cell.beta_deg, cell.gamma_deg])
    return np.array([
        [cell.a**2, cell.a*cell.b*np.cos(gamma), cell.a*cell.c*np.cos(beta)],
        [cell.a*cell.b*np.cos(gamma), cell.b**2, cell.b*cell.c*np.cos(alpha)],
        [cell.a*cell.c*np.cos(beta), cell.b*cell.c*np.cos(alpha), cell.c**2],
    ], dtype=float)


def evaluate_stage(parent: GeneralCell, product: GeneralCell, correspondence: Array) -> StageResult:
    ma = metric_from_cell(parent)
    mm = metric_from_cell(product)
    c = np.asarray(correspondence, float)
    if c.shape != (3, 3) or abs(np.linalg.det(c)) < 1e-12:
        raise ValueError("Correspondence must be an invertible 3×3 matrix mapping parent directions to product coordinates.")
    cmc = c.T @ mm @ c - ma
    cmc = 0.5 * (cmc + cmc.T)
    deg = cmc_degeneracy(cmc, rtol=1e-6)

    a_half = symmetric_sqrt(ma)
    a_inv_half = np.linalg.inv(a_half)
    target = c.T @ mm @ c
    gram = a_inv_half @ target @ a_inv_half
    U = symmetric_sqrt(gram)
    vals = np.linalg.eigvalsh(U)
    return StageResult(
        parent_metric=ma,
        product_metric=mm,
        correspondence=c,
        cmc=cmc,
        cmc_relative_zero=deg.relative_zero,
        cmc_order=deg.order,
        lambda1=float(vals[0]),
        lambda2=float(vals[1]),
        lambda3=float(vals[2]),
        det_u=float(np.linalg.det(U)),
    )


def evaluate_chain(stages: list[tuple[str, GeneralCell, GeneralCell, Array]]) -> pd.DataFrame:
    rows = []
    for name, parent, product, correspondence in stages:
        r = evaluate_stage(parent, product, correspondence)
        rows.append({
            "stage": name,
            "lambda1": r.lambda1,
            "lambda2": r.lambda2,
            "lambda3": r.lambda3,
            "abs_lambda2_minus_1": abs(r.lambda2 - 1.0),
            "det_U": r.det_u,
            "cmc_relative_zero": r.cmc_relative_zero,
            "cmc_order": r.cmc_order,
        })
    return pd.DataFrame(rows)
