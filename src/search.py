from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import differential_evolution

from .core import (
    LatticeInput,
    cmc_degeneracy,
    cmc_matrix,
    normalized_metrics,
    rank_matches,
)


@dataclass(frozen=True)
class SearchResult:
    a_ratio: float
    b_ratio: float
    c_ratio: float
    beta_deg: float
    a_b19p: float
    b_b19p: float
    c_b19p: float
    epsilon: float
    habit_plane: np.ndarray
    twin_label: str
    parent_element: str
    inequality_margin: float


def search_c1(
    a_b2: float,
    twin_label: str,
    a_bounds: tuple[float, float] = (0.82, 1.05),
    c_bounds: tuple[float, float] = (1.40, 1.80),
    beta_bounds: tuple[float, float] = (94.0, 103.0),
    seed: int = 42,
    maxiter: int = 55,
) -> SearchResult:
    """Search a, c, beta while enforcing C1 with b = sqrt(2)."""

    sqrt2 = float(np.sqrt(2.0))

    def objective(x: np.ndarray) -> float:
        a, c, beta = [float(v) for v in x]
        inp = LatticeInput(a_b2, a * a_b2, sqrt2 * a_b2, c * a_b2, beta)
        ma, mm = normalized_metrics(inp)
        deg = cmc_degeneracy(cmc_matrix(ma, mm), rtol=1e-5)
        if not deg.habit_planes:
            return 50.0 + deg.relative_zero * 1000.0
        beta_r = np.deg2rad(beta)
        margin = 2 * a * a + c * c - a * a * c * c * np.sin(beta_r) ** 2 - 2
        penalty = 100.0 * max(0.0, -margin)
        matches = [m for m in rank_matches(inp) if m.twin.label == twin_label]
        if not matches:
            return 50.0 + penalty
        return float(matches[0].epsilon + penalty)

    opt = differential_evolution(
        objective,
        bounds=[a_bounds, c_bounds, beta_bounds],
        seed=seed,
        polish=True,
        updating="immediate",
        workers=1,
        maxiter=maxiter,
        popsize=10,
        tol=1e-8,
    )
    a, c, beta = [float(v) for v in opt.x]
    inp = LatticeInput(a_b2, a * a_b2, sqrt2 * a_b2, c * a_b2, beta)
    matches = [m for m in rank_matches(inp) if m.twin.label == twin_label]
    if not matches:
        raise RuntimeError("Search ended without a valid habit-plane/twin match.")
    best = matches[0]
    beta_r = np.deg2rad(beta)
    margin = float(2 * a * a + c * c - a * a * c * c * np.sin(beta_r) ** 2 - 2)
    return SearchResult(
        a_ratio=a,
        b_ratio=sqrt2,
        c_ratio=c,
        beta_deg=beta,
        a_b19p=a * a_b2,
        b_b19p=sqrt2 * a_b2,
        c_b19p=c * a_b2,
        epsilon=best.epsilon,
        habit_plane=best.habit_plane,
        twin_label=best.twin.label,
        parent_element=best.twin.parent_element,
        inequality_margin=margin,
    )
