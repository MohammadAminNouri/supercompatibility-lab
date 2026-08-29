from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .core import C_A_TO_M, Array
from .symmetry import commuting_variant_pairs, cubic_point_group, stretch_variants


EXTREME_MONOCLINIC_II_TARGET_1 = np.array([
    [1.0, 0.0, 0.0],
    [0.0, 3.0/np.sqrt(2.0), -1.0/np.sqrt(2.0)],
    [0.0, -1.0/np.sqrt(2.0), 1.0/np.sqrt(2.0)],
])
EXTREME_MONOCLINIC_II_TARGET_2 = np.array([
    [1.0, 0.0, 0.0],
    [0.0, 1.0/np.sqrt(2.0), 1.0/np.sqrt(2.0)],
    [0.0, 1.0/np.sqrt(2.0), 3.0/np.sqrt(2.0)],
])


@dataclass(frozen=True)
class FrontierDiagnostic:
    n_stretch_variants: int
    commuting_pairs: list[tuple[int, int, float]]
    monoclinic_axis_parent_direction: Array
    built_in_is_monoclinic_ii: bool
    extreme_target_distance_if_applicable: float | None
    note: str


def _canonical_direction(v: Array) -> Array:
    x = np.asarray(v, float)
    nz = np.abs(x[np.abs(x) > 1e-10])
    if len(nz) == 0:
        return x
    y = x / np.min(nz)
    y = np.rint(y).astype(int).astype(float)
    for z in y:
        if z != 0:
            if z < 0:
                y *= -1
            break
    return y


def monoclinic_unique_axis_in_parent() -> Array:
    """Parent direction that corresponds to the B19′ unique b axis [010]."""
    e_b_m = np.array([0.0, 1.0, 0.0])
    # u_M = C_M→A u_A, hence u_A = C_A→M u_M.
    return _canonical_direction(C_A_TO_M @ e_b_m)


def _extreme_distance(U: Array) -> float:
    targets = (EXTREME_MONOCLINIC_II_TARGET_1, EXTREME_MONOCLINIC_II_TARGET_2)
    best = float("inf")
    for g in cubic_point_group():
        ug = g @ U @ g.T
        for t in targets:
            scale = max(float(np.linalg.norm(t, ord="fro")), 1e-12)
            best = min(best, float(np.linalg.norm(ug - t, ord="fro") / scale))
    return best


def frontier_diagnostic(U: Array, commutation_tol: float = 1e-9) -> FrontierDiagnostic:
    axis = monoclinic_unique_axis_in_parent()
    # The 2026 extreme-compatibility theory defines cubic→monoclinic-II by the
    # monoclinic symmetry axis corresponding to a cubic <100> direction.
    nz = np.flatnonzero(np.abs(axis) > 1e-8)
    is_mii = len(nz) == 1
    pairs = commuting_variant_pairs(U, tol=commutation_tol)
    distance = _extreme_distance(U) if is_mii else None
    if is_mii:
        note = (
            "Built-in correspondence is in the cubic→monoclinic-II applicability class. "
            "Distance is a stretch-space diagnostic to the two published exact extreme-compatible targets; "
            "it is not a substitute for checking the full theorem hypotheses."
        )
    else:
        note = (
            "The built-in B2→B19′ correspondence maps the monoclinic unique b axis to a cubic <110> direction, "
            "not <100>. Therefore the published cubic→monoclinic-II extreme-compatibility target tensors are "
            "not applied as a certification criterion here. Commuting-variant diagnostics are still reported."
        )
    return FrontierDiagnostic(
        n_stretch_variants=len(stretch_variants(U)),
        commuting_pairs=pairs,
        monoclinic_axis_parent_direction=axis,
        built_in_is_monoclinic_ii=is_mii,
        extreme_target_distance_if_applicable=distance,
        note=note,
    )
