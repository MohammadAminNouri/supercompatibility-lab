from __future__ import annotations

"""Round-trip parent/daughter orientation-cycle analysis.

This module is deliberately orientation-focused.  It takes a parent reconstruction and
asks what daughter orientations are crystallographically allowed if the reconstructed
parent transforms again under the *same explicitly supplied orientation relationship*.

For the NiTi B2 <-> B19' workflow a metric-aware helper is provided for the commonly
reported natural/AQ orientation relationship

    (010)_B19' || (110)_B2
    [101]_B19' || [-1 1 1]_B2

The helper converts direct-lattice directions and reciprocal-lattice plane normals to
orthonormal Cartesian crystal frames before constructing the daughter->parent rotation.
This avoids treating monoclinic Miller indices as if they were Euclidean Cartesian
components.

Important scope guard
---------------------
A B19' -> B2 -> B19' round trip is an internal crystallographic consistency test and a
library of allowed re-transformation orientations.  It does *not* by itself predict
which B19' variant will nucleate on a later thermal/mechanical cycle.  Variant selection
also depends on stress, interfaces, defects, thermal path and other microstructural
physics.  Independent second-cycle EBSD can be matched to the regenerated library using
the functions below.
"""

from dataclasses import dataclass
from io import BytesIO
import json
from typing import Iterable, Mapping, Sequence
from zipfile import ZIP_DEFLATED, ZipFile

import numpy as np
import pandas as pd

from .core import C_A_TO_M, C_M_TO_A
from .reconstruction import (
    ReconstructionResult,
    matrix_to_bunge_euler,
    misorientation_deg,
    quaternion_wxyz,
    unique_child_variants,
)

Array = np.ndarray


@dataclass(frozen=True)
class NiTiAQLattice:
    """Lattice parameters used only to construct the metric-aware NiTi AQ OR.

    Lengths may be supplied in any common unit because only direction ratios enter the
    orientation calculation; all four lengths must nevertheless use the same unit.
    """

    a_b2: float = 3.010
    a_b19p: float = 2.898
    b_b19p: float = 4.108
    c_b19p: float = 4.646
    beta_deg: float = 97.78

    def validate(self) -> None:
        vals = np.asarray([self.a_b2, self.a_b19p, self.b_b19p, self.c_b19p], float)
        if not np.all(np.isfinite(vals)) or not np.isfinite(self.beta_deg):
            raise ValueError("All NiTi AQ lattice parameters must be finite numbers.")
        if np.any(vals <= 0):
            raise ValueError("All NiTi AQ lattice lengths must be positive.")
        if not (0.0 < float(self.beta_deg) < 180.0):
            raise ValueError("The B19' monoclinic angle beta must lie strictly between 0 and 180 degrees.")


def _project_so3(m: Array) -> Array:
    u, _, vt = np.linalg.svd(np.asarray(m, float))
    r = u @ vt
    if np.linalg.det(r) < 0:
        u[:, -1] *= -1
        r = u @ vt
    return r


def _crystal_basis_b2(a: float) -> Array:
    return float(a) * np.eye(3)


def _crystal_basis_b19p(a: float, b: float, c: float, beta_deg: float) -> Array:
    """Columns are monoclinic direct-lattice basis vectors a,b,c in a Cartesian frame."""
    beta = np.deg2rad(float(beta_deg))
    return np.array(
        [
            [float(a), 0.0, float(c) * np.cos(beta)],
            [0.0, float(b), 0.0],
            [0.0, 0.0, float(c) * np.sin(beta)],
        ],
        dtype=float,
    )


def _unit(v: Array, label: str) -> Array:
    x = np.asarray(v, float)
    n = float(np.linalg.norm(x))
    if n <= 1e-14:
        raise ValueError(f"{label} cannot be the zero vector.")
    return x / n


def _cartesian_direction(basis: Array, uvw: Sequence[float]) -> Array:
    return _unit(np.asarray(basis, float) @ np.asarray(uvw, float), "crystallographic direction")


def _cartesian_plane_normal(basis: Array, hkl: Sequence[float]) -> Array:
    # Reciprocal-basis vector; the common 2*pi factor is irrelevant for a direction.
    return _unit(np.linalg.inv(np.asarray(basis, float)).T @ np.asarray(hkl, float), "plane normal")


def _orthonormal_triad(normal: Array, direction: Array) -> Array:
    n = _unit(normal, "plane normal")
    d = np.asarray(direction, float)
    # Remove only numerical/non-ideal leakage. For a physically valid plane/direction
    # parallelism this component is zero within floating-point precision.
    d = d - n * float(n @ d)
    d = _unit(d, "in-plane direction")
    t = _unit(np.cross(n, d), "second in-plane direction")
    return np.column_stack([d, t, n])


def niti_aq_orientation_relationship(lattice: NiTiAQLattice = NiTiAQLattice()) -> Array:
    """Return the metric-aware B19'->B2 rotation for the natural/AQ NiTi OR.

    The implemented parallelisms are

        (010)_B19' || (110)_B2
        [101]_B19' || [-1 1 1]_B2

    and the result is a proper rotation R_cp mapping B19' crystal-Cartesian coordinates
    into B2 crystal-Cartesian coordinates.  The lattice metric matters because [101] in
    a monoclinic basis is not the Euclidean vector (1,0,1).
    """
    lattice.validate()
    bp = _crystal_basis_b2(lattice.a_b2)
    bc = _crystal_basis_b19p(lattice.a_b19p, lattice.b_b19p, lattice.c_b19p, lattice.beta_deg)

    n_parent = _cartesian_plane_normal(bp, [1.0, 1.0, 0.0])
    d_parent = _cartesian_direction(bp, [-1.0, 1.0, 1.0])
    n_child = _cartesian_plane_normal(bc, [0.0, 1.0, 0.0])
    d_child = _cartesian_direction(bc, [1.0, 0.0, 1.0])

    parent_triad = _orthonormal_triad(n_parent, d_parent)
    child_triad = _orthonormal_triad(n_child, d_child)
    return _project_so3(parent_triad @ child_triad.T)



def niti_ct_otsuka_ren_orientation_relationship(lattice: NiTiAQLattice = NiTiAQLattice()) -> Array:
    """Return a model-derived B19'->B2 initial OR from CT/Otsuka--Ren correspondence.

    The Otsuka--Ren correspondence used by the built-in NiTi model is the exact
    crystallographic correspondence matrix stored in :mod:`src.core`.  A
    correspondence is *not*, by itself, a unique experimental orientation
    relationship.  To bridge the metric/correspondence model to the orientation
    reconstruction engine in a reproducible way, we form the two-point lattice
    deformation

        F_{M<-A} = B_M C^{M->A} B_A^{-1}

    from orthonormal Cartesian basis matrices B_A (B2) and B_M (B19'), take its
    right polar decomposition F = R_{A->M} U, and return

        R_cp = R_{A->M}^T,

    the daughter(B19')->parent(B2) proper rotation required by the reconstruction
    engine.  This is therefore a *model-derived initial OR*, not a claim that CT
    or the Otsuka--Ren correspondence uniquely determines the experimentally
    measured OR.  When an experimental OR is available it should be supplied or
    used to refine this starting value.
    """
    lattice.validate()
    b_parent = _crystal_basis_b2(lattice.a_b2)
    b_child = _crystal_basis_b19p(lattice.a_b19p, lattice.b_b19p, lattice.c_b19p, lattice.beta_deg)
    f_child_from_parent = b_child @ np.asarray(C_M_TO_A, float) @ np.linalg.inv(b_parent)
    u, _, vt = np.linalg.svd(f_child_from_parent)
    r_parent_to_child = u @ vt
    if np.linalg.det(r_parent_to_child) < 0:
        u[:, -1] *= -1
        r_parent_to_child = u @ vt
    return _project_so3(r_parent_to_child.T)


def niti_ct_otsuka_ren_diagnostics(
    lattice: NiTiAQLattice = NiTiAQLattice(),
    r_child_to_parent: Array | None = None,
) -> dict[str, object]:
    """Return auditable matrices/diagnostics for the CT/Otsuka--Ren bridge.

    The output deliberately exposes the exact correspondence, deformation,
    principal stretches and distance to the separately implemented natural/AQ OR.
    No single scalar is promoted to a pass/fail verdict.
    """
    lattice.validate()
    b_parent = _crystal_basis_b2(lattice.a_b2)
    b_child = _crystal_basis_b19p(lattice.a_b19p, lattice.b_b19p, lattice.c_b19p, lattice.beta_deg)
    f = b_child @ np.asarray(C_M_TO_A, float) @ np.linalg.inv(b_parent)
    rcp = niti_ct_otsuka_ren_orientation_relationship(lattice) if r_child_to_parent is None else _project_so3(r_child_to_parent)
    stretches = np.linalg.svd(f, compute_uv=False)
    stretches = np.sort(np.asarray(stretches, float))
    aq = niti_aq_orientation_relationship(lattice)
    rel = _project_so3(rcp @ aq.T)
    c = float(np.clip((np.trace(rel) - 1.0) / 2.0, -1.0, 1.0))
    return {
        "correspondence_A_to_M": np.asarray(C_A_TO_M, float).copy(),
        "correspondence_M_to_A": np.asarray(C_M_TO_A, float).copy(),
        "F_child_from_parent": np.asarray(f, float).copy(),
        "principal_stretches_sorted": stretches,
        "det_F": float(np.linalg.det(f)),
        "det_R_cp": float(np.linalg.det(rcp)),
        "orthogonality_frobenius_residual": float(np.linalg.norm(rcp.T @ rcp - np.eye(3), ord="fro")),
        "misorientation_to_natural_AQ_OR_deg": float(np.degrees(np.arccos(c))),
    }

def niti_aq_parallelism_residuals(lattice: NiTiAQLattice, r_child_to_parent: Array | None = None) -> dict[str, float]:
    """Angular residuals (degrees) of the two defining AQ parallelisms."""
    lattice.validate()
    r = niti_aq_orientation_relationship(lattice) if r_child_to_parent is None else _project_so3(r_child_to_parent)
    bp = _crystal_basis_b2(lattice.a_b2)
    bc = _crystal_basis_b19p(lattice.a_b19p, lattice.b_b19p, lattice.c_b19p, lattice.beta_deg)
    np_ = _cartesian_plane_normal(bp, [1, 1, 0])
    dp = _cartesian_direction(bp, [-1, 1, 1])
    nc = _cartesian_plane_normal(bc, [0, 1, 0])
    dc = _cartesian_direction(bc, [1, 0, 1])

    def unsigned_angle(a: Array, b: Array) -> float:
        c = float(np.clip(abs(float(_unit(a, "a") @ _unit(b, "b"))), -1.0, 1.0))
        return float(np.degrees(np.arccos(c)))

    return {
        "plane_parallelism_residual_deg": unsigned_angle(r @ nc, np_),
        "direction_parallelism_residual_deg": unsigned_angle(r @ dc, dp),
        "det_R_cp": float(np.linalg.det(r)),
        "orthogonality_frobenius_residual": float(np.linalg.norm(r.T @ r - np.eye(3), ord="fro")),
    }


def _variant_rows(parent_id: int, variants: list[Array]) -> list[dict]:
    rows: list[dict] = []
    for vid, g in enumerate(variants, 1):
        p1, P, p2 = matrix_to_bunge_euler(g)
        qw, qx, qy, qz = quaternion_wxyz(g)
        row = {
            "parent_id": int(parent_id),
            "regenerated_daughter_variant_id": int(vid),
            "daughter_Bunge_phi1_deg": p1,
            "daughter_Bunge_Phi_deg": P,
            "daughter_Bunge_phi2_deg": p2,
            "daughter_quaternion_qw": qw,
            "daughter_quaternion_qx": qx,
            "daughter_quaternion_qy": qy,
            "daughter_quaternion_qz": qz,
        }
        for i in range(3):
            for j in range(3):
                row[f"g{i+1}{j+1}"] = float(g[i, j])
        rows.append(row)
    return rows


def regenerated_variant_library(
    parent_orientations: Mapping[int, Array],
    r_child_to_parent: Array,
    parent_sym: Iterable[Array],
    child_sym: Iterable[Array],
) -> tuple[pd.DataFrame, dict[int, list[Array]]]:
    """Generate every symmetry-distinct daughter orientation from every parent."""
    cache: dict[int, list[Array]] = {}
    rows: list[dict] = []
    for pid in sorted(int(x) for x in parent_orientations):
        variants = unique_child_variants(
            np.asarray(parent_orientations[pid], float),
            np.asarray(r_child_to_parent, float),
            parent_sym,
            child_sym,
        )
        cache[pid] = variants
        rows.extend(_variant_rows(pid, variants))
    return pd.DataFrame(rows), cache


def _fit_quality(x: float, strong_deg: float, acceptable_deg: float, review_deg: float) -> str:
    if x <= strong_deg:
        return "strong closure"
    if x <= acceptable_deg:
        return "acceptable closure"
    if x <= review_deg:
        return "review"
    return "poor closure"


def cycle_closure_table(
    grain_ids: Sequence[int],
    measured_child_orientations: Sequence[Array],
    reconstruction: ReconstructionResult,
    r_child_to_parent: Array,
    parent_sym: Iterable[Array],
    child_sym: Iterable[Array],
    *,
    strong_deg: float = 1.0,
    acceptable_deg: float = 2.5,
    review_deg: float = 5.0,
) -> pd.DataFrame:
    """Close measured daughter -> reconstructed parent -> regenerated daughter.

    For daughter grain i assigned to reconstructed parent P_i, the closure residual is

        delta_i = min_k d_D(g_D,i^meas, g_D,P_i^(k))

    where d_D is the symmetry-reduced daughter disorientation and k enumerates every
    symmetry-distinct regenerated daughter orientation of that parent.
    """
    if not (0 <= strong_deg <= acceptable_deg <= review_deg):
        raise ValueError("Closure thresholds must satisfy 0 <= strong <= acceptable <= review.")
    gids = np.asarray(grain_ids, int)
    measured = list(measured_child_orientations)
    if len(gids) != len(measured) or len(gids) != len(reconstruction.parent_ids):
        raise ValueError("grain_ids, measured daughter orientations and reconstruction assignments must have equal length.")

    _, cache = regenerated_variant_library(reconstruction.parent_orientations, r_child_to_parent, parent_sym, child_sym)
    rows: list[dict] = []
    for i, (gid, g_meas) in enumerate(zip(gids, measured)):
        pid = int(reconstruction.parent_ids[i])
        variants = cache[pid]
        distances = np.asarray([misorientation_deg(np.asarray(g_meas, float), g, child_sym) for g in variants], float)
        order = np.argsort(distances)
        best = int(order[0])
        second = float(distances[order[1]]) if len(order) > 1 else np.nan
        best_fit = float(distances[best])
        row = {
            "daughter_grain_id": int(gid),
            "reconstructed_parent_id": pid,
            "observed_regenerated_branch_id": best + 1,
            "cycle_closure_misorientation_deg": best_fit,
            "second_best_regenerated_branch_fit_deg": second,
            "regenerated_branch_separation_deg": float(second - best_fit) if np.isfinite(second) else np.nan,
            "cycle_closure_quality": _fit_quality(best_fit, strong_deg, acceptable_deg, review_deg),
            "original_parent_reconstruction_OR_fit_deg": float(reconstruction.fit_deg[i]),
            "original_parent_reconstruction_support_score": float(reconstruction.confidence[i]),
            "allowed_regenerated_branches_for_parent": int(len(variants)),
        }
        rows.append(row)
    return pd.DataFrame(rows)


def parent_cycle_summary(closure: pd.DataFrame, variant_library: pd.DataFrame) -> pd.DataFrame:
    """Parent-level summary of round-trip closure and currently observed branch coverage."""
    rows: list[dict] = []
    for pid, grp in closure.groupby("reconstructed_parent_id", sort=True):
        vals = pd.to_numeric(grp["cycle_closure_misorientation_deg"], errors="coerce").to_numpy(float)
        allowed = int((variant_library["parent_id"] == int(pid)).sum())
        observed = int(grp["observed_regenerated_branch_id"].nunique())
        rows.append({
            "parent_id": int(pid),
            "measured_daughter_grains": int(len(grp)),
            "allowed_regenerated_daughter_branches": allowed,
            "currently_observed_branches": observed,
            "observed_branch_coverage_pct": float(100.0 * observed / allowed) if allowed else np.nan,
            "mean_cycle_closure_deg": float(np.nanmean(vals)),
            "median_cycle_closure_deg": float(np.nanmedian(vals)),
            "P95_cycle_closure_deg": float(np.nanpercentile(vals, 95)),
            "max_cycle_closure_deg": float(np.nanmax(vals)),
            "fraction_closing_within_1deg": float(np.mean(vals <= 1.0)),
            "fraction_closing_within_2_5deg": float(np.mean(vals <= 2.5)),
            "fraction_closing_within_5deg": float(np.mean(vals <= 5.0)),
        })
    return pd.DataFrame(rows)


def observed_branch_occupancy(closure: pd.DataFrame, source_df: pd.DataFrame | None = None) -> pd.DataFrame:
    """Observed regenerated-branch counts, optionally area weighted."""
    df = closure.copy()
    area_col = None
    if source_df is not None and "grain_id" in source_df.columns:
        for c in ("area_um2", "area_est_um2", "area_est_mapunit2", "area"):
            if c in source_df.columns:
                area_col = c
                break
        if area_col:
            df = df.merge(
                source_df[["grain_id", area_col]].rename(columns={"grain_id": "daughter_grain_id"}),
                on="daughter_grain_id",
                how="left",
            )
    rows: list[dict] = []
    for (pid, vid), grp in df.groupby(["reconstructed_parent_id", "observed_regenerated_branch_id"], sort=True):
        parent_total = int((df["reconstructed_parent_id"] == pid).sum())
        row = {
            "parent_id": int(pid),
            "observed_regenerated_branch_id": int(vid),
            "daughter_grains": int(len(grp)),
            "grain_fraction_within_parent_pct": float(100.0 * len(grp) / parent_total) if parent_total else np.nan,
            "mean_cycle_closure_deg": float(pd.to_numeric(grp["cycle_closure_misorientation_deg"], errors="coerce").mean()),
        }
        if area_col:
            area = pd.to_numeric(grp[area_col], errors="coerce").fillna(0.0)
            parent_area = pd.to_numeric(df.loc[df["reconstructed_parent_id"] == pid, area_col], errors="coerce").fillna(0.0).sum()
            row["daughter_area"] = float(area.sum())
            row["area_fraction_within_parent_pct"] = float(100.0 * area.sum() / parent_area) if parent_area > 0 else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def branch_switch_catalog(
    parent_orientations: Mapping[int, Array],
    r_child_to_parent: Array,
    parent_sym: Iterable[Array],
    child_sym: Iterable[Array],
) -> pd.DataFrame:
    """Pairwise orientation changes between all regenerated daughter branches.

    This is a geometric switch catalogue, not a nucleation-probability model.  Every row
    states how much daughter orientation would change if a parent transformed from one
    regenerated branch to another.
    """
    _, cache = regenerated_variant_library(parent_orientations, r_child_to_parent, parent_sym, child_sym)
    rows: list[dict] = []
    for pid, variants in cache.items():
        for i, ga in enumerate(variants, 1):
            for j, gb in enumerate(variants, 1):
                ang = 0.0 if i == j else float(misorientation_deg(ga, gb, child_sym))
                rows.append({
                    "parent_id": int(pid),
                    "from_regenerated_branch_id": int(i),
                    "to_regenerated_branch_id": int(j),
                    "daughter_orientation_change_deg": ang,
                    "transition_class": "same regenerated branch" if i == j else "alternative regenerated branch",
                })
    return pd.DataFrame(rows)


def match_new_cycle_daughters(
    grain_ids: Sequence[int],
    new_child_orientations: Sequence[Array],
    parent_orientations: Mapping[int, Array],
    r_child_to_parent: Array,
    parent_sym: Iterable[Array],
    child_sym: Iterable[Array],
    known_parent_ids: Sequence[int | float | None] | None = None,
) -> pd.DataFrame:
    """Match independently measured new-cycle daughter grains to the regenerated library.

    If ``known_parent_ids`` is supplied, each grain is tested only against that parent.
    Otherwise every reconstructed parent is considered and the best/second-best parent
    separation is reported so orientation-only parent assignment ambiguity remains visible.
    """
    gids = np.asarray(grain_ids, int)
    measured = list(new_child_orientations)
    if len(gids) != len(measured):
        raise ValueError("New-cycle grain IDs and orientations must have equal length.")
    if known_parent_ids is not None and len(known_parent_ids) != len(gids):
        raise ValueError("known_parent_ids must have one entry per new-cycle daughter grain.")

    _, cache = regenerated_variant_library(parent_orientations, r_child_to_parent, parent_sym, child_sym)
    parent_ids = sorted(cache)
    rows: list[dict] = []
    for idx, (gid, g) in enumerate(zip(gids, measured)):
        known = None
        if known_parent_ids is not None:
            val = known_parent_ids[idx]
            if val is not None and not (isinstance(val, float) and np.isnan(val)):
                known = int(val)
                if known not in cache:
                    raise ValueError(f"New-cycle grain {gid} references unknown parent_id={known}.")
        candidates = [known] if known is not None else parent_ids
        parent_best: list[tuple[float, int, int, float]] = []
        for pid in candidates:
            distances = np.asarray([misorientation_deg(np.asarray(g, float), gv, child_sym) for gv in cache[int(pid)]], float)
            order = np.argsort(distances)
            b = int(order[0])
            second_variant = float(distances[order[1]]) if len(order) > 1 else np.nan
            parent_best.append((float(distances[b]), int(pid), b + 1, second_variant))
        parent_best.sort(key=lambda x: (x[0], x[1], x[2]))
        best_fit, best_pid, best_vid, second_variant_fit = parent_best[0]
        second_parent_fit = float(parent_best[1][0]) if len(parent_best) > 1 else np.nan
        rows.append({
            "new_cycle_daughter_grain_id": int(gid),
            "matched_reconstructed_parent_id": int(best_pid),
            "matched_regenerated_branch_id": int(best_vid),
            "new_cycle_OR_library_fit_deg": float(best_fit),
            "second_best_branch_fit_within_parent_deg": float(second_variant_fit),
            "best_vs_second_branch_separation_deg": float(second_variant_fit - best_fit) if np.isfinite(second_variant_fit) else np.nan,
            "second_best_parent_fit_deg": second_parent_fit,
            "best_vs_second_parent_separation_deg": float(second_parent_fit - best_fit) if np.isfinite(second_parent_fit) else np.nan,
            "parent_assignment_basis": "user-supplied parent_id" if known is not None else "orientation-only best match across reconstructed parents",
        })
    return pd.DataFrame(rows)


def cycle_evidence_zip(
    *,
    variant_library: pd.DataFrame,
    closure: pd.DataFrame,
    parent_summary: pd.DataFrame,
    occupancy: pd.DataFrame,
    switch_catalog: pd.DataFrame,
    metadata: Mapping[str, object],
    new_cycle_matches: pd.DataFrame | None = None,
) -> bytes:
    """Create a self-contained academic CSV/JSON evidence bundle."""
    bio = BytesIO()
    with ZipFile(bio, "w", ZIP_DEFLATED) as z:
        z.writestr("regenerated_daughter_variant_library.csv", variant_library.to_csv(index=False))
        z.writestr("measured_round_trip_cycle_closure.csv", closure.to_csv(index=False))
        z.writestr("parent_cycle_summary.csv", parent_summary.to_csv(index=False))
        z.writestr("observed_branch_occupancy.csv", occupancy.to_csv(index=False))
        z.writestr("regenerated_branch_switch_catalog.csv", switch_catalog.to_csv(index=False))
        if new_cycle_matches is not None:
            z.writestr("independent_new_cycle_daughter_matches.csv", new_cycle_matches.to_csv(index=False))
        z.writestr("cycle_metadata.json", json.dumps(dict(metadata), indent=2, sort_keys=True))
        z.writestr(
            "README.txt",
            "B19' / parent / regenerated-B19' orientation-cycle evidence bundle\n\n"
            "The round-trip closure residual is an internal crystallographic consistency metric because the same supplied OR is used for reconstruction and regeneration. "
            "It must not be described as independent validation. Independent evidence can come from retained-parent measurements or a separately measured later-cycle daughter EBSD map.\n\n"
            "The branch-switch catalogue contains crystallographically allowed orientation changes only. It does not assign nucleation probabilities or predict future volume fractions/morphology.\n",
        )
    return bio.getvalue()
