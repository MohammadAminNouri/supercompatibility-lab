from __future__ import annotations

"""Academic diagnostics and exports for parent/daughter reconstruction.

The numerical reconstruction algorithms live in :mod:`src.reconstruction`.  This module
adds diagnostics that are useful when judging whether a reconstruction is defensible:
method-to-method cluster agreement, prior-parent boundary agreement, matched parent
orientation agreement, variant/operator statistics, and a reproducible export bundle.

None of the agreement metrics is promoted to a universal 'best method' score.  They answer
different questions and are intended to be inspected together.
"""

from dataclasses import dataclass
from io import BytesIO
import json
from typing import Iterable
from zipfile import ZIP_DEFLATED, ZipFile

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

from .reconstruction import (
    ReconstructionResult,
    matrix_to_bunge_euler,
    misorientation_deg,
    quaternion_wxyz,
    theoretical_child_operators,
    two_sided_operator_distance_deg,
)


@dataclass(frozen=True)
class AgreementSummary:
    ari: pd.DataFrame
    nmi: pd.DataFrame
    boundary_jaccard: pd.DataFrame
    matched_parent_orientation_deg: pd.DataFrame


def _names(results: dict[str, ReconstructionResult]) -> list[str]:
    return list(results.keys())


def clustering_agreement_matrices(results: dict[str, ReconstructionResult]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return ARI and NMI matrices for the reconstructed parent partitions."""
    names = _names(results)
    ari = np.eye(len(names), dtype=float)
    nmi = np.eye(len(names), dtype=float)
    for i, a in enumerate(names):
        for j in range(i + 1, len(names)):
            b = names[j]
            la = np.asarray(results[a].parent_ids, int)
            lb = np.asarray(results[b].parent_ids, int)
            if len(la) != len(lb):
                raise ValueError("All compared reconstructions must contain the same daughter grains in the same order.")
            ari_v = float(adjusted_rand_score(la, lb))
            nmi_v = float(normalized_mutual_info_score(la, lb))
            ari[i, j] = ari[j, i] = ari_v
            nmi[i, j] = nmi[j, i] = nmi_v
    return pd.DataFrame(ari, index=names, columns=names), pd.DataFrame(nmi, index=names, columns=names)


def _boundary_flags(result: ReconstructionResult, edges: Iterable[tuple[int, int]]) -> np.ndarray:
    p = np.asarray(result.parent_ids, int)
    flags = []
    for i, j in edges:
        i, j = int(i), int(j)
        if not (0 <= i < len(p) and 0 <= j < len(p)):
            raise ValueError("Boundary comparison received an out-of-range adjacency index.")
        flags.append(bool(p[i] != p[j]))
    return np.asarray(flags, bool)


def boundary_jaccard_matrix(results: dict[str, ReconstructionResult], edges: Iterable[tuple[int, int]]) -> pd.DataFrame:
    """Jaccard agreement of reconstructed prior-parent boundary edges.

    Only experimentally supplied/constructed daughter-grain adjacency edges are considered.
    A value of 1 means the two methods mark exactly the same adjacency edges as prior-parent
    boundaries.  If neither method marks any boundary, the agreement is defined as 1.
    """
    edges = list(edges)
    names = _names(results)
    out = np.eye(len(names), dtype=float)
    flags = {n: _boundary_flags(results[n], edges) for n in names}
    for i, a in enumerate(names):
        for j in range(i + 1, len(names)):
            b = names[j]
            aa, bb = flags[a], flags[b]
            union = int(np.sum(aa | bb))
            score = 1.0 if union == 0 else float(np.sum(aa & bb) / union)
            out[i, j] = out[j, i] = score
    return pd.DataFrame(out, index=names, columns=names)


def boundary_consensus_table(
    results: dict[str, ReconstructionResult],
    edges: Iterable[tuple[int, int]],
    grain_ids: Iterable[int],
) -> pd.DataFrame:
    """One row per daughter adjacency edge with method votes for a prior-parent boundary."""
    edges = list(edges)
    gids = np.asarray(list(grain_ids), int)
    names = _names(results)
    rows: list[dict] = []
    for edge_id, (i, j) in enumerate(edges, 1):
        row = {
            "edge_id": edge_id,
            "daughter_grain_id_1": int(gids[int(i)]),
            "daughter_grain_id_2": int(gids[int(j)]),
        }
        votes = 0
        for name in names:
            is_boundary = bool(results[name].parent_ids[int(i)] != results[name].parent_ids[int(j)])
            row[f"boundary__{name}"] = is_boundary
            votes += int(is_boundary)
        row["methods_calling_parent_boundary"] = votes
        row["boundary_consensus_fraction"] = float(votes / len(names)) if names else np.nan
        row["consensus_label"] = (
            "unanimous boundary" if votes == len(names) and names else
            "unanimous same parent" if votes == 0 else
            "method disagreement"
        )
        rows.append(row)
    return pd.DataFrame(rows)


def _contingency(a: np.ndarray, b: np.ndarray) -> tuple[list[int], list[int], np.ndarray]:
    ua = sorted(set(int(x) for x in a.tolist()))
    ub = sorted(set(int(x) for x in b.tolist()))
    mat = np.zeros((len(ua), len(ub)), dtype=int)
    ia = {v: i for i, v in enumerate(ua)}
    ib = {v: i for i, v in enumerate(ub)}
    for x, y in zip(a, b):
        mat[ia[int(x)], ib[int(y)]] += 1
    return ua, ub, mat


def matched_parent_orientation_table(
    result_a: ReconstructionResult,
    result_b: ReconstructionResult,
    parent_sym: Iterable[np.ndarray],
    method_a: str = "method A",
    method_b: str = "method B",
) -> pd.DataFrame:
    """Match parent clusters by child-grain overlap, then compare their orientations.

    Hungarian assignment maximizes overlap.  This avoids comparing arbitrary parent ID
    numbers.  Rows with zero overlap are omitted.
    """
    a = np.asarray(result_a.parent_ids, int)
    b = np.asarray(result_b.parent_ids, int)
    if len(a) != len(b):
        raise ValueError("Compared reconstructions must contain the same daughter grains.")
    ua, ub, cont = _contingency(a, b)
    if cont.size == 0:
        return pd.DataFrame()
    r, c = linear_sum_assignment(-cont)
    rows = []
    for ri, ci in zip(r, c):
        overlap = int(cont[ri, ci])
        if overlap <= 0:
            continue
        pa, pb = ua[ri], ub[ci]
        ang = misorientation_deg(result_a.parent_orientations[pa], result_b.parent_orientations[pb], parent_sym)
        size_a = int(np.sum(a == pa))
        size_b = int(np.sum(b == pb))
        rows.append({
            "method_A": method_a,
            "method_B": method_b,
            "parent_id_A": pa,
            "parent_id_B": pb,
            "shared_daughter_grains": overlap,
            "fraction_of_parent_A_shared": float(overlap / size_a),
            "fraction_of_parent_B_shared": float(overlap / size_b),
            "parent_orientation_disagreement_deg": float(ang),
        })
    return pd.DataFrame(rows)


def matched_parent_orientation_matrix(
    results: dict[str, ReconstructionResult],
    parent_sym: Iterable[np.ndarray],
) -> pd.DataFrame:
    """Pairwise overlap-weighted mean misorientation of matched reconstructed parents."""
    names = _names(results)
    out = np.zeros((len(names), len(names)), dtype=float)
    for i, a in enumerate(names):
        for j in range(i + 1, len(names)):
            b = names[j]
            tbl = matched_parent_orientation_table(results[a], results[b], parent_sym, a, b)
            if tbl.empty:
                val = np.nan
            else:
                w = tbl["shared_daughter_grains"].to_numpy(float)
                x = tbl["parent_orientation_disagreement_deg"].to_numpy(float)
                val = float(np.average(x, weights=w))
            out[i, j] = out[j, i] = val
    return pd.DataFrame(out, index=names, columns=names)


def build_agreement_summary(
    results: dict[str, ReconstructionResult],
    edges: Iterable[tuple[int, int]],
    parent_sym: Iterable[np.ndarray],
) -> AgreementSummary:
    ari, nmi = clustering_agreement_matrices(results)
    return AgreementSummary(
        ari=ari,
        nmi=nmi,
        boundary_jaccard=boundary_jaccard_matrix(results, edges),
        matched_parent_orientation_deg=matched_parent_orientation_matrix(results, parent_sym),
    )


def variant_frequency_table(result: ReconstructionResult, source_df: pd.DataFrame | None = None) -> pd.DataFrame:
    """Global selected candidate-variant frequencies, optionally area-weighted."""
    df = result.table[["grain_id", "reconstructed_parent_id", "variant_candidate_id"]].copy()
    area_col = None
    if source_df is not None:
        for c in ("area_um2", "area_est_um2", "area"):
            if c in source_df.columns:
                area_col = c
                break
        keep = ["grain_id"] + ([area_col] if area_col else [])
        df = df.merge(source_df[keep], on="grain_id", how="left")
    rows = []
    total_n = len(df)
    total_area = float(pd.to_numeric(df[area_col], errors="coerce").fillna(0).sum()) if area_col else np.nan
    for vid, grp in df.groupby("variant_candidate_id", sort=True):
        row = {
            "candidate_variant_id": int(vid),
            "daughter_grains": int(len(grp)),
            "grain_fraction_pct": float(100.0 * len(grp) / total_n) if total_n else np.nan,
            "parents_containing_variant": int(grp["reconstructed_parent_id"].nunique()),
        }
        if area_col:
            area = float(pd.to_numeric(grp[area_col], errors="coerce").fillna(0).sum())
            row["daughter_area"] = area
            row["area_fraction_pct"] = float(100.0 * area / total_area) if total_area > 0 else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def operator_edge_table(
    child_orientations: list[np.ndarray],
    edges: Iterable[tuple[int, int]],
    grain_ids: Iterable[int],
    r_child_to_parent: np.ndarray,
    parent_sym: Iterable[np.ndarray],
    child_sym: Iterable[np.ndarray],
    operator_tolerance_deg: float,
    result: ReconstructionResult | None = None,
) -> pd.DataFrame:
    """ARPGE-style nearest theoretical operator classification on every neighbor edge."""
    edges = list(edges)
    gids = np.asarray(list(grain_ids), int)
    ops = theoretical_child_operators(r_child_to_parent, parent_sym, child_sym)
    rows = []
    for k, (i, j) in enumerate(edges, 1):
        i, j = int(i), int(j)
        measured = child_orientations[i].T @ child_orientations[j]
        residuals = [two_sided_operator_distance_deg(measured, op, child_sym) for op in ops]
        op_idx = int(np.argmin(residuals)) + 1
        residual = float(np.min(residuals))
        row = {
            "edge_id": k,
            "daughter_grain_id_1": int(gids[i]),
            "daughter_grain_id_2": int(gids[j]),
            "nearest_theoretical_operator_id": op_idx,
            "operator_residual_deg": residual,
            "within_operator_tolerance": bool(residual <= operator_tolerance_deg),
            "operator_tolerance_deg": float(operator_tolerance_deg),
        }
        if result is not None:
            row["reconstructed_same_parent"] = bool(result.parent_ids[i] == result.parent_ids[j])
            row["reconstructed_parent_id_1"] = int(result.parent_ids[i])
            row["reconstructed_parent_id_2"] = int(result.parent_ids[j])
        rows.append(row)
    return pd.DataFrame(rows)


def operator_frequency_table(edge_table: pd.DataFrame, accepted_only: bool = True) -> pd.DataFrame:
    if edge_table.empty:
        return pd.DataFrame(columns=["operator_id", "edges", "edge_fraction_pct"])
    df = edge_table
    if accepted_only and "within_operator_tolerance" in df.columns:
        df = df[df["within_operator_tolerance"]]
    total = len(df)
    rows = []
    for op, grp in df.groupby("nearest_theoretical_operator_id", sort=True):
        rows.append({
            "operator_id": int(op),
            "edges": int(len(grp)),
            "edge_fraction_pct": float(100.0 * len(grp) / total) if total else np.nan,
            "mean_operator_residual_deg": float(grp["operator_residual_deg"].mean()),
            "median_operator_residual_deg": float(grp["operator_residual_deg"].median()),
        })
    return pd.DataFrame(rows)


def parent_orientation_table(result: ReconstructionResult) -> pd.DataFrame:
    rows = []
    for pid, g in sorted(result.parent_orientations.items()):
        p1, P, p2 = matrix_to_bunge_euler(g)
        qw, qx, qy, qz = quaternion_wxyz(g)
        rows.append({
            "parent_id": int(pid),
            "Bunge_phi1_deg": p1,
            "Bunge_Phi_deg": P,
            "Bunge_phi2_deg": p2,
            "quaternion_qw": qw,
            "quaternion_qx": qx,
            "quaternion_qy": qy,
            "quaternion_qz": qz,
        })
    return pd.DataFrame(rows)


def academic_export_zip(
    *,
    source_df: pd.DataFrame,
    adjacency_df: pd.DataFrame,
    results: dict[str, ReconstructionResult],
    comparison: pd.DataFrame,
    agreement: AgreementSummary | None,
    parent_tables: dict[str, pd.DataFrame],
    daughter_tables: dict[str, pd.DataFrame],
    variant_tables: dict[str, pd.DataFrame],
    operator_tables: dict[str, pd.DataFrame] | None,
    metadata: dict,
    methods_text: str,
    boundary_consensus: pd.DataFrame | None = None,
) -> bytes:
    """Create a self-contained reconstruction evidence bundle for supplementary files."""
    bio = BytesIO()
    with ZipFile(bio, "w", ZIP_DEFLATED) as z:
        z.writestr("README.txt", methods_text.strip() + "\n")
        z.writestr("metadata.json", json.dumps(metadata, indent=2, default=str))
        z.writestr("input/daughter_grains.csv", source_df.to_csv(index=False))
        z.writestr("input/daughter_adjacency.csv", adjacency_df.to_csv(index=False))
        z.writestr("comparison/method_comparison.csv", comparison.to_csv(index=False))
        if agreement is not None:
            z.writestr("comparison/ARI_clustering_agreement.csv", agreement.ari.to_csv())
            z.writestr("comparison/NMI_clustering_agreement.csv", agreement.nmi.to_csv())
            z.writestr("comparison/boundary_Jaccard_agreement.csv", agreement.boundary_jaccard.to_csv())
            z.writestr("comparison/matched_parent_orientation_disagreement_deg.csv", agreement.matched_parent_orientation_deg.to_csv())
        if boundary_consensus is not None:
            z.writestr("comparison/prior_parent_boundary_consensus.csv", boundary_consensus.to_csv(index=False))
        for method in results:
            safe = "".join(c if c.isalnum() else "_" for c in method).strip("_")
            z.writestr(f"methods/{safe}/parent_summary.csv", parent_tables[method].to_csv(index=False))
            z.writestr(f"methods/{safe}/daughter_assignments.csv", daughter_tables[method].to_csv(index=False))
            z.writestr(f"methods/{safe}/variant_frequency.csv", variant_tables[method].to_csv(index=False))
            z.writestr(f"methods/{safe}/diagnostics.json", json.dumps(results[method].diagnostics, indent=2, default=float))
            z.writestr(f"methods/{safe}/parent_orientations.csv", parent_orientation_table(results[method]).to_csv(index=False))
            if operator_tables and method in operator_tables:
                z.writestr(f"methods/{safe}/operator_edges.csv", operator_tables[method].to_csv(index=False))
                z.writestr(f"methods/{safe}/operator_frequency.csv", operator_frequency_table(operator_tables[method]).to_csv(index=False))
    return bio.getvalue()
