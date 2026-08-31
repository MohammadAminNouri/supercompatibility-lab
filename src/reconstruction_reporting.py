from __future__ import annotations

"""Reporting, validation and method-selection helpers for parent/daughter reconstruction.

This module deliberately separates *reconstruction mathematics* (src.reconstruction)
from *research reporting*.  Scores reported here are descriptive diagnostics unless an
algorithm explicitly defines a probability.  In particular ``confidence`` produced by
src.reconstruction is a heuristic support score, not a calibrated probability.
"""

from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd

from .reconstruction import ReconstructionResult, matrix_to_bunge_euler, misorientation_deg, quaternion_wxyz


METHOD_REQUIREMENTS = {
    "Neighbor voting": {
        "minimum": "daughter grain orientations + neighbor graph + initial OR + parent/daughter symmetry",
        "uses_neighbors": True,
        "main_controls": "OR-fit angular width; parent-merge tolerance",
        "best_for": "fast local consensus and a transparent baseline",
        "primary_output": "parent assignment per daughter grain + local support",
        "limitation": "can propagate local mistakes and depends strongly on neighbor topology",
    },
    "Grain graph + Markov clustering": {
        "minimum": "daughter grain orientations + neighbor graph + initial OR + parent/daughter symmetry",
        "uses_neighbors": True,
        "main_controls": "OR-fit angular width; graph inflation",
        "best_for": "partitioning a connected daughter-grain network into prior-parent clusters",
        "primary_output": "parent clusters + per-grain fit/support",
        "limitation": "cluster granularity depends on graph quality and inflation",
    },
    "Variant graph probability propagation": {
        "minimum": "daughter grain orientations + neighbor graph + initial OR + parent/daughter symmetry",
        "uses_neighbors": True,
        "main_controls": "candidate angular width; graph inflation",
        "best_for": "ambiguous cases where each daughter grain has several plausible parent candidates",
        "primary_output": "candidate-aware parent clusters + selected variant candidate",
        "limitation": "computationally heavier; scores are algorithmic support, not calibrated probabilities",
    },
    "Nucleation + growth": {
        "minimum": "daughter grain orientations + neighbor graph + initial OR + parent/daughter symmetry",
        "uses_neighbors": True,
        "main_controls": "strict nucleation tolerance; looser growth tolerance",
        "best_for": "maps containing high-confidence seeds that can be expanded spatially",
        "primary_output": "seeded parent regions + grown daughter assignments",
        "limitation": "weak or sparse seeds can under-reconstruct; loose growth can over-grow",
    },
    "Operator / groupoid consistency": {
        "minimum": "daughter grain orientations + neighbor graph + initial OR + parent/daughter symmetry",
        "uses_neighbors": True,
        "main_controls": "daughter–daughter operator tolerance; common-parent consistency tolerance",
        "best_for": "explicit crystallographic checking of daughter–daughter transformation operators",
        "primary_output": "operator-consistent parent groups + fit/support",
        "limitation": "sensitive to OR accuracy and angular tolerance; not a binary clone of ARPGE",
    },
}


@dataclass(frozen=True)
class InputAudit:
    n_grains: int
    n_edges: int
    has_xy: bool
    has_area: bool
    has_phase: bool
    graph_kind: str
    connected_fraction: float
    isolated_grains: int
    warnings: tuple[str, ...]


def method_requirements_table() -> pd.DataFrame:
    rows = []
    for method, d in METHOD_REQUIREMENTS.items():
        rows.append({
            "method": method,
            "minimum input": d["minimum"],
            "main controls": d["main_controls"],
            "best used for": d["best_for"],
            "main output": d["primary_output"],
            "important limitation": d["limitation"],
        })
    return pd.DataFrame(rows)


def audit_reconstruction_input(df: pd.DataFrame, edges: Iterable[tuple[int, int]], graph_kind: str) -> InputAudit:
    n = len(df)
    edges = list(edges)
    deg = np.zeros(n, dtype=int)
    if "grain_id" in df.columns:
        ids = list(df["grain_id"].astype(int))
        idx = {g: i for i, g in enumerate(ids)}
        # edges may be index pairs or grain-id pairs; infer safely.
        for a, b in edges:
            if a in idx and b in idx:
                ia, ib = idx[a], idx[b]
            elif 0 <= int(a) < n and 0 <= int(b) < n:
                ia, ib = int(a), int(b)
            else:
                continue
            if ia != ib:
                deg[ia] += 1
                deg[ib] += 1
    connected = float(np.mean(deg > 0)) if n else 0.0
    isolated = int(np.sum(deg == 0)) if n else 0
    warnings: list[str] = []
    if n < 2:
        warnings.append("At least two daughter grains are needed for neighbor-based parent reconstruction.")
    if not edges:
        warnings.append("No neighbor graph is available; the five implemented reconstruction methods cannot use spatial consensus.")
    if graph_kind.lower().startswith("approx"):
        warnings.append("Neighbor topology is approximate (centroid k-NN), not a measured/shared-boundary graph; report this explicitly.")
    if isolated:
        warnings.append(f"{isolated} daughter grain(s) have no graph neighbor and therefore receive weak/no spatial support.")
    if "area" in df.columns and pd.to_numeric(df["area"], errors="coerce").isna().any():
        warnings.append("Some grain-area values are non-numeric; area-weighted parent summaries will be incomplete.")
    return InputAudit(
        n_grains=n,
        n_edges=len(edges),
        has_xy={"x", "y"}.issubset(df.columns),
        has_area="area" in df.columns,
        has_phase="phase" in df.columns,
        graph_kind=graph_kind,
        connected_fraction=connected,
        isolated_grains=isolated,
        warnings=tuple(warnings),
    )


def recommend_methods(audit: InputAudit) -> pd.DataFrame:
    rows = []
    for method, d in METHOD_REQUIREMENTS.items():
        if audit.n_edges == 0:
            status, reason = "not ready", "requires a daughter-grain neighbor graph"
        elif audit.n_grains < 3:
            status, reason = "weak evidence", "too few daughter grains for robust spatial consensus"
        elif audit.graph_kind.lower().startswith("approx") and method in {
            "Grain graph + Markov clustering", "Variant graph probability propagation", "Operator / groupoid consistency"
        }:
            status, reason = "usable with caution", "graph is approximate; true shared-boundary adjacency is preferred"
        else:
            status, reason = "ready", d["best_for"]
        rows.append({"method": method, "readiness": status, "why": reason})
    return pd.DataFrame(rows)


def fit_quality_label(fit_deg: float, excellent_deg: float = 1.0, acceptable_deg: float = 2.5, weak_deg: float = 5.0) -> str:
    x = float(fit_deg)
    if x <= excellent_deg:
        return "excellent"
    if x <= acceptable_deg:
        return "acceptable"
    if x <= weak_deg:
        return "weak"
    return "poor"


def support_quality_label(score: float, high: float = 0.75, medium: float = 0.45) -> str:
    s = float(score)
    if s >= high:
        return "high"
    if s >= medium:
        return "medium"
    return "low"


def daughter_assignment_table(
    result: ReconstructionResult,
    source_df: pd.DataFrame,
    excellent_deg: float = 1.0,
    acceptable_deg: float = 2.5,
    weak_deg: float = 5.0,
) -> pd.DataFrame:
    out = source_df.copy().merge(result.table, on="grain_id", how="left")
    out["fit_quality"] = [fit_quality_label(v, excellent_deg, acceptable_deg, weak_deg) for v in out["fit_deg"]]
    out["support_quality"] = [support_quality_label(v) for v in out["confidence"]]
    rename = {
        "reconstructed_parent_id": "parent_id",
        "variant_candidate_id": "daughter_variant_candidate_id",
        "fit_deg": "OR_fit_misorientation_deg",
        "confidence": "support_score_0_to_1",
        "parent_phi1_deg": "parent_Bunge_phi1_deg",
        "parent_Phi_deg": "parent_Bunge_Phi_deg",
        "parent_phi2_deg": "parent_Bunge_phi2_deg",
        "parent_qw": "parent_quaternion_qw",
        "parent_qx": "parent_quaternion_qx",
        "parent_qy": "parent_quaternion_qy",
        "parent_qz": "parent_quaternion_qz",
    }
    return out.rename(columns=rename)


def parent_summary_table(result: ReconstructionResult, source_df: pd.DataFrame) -> pd.DataFrame:
    merged = source_df.copy().merge(result.table, on="grain_id", how="left")
    rows = []
    for pid, group in merged.groupby("reconstructed_parent_id", sort=True):
        mat = result.parent_orientations[int(pid)]
        p1, P, p2 = matrix_to_bunge_euler(mat)
        qw, qx, qy, qz = quaternion_wxyz(mat)
        fit = pd.to_numeric(group["fit_deg"], errors="coerce")
        sup = pd.to_numeric(group["confidence"], errors="coerce")
        row = {
            "parent_id": int(pid),
            "supporting_daughter_grains": int(len(group)),
            "mean_OR_fit_deg": float(fit.mean()),
            "median_OR_fit_deg": float(fit.median()),
            "max_OR_fit_deg": float(fit.max()),
            "mean_support_score": float(sup.mean()),
            "minimum_support_score": float(sup.min()),
            "parent_Bunge_phi1_deg": p1,
            "parent_Bunge_Phi_deg": P,
            "parent_Bunge_phi2_deg": p2,
            "parent_quaternion_qw": qw,
            "parent_quaternion_qx": qx,
            "parent_quaternion_qy": qy,
            "parent_quaternion_qz": qz,
        }
        if "area" in group.columns:
            area = pd.to_numeric(group["area"], errors="coerce")
            row["summed_daughter_area"] = float(area.sum(min_count=1)) if area.notna().any() else np.nan
        if {"x", "y"}.issubset(group.columns):
            x = pd.to_numeric(group["x"], errors="coerce")
            y = pd.to_numeric(group["y"], errors="coerce")
            if "area" in group.columns:
                w = pd.to_numeric(group["area"], errors="coerce").fillna(0).clip(lower=0).to_numpy(float)
                if np.sum(w) > 0:
                    row["parent_centroid_x"] = float(np.average(x.to_numpy(float), weights=w))
                    row["parent_centroid_y"] = float(np.average(y.to_numpy(float), weights=w))
                else:
                    row["parent_centroid_x"] = float(x.mean())
                    row["parent_centroid_y"] = float(y.mean())
            else:
                row["parent_centroid_x"] = float(x.mean())
                row["parent_centroid_y"] = float(y.mean())
        rows.append(row)
    return pd.DataFrame(rows)


def method_comparison_row(result: ReconstructionResult, runtime_s: float, source_df: pd.DataFrame | None = None, validation_accuracy: float | None = None) -> dict:
    fit = np.asarray(result.fit_deg, float)
    support = np.asarray(result.confidence, float)
    row = {
        "method": result.method,
        "parent_grains_reconstructed": int(result.diagnostics["n_reconstructed_parents"]),
        "mean_OR_fit_deg": float(np.mean(fit)),
        "median_OR_fit_deg": float(np.median(fit)),
        "P95_OR_fit_deg": float(np.percentile(fit, 95)),
        "poor_fit_fraction_gt5deg": float(np.mean(fit > 5.0)),
        "mean_support_score_0_to_1": float(np.mean(support)),
        "low_support_fraction_lt0_45": float(np.mean(support < 0.45)),
        "runtime_s": float(runtime_s),
    }
    if validation_accuracy is not None:
        row["known_truth_clustering_accuracy"] = float(validation_accuracy)
    return row


def reconstruction_quality_summary(result: ReconstructionResult) -> dict[str, float | int | str]:
    fit = np.asarray(result.fit_deg, float)
    support = np.asarray(result.confidence, float)
    return {
        "daughter_grains": int(len(fit)),
        "reconstructed_parents": int(result.diagnostics["n_reconstructed_parents"]),
        "mean_fit_deg": float(np.mean(fit)),
        "p95_fit_deg": float(np.percentile(fit, 95)),
        "fraction_fit_le_2_5deg": float(np.mean(fit <= 2.5)),
        "fraction_fit_gt_5deg": float(np.mean(fit > 5.0)),
        "mean_support": float(np.mean(support)),
        "fraction_support_lt_0_45": float(np.mean(support < 0.45)),
        "interpretation": (
            "strong internal consistency" if np.mean(fit <= 2.5) >= 0.9 and np.mean(support < 0.45) <= 0.1
            else "mixed consistency; inspect weak grains and boundaries"
        ),
    }


def retained_parent_anchor_table(result: ReconstructionResult, anchors: dict[int, np.ndarray], parent_sym: Iterable[np.ndarray]) -> pd.DataFrame:
    """Compare reconstructed parents with independently measured retained-parent anchors.

    This is a *validation* table.  Anchors do not constrain the reconstruction in the
    current implementation, so the reported angular distance remains independent evidence.
    """
    rows = []
    for pid, gp in sorted(result.parent_orientations.items()):
        candidates = [(aid, misorientation_deg(gp, ga, parent_sym)) for aid, ga in anchors.items()]
        if not candidates:
            continue
        aid, ang = min(candidates, key=lambda x: x[1])
        rows.append({
            "reconstructed_parent_id": int(pid),
            "nearest_retained_parent_anchor_id": int(aid),
            "anchor_misorientation_deg": float(ang),
        })
    return pd.DataFrame(rows)


COLUMN_MEANINGS = {
    "grain_id": "Unique daughter-grain identifier. Dimensionless label; it has no crystallographic meaning.",
    "x": "Daughter-grain centroid x coordinate in the map. Use the same length unit for x and y, normally µm.",
    "y": "Daughter-grain centroid y coordinate in the map. Use the same length unit for x and y, normally µm.",
    "area": "Measured daughter-grain area, normally µm². Optional; used only for area-weighted summaries.",
    "phase": "Measured daughter phase label/name. Optional metadata; the reconstruction symmetry is selected separately.",
    "phi1_deg": "First Bunge Euler angle φ₁ of the measured daughter orientation, in degrees.",
    "Phi_deg": "Second Bunge Euler angle Φ of the measured daughter orientation, in degrees.",
    "phi2_deg": "Third Bunge Euler angle φ₂ of the measured daughter orientation, in degrees.",
    "qw": "Scalar component w of a unit quaternion describing daughter orientation; dimensionless.",
    "qx": "Quaternion x component; dimensionless.",
    "qy": "Quaternion y component; dimensionless.",
    "qz": "Quaternion z component; dimensionless.",
    "parent_id": "Reconstructed parent-grain cluster identifier. Arbitrary label; compare orientations/regions, not the numeric ID itself.",
    "daughter_variant_candidate_id": "Index of the symmetry-related daughter/parent candidate selected for that daughter grain. Dimensionless label.",
    "OR_fit_misorientation_deg": "Angular mismatch between the daughter grain and the orientation predicted from its reconstructed parent through the selected/refined OR, in degrees. Smaller is better.",
    "support_score_0_to_1": "Heuristic algorithmic support score from 0 to 1. It is NOT a calibrated probability.",
    "fit_quality": "Human-readable class derived from user-visible angular thresholds; descriptive, not a theorem.",
    "support_quality": "Human-readable class derived from the heuristic support score; descriptive, not a probability statement.",
    "parent_Bunge_phi1_deg": "First Bunge Euler angle φ₁ of the reconstructed parent orientation, in degrees.",
    "parent_Bunge_Phi_deg": "Second Bunge Euler angle Φ of the reconstructed parent orientation, in degrees.",
    "parent_Bunge_phi2_deg": "Third Bunge Euler angle φ₂ of the reconstructed parent orientation, in degrees.",
}


def column_dictionary(columns: Iterable[str]) -> pd.DataFrame:
    rows = []
    for c in columns:
        rows.append({"column": c, "meaning": COLUMN_MEANINGS.get(c, "See the local table caption/documentation for this derived field.")})
    return pd.DataFrame(rows)
