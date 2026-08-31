from __future__ import annotations

"""High-clarity Streamlit workbench for parent/daughter orientation reconstruction.

This module wraps the numerical reconstruction engine in :mod:`src.reconstruction`.
It deliberately keeps four things separate:

1. experimental input / segmentation,
2. crystallographic orientation relationship (OR),
3. reconstruction algorithm,
4. interpretation / confidence of the output.

The goal is to make every input, output and threshold explicit enough to report in
an academic methods section.  Heuristic confidence values are always labelled as
heuristic support scores rather than probabilities.
"""

from dataclasses import dataclass
from io import StringIO
from pathlib import Path
import hashlib
import json
import time
from typing import Iterable

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation
from sklearn.metrics import adjusted_rand_score

from src.reconstruction import (
    ORPreset,
    ReconstructionResult,
    approximate_knn_edges,
    bunge_euler_to_matrix,
    edges_from_dataframe,
    grain_graph_reconstruction,
    matrix_from_quaternion_wxyz,
    matrix_to_bunge_euler,
    misorientation_deg,
    neighbor_voting_reconstruction,
    nucleation_growth_reconstruction,
    operator_groupoid_reconstruction,
    orientation_relationship_presets,
    orientations_from_dataframe,
    quaternion_wxyz,
    refine_orientation_relationship,
    rotation_from_parallelisms,
    symmetry_group,
    synthetic_parent_reconstruction_demo,
    synthetic_parent_reconstruction_demo_from_or,
    unique_child_variants,
    variant_graph_reconstruction,
)


from src.retransformation import (
    NiTiAQLattice,
    branch_switch_catalog,
    cycle_closure_table,
    cycle_evidence_zip,
    match_new_cycle_daughters,
    niti_aq_orientation_relationship,
    niti_aq_parallelism_residuals,
    niti_ct_otsuka_ren_diagnostics,
    niti_ct_otsuka_ren_orientation_relationship,
    observed_branch_occupancy,
    parent_cycle_summary,
    regenerated_variant_library,
)
from src.core import C_A_TO_M, C_M_TO_A

from src.reconstruction_academic import (
    academic_export_zip,
    known_truth_validation_metrics,
    boundary_consensus_table,
    build_agreement_summary,
    operator_edge_table,
    operator_frequency_table,
    variant_frequency_table,
)
from src.reconstruction_reporting import (
    audit_reconstruction_input,
    column_dictionary,
    daughter_assignment_table,
    method_requirements_table,
    parent_summary_table,
    recommend_methods,
    reconstruction_quality_summary,
)


METHODS = [
    "Neighbor voting",
    "Grain graph + Markov clustering",
    "Variant graph probability propagation",
    "Nucleation + growth",
    "Operator / groupoid consistency",
]

SYMMETRIES = ["cubic", "hexagonal", "orthorhombic", "monoclinic", "triclinic"]


@dataclass(frozen=True)
class PreparedDataset:
    name: str
    grains: pd.DataFrame
    adjacency: pd.DataFrame
    convention: str
    provenance: str
    graph_kind: str = "measured/shared-boundary"
    length_unit: str = "map unit"
    input_notes: str = ""


# ---------------------------------------------------------------------------
# Small pure helpers (also useful in tests / batch scripts)
# ---------------------------------------------------------------------------

def _quality_label(fit_deg: float, support: float) -> str:
    """Human-readable heuristic quality band; not a calibrated probability."""
    if fit_deg <= 1.0 and support >= 0.80:
        return "Very strong"
    if fit_deg <= 2.5 and support >= 0.60:
        return "Strong"
    if fit_deg <= 5.0 and support >= 0.35:
        return "Review"
    return "Weak / ambiguous"


def _method_explanation(name: str) -> tuple[str, str, str]:
    if name == "Neighbor voting":
        return (
            "Local consensus",
            "Each daughter grain compares its candidate parent orientations with candidates supported by neighboring daughter grains.",
            "Needs a trustworthy neighbor graph; best for locally coherent parent regions.",
        )
    if name == "Grain graph + Markov clustering":
        return (
            "Grain graph clustering",
            "Daughter grains are graph nodes; edge weights encode whether neighboring grains can share a crystallographically consistent parent.",
            "Sensitive to graph topology and the Markov inflation parameter.",
        )
    if name == "Variant graph probability propagation":
        return (
            "Candidate-variant graph",
            "Each grain keeps several possible parent candidates; compatibility is propagated between candidate nodes before a parent is selected.",
            "Usually the richest graph model here, but also the most expensive.",
        )
    if name == "Nucleation + growth":
        return (
            "Seed then grow",
            "Strict high-confidence parent nuclei are identified first, then expanded to neighboring daughter grains using a looser growth tolerance.",
            "Useful when reliable local seeds exist; results depend on nucleation and growth thresholds.",
        )
    return (
        "Transformation-operator consistency",
        "Measured daughter-to-daughter misorientations are checked against theoretical transformation operators, followed by common-parent consistency.",
        "Most explicitly crystallographic of the five routes; operator and parent-consistency tolerances must be reported.",
    )


def _parse_vec(text: str, label: str) -> np.ndarray:
    vals = [float(x) for x in text.replace(",", " ").replace("[", " ").replace("]", " ").replace("(", " ").replace(")", " ").split()]
    if len(vals) != 3:
        raise ValueError(f"{label} must contain exactly three numbers, e.g. 1 1 1.")
    a = np.asarray(vals, float)
    if np.linalg.norm(a) <= 1e-14:
        raise ValueError(f"{label} cannot be the zero vector.")
    return a


def _parse_matrix(text: str, label: str = "OR matrix") -> np.ndarray:
    rows = []
    for row in text.strip().split(";"):
        if row.strip():
            rows.append([float(x) for x in row.replace(",", " ").split()])
    a = np.asarray(rows, float)
    if a.shape != (3, 3):
        raise ValueError(f"{label} must have exactly 3 rows × 3 columns. Separate rows with semicolons.")
    u, _, vt = np.linalg.svd(a)
    r = u @ vt
    if np.linalg.det(r) < 0:
        u[:, -1] *= -1
        r = u @ vt
    return r


def _coerce_grain_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Map common export aliases to the names used by the reconstruction engine."""
    out = df.copy()
    aliases = {
        "grain": "grain_id",
        "grainid": "grain_id",
        "grain_id": "grain_id",
        "phi1": "phi1_deg",
        "phi_1": "phi1_deg",
        "phi1_deg": "phi1_deg",
        "phi": "Phi_deg",
        "phi_deg": "Phi_deg",
        "phi2": "phi2_deg",
        "phi_2": "phi2_deg",
        "phi2_deg": "phi2_deg",
        "area": "area_um2",
        "area_um2": "area_um2",
        "phaseid": "phase",
    }
    rename = {}
    for c in out.columns:
        key = str(c).strip().lower().replace(" ", "").replace("-", "_")
        if key in aliases:
            rename[c] = aliases[key]
    out = out.rename(columns=rename)
    return out


def _read_ang_bytes(raw: bytes) -> pd.DataFrame:
    """Read a conventional EDAX/TSL ``.ang`` text map conservatively.

    Standard ANG v5 columns are ``phi1 PHI phi2 x y IQ CI phase SEM fit`` and the
    Euler angles are radians.  We preserve the quality columns when present.  For a
    non-standard file whose Euler values are clearly outside a radian range, we keep
    the data but mark that degrees were inferred in ``DataFrame.attrs`` instead of
    silently applying a radian conversion.
    """
    text = raw.decode("utf-8", errors="replace")
    header = [ln.strip() for ln in text.splitlines() if ln.strip().startswith("#")]
    rows: list[list[float]] = []
    widths: list[int] = []
    for line_no, line in enumerate(text.splitlines(), 1):
        sline = line.strip()
        if not sline or sline.startswith("#"):
            continue
        parts = sline.replace(",", " ").split()
        try:
            vals = [float(x) for x in parts]
        except ValueError:
            continue
        if len(vals) >= 8:
            rows.append(vals)
            widths.append(len(vals))
    if not rows:
        raise ValueError("No standard numeric ANG rows with at least 8 columns were found.")
    minw = min(widths)
    if minw < 8:
        raise ValueError("ANG rows do not contain the required phi1/PHI/phi2/x/y/IQ/CI/phase fields.")
    maxw = min(max(widths), 10)
    arr = np.full((len(rows), maxw), np.nan, dtype=float)
    for i, vals in enumerate(rows):
        arr[i, : min(len(vals), maxw)] = vals[:maxw]

    e = arr[:, :3]
    finite_e = e[np.isfinite(e)]
    if finite_e.size == 0:
        raise ValueError("ANG Euler-angle columns are empty/non-finite.")
    max_abs = float(np.max(np.abs(finite_e)))
    # EDAX ANG is radians.  Values clearly beyond 2*pi are treated as a non-standard
    # degree export; ambiguous cases stay on the documented ANG-radian route.
    inferred_degrees = max_abs > (2.0 * np.pi + 0.25)
    if inferred_degrees:
        euler_deg = e.copy()
        angle_note = "Non-standard ANG: Euler values exceeded 2π, so degrees were inferred. Verify against the acquisition software before publication."
    else:
        euler_deg = np.degrees(e)
        angle_note = "Standard ANG convention: Euler angles interpreted as radians and converted to degrees."

    out = pd.DataFrame({
        "phi1_deg": euler_deg[:, 0],
        "Phi_deg": euler_deg[:, 1],
        "phi2_deg": euler_deg[:, 2],
        "x": arr[:, 3],
        "y": arr[:, 4],
        "IQ": arr[:, 5],
        "CI": arr[:, 6],
        "phase": arr[:, 7],
    })
    if maxw >= 9:
        out["SEM_signal"] = arr[:, 8]
    if maxw >= 10:
        out["fit"] = arr[:, 9]
    # Phase 0 is conventionally non-indexed.  Do not remove it here: the UI shows the
    # indexed fraction and lets the user select/filter the daughter phase explicitly.
    out.attrs.update({
        "format": "ANG",
        "source_header": "\n".join(header[:80]),
        "angle_note": angle_note,
        "euler_input_unit": "degrees (inferred)" if inferred_degrees else "radians (standard ANG)",
        "standard_columns": "phi1 PHI phi2 x y IQ CI phase SEM fit",
    })
    return out


def _read_ctf_bytes(raw: bytes) -> pd.DataFrame:
    """Read Oxford/HKL Channel Text File point data with explicit field validation."""
    text = raw.decode("utf-8", errors="replace")
    lines = text.splitlines()
    header_idx = None
    for i, line in enumerate(lines):
        low = line.strip().lower()
        if low.startswith("phase") and all(k in low for k in ("x", "y", "euler1", "euler2", "euler3")):
            header_idx = i
            break
    if header_idx is None:
        raise ValueError("Could not find a CTF point-data header containing Phase, X, Y, Euler1, Euler2 and Euler3.")

    block = "\n".join(lines[header_idx:])
    # Standard CTF is tab delimited.  Use tab first to avoid accidental splitting of
    # textual header fields, then fall back to generic whitespace for permissive exports.
    try:
        data = pd.read_csv(StringIO(block), sep="\t")
        if len(data.columns) < 6:
            raise ValueError
    except Exception:
        data = pd.read_csv(StringIO(block), sep=r"\s+", engine="python")
    lower = {str(c).strip().lower(): c for c in data.columns}
    needed = ["phase", "x", "y", "euler1", "euler2", "euler3"]
    missing = [c for c in needed if c not in lower]
    if missing:
        raise ValueError(f"CTF file is missing expected columns: {', '.join(missing)}")

    out = pd.DataFrame({
        "phase": pd.to_numeric(data[lower["phase"]], errors="coerce"),
        "x": pd.to_numeric(data[lower["x"]], errors="coerce"),
        "y": pd.to_numeric(data[lower["y"]], errors="coerce"),
        "phi1_deg": pd.to_numeric(data[lower["euler1"]], errors="coerce"),
        "Phi_deg": pd.to_numeric(data[lower["euler2"]], errors="coerce"),
        "phi2_deg": pd.to_numeric(data[lower["euler3"]], errors="coerce"),
    })
    optional = {
        "bands": "Bands", "error": "Error", "mad": "MAD", "bc": "BC", "bs": "BS"
    }
    for key, nice in optional.items():
        if key in lower:
            out[nice] = pd.to_numeric(data[lower[key]], errors="coerce")
    out = out.dropna(subset=["phase", "x", "y", "phi1_deg", "Phi_deg", "phi2_deg"]).reset_index(drop=True)
    if out.empty:
        raise ValueError("CTF point table was found, but no complete numeric orientation rows remained after parsing.")

    header_meta: dict[str, str] = {}
    for line in lines[:header_idx]:
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            header_meta[parts[0].strip()] = parts[1].strip()
    out.attrs.update({
        "format": "CTF",
        "source_header": "\n".join(lines[:header_idx][:120]),
        "angle_note": "Standard CTF convention: Euler1/Euler2/Euler3 interpreted as degrees.",
        "euler_input_unit": "degrees",
        "header_metadata": header_meta,
    })
    return out


def _raw_ebsd_audit(df: pd.DataFrame) -> pd.DataFrame:
    """Compact import audit shown before segmentation; never silently 'fixes' the data."""
    rows: list[dict[str, object]] = []
    n = len(df)
    phase = pd.to_numeric(df.get("phase", pd.Series(np.ones(n))), errors="coerce")
    indexed = phase.fillna(0) > 0
    rows.extend([
        {"check": "rows parsed", "value": n, "interpretation": "EBSD measurements read from the file"},
        {"check": "indexed rows", "value": int(indexed.sum()), "interpretation": "phase > 0"},
        {"check": "indexed fraction", "value": f"{100*float(indexed.mean()):.2f}%" if n else "n/a", "interpretation": "non-indexed phase 0 is excluded before daughter segmentation"},
        {"check": "Euler input convention", "value": df.attrs.get("euler_input_unit", "unknown"), "interpretation": df.attrs.get("angle_note", "verify orientation convention")},
    ])
    if {"x", "y"}.issubset(df.columns) and n > 1:
        xy = df[["x", "y"]].to_numpy(float)
        good = np.isfinite(xy).all(axis=1)
        if np.sum(good) > 1:
            tree = cKDTree(xy[good])
            d, _ = tree.query(xy[good], k=2)
            nn = d[:, 1]
            nn = nn[np.isfinite(nn) & (nn > 0)]
            if len(nn):
                rows.append({"check": "median nearest-point spacing", "value": f"{float(np.median(nn)):.6g}", "interpretation": "same length unit as x/y; used to infer local pixel connectivity"})
    for q in ("CI", "IQ", "MAD", "BC", "Bands", "fit"):
        if q in df.columns:
            vals = pd.to_numeric(df[q], errors="coerce")
            rows.append({"check": f"{q} available", "value": f"median={float(vals.median()):.4g}", "interpretation": "optional vendor quality field; filtering is explicit, never automatic"})
    return pd.DataFrame(rows)


def _read_uploaded_orientation_file(uploaded) -> tuple[pd.DataFrame, str]:
    name = uploaded.name
    raw = uploaded.getvalue()
    suffix = Path(name).suffix.lower()
    if suffix == ".ang":
        return _read_ang_bytes(raw), "raw .ang EBSD point table"
    if suffix == ".ctf":
        return _read_ctf_bytes(raw), "raw .ctf EBSD point table"
    if suffix in {".csv", ".txt", ".tsv"}:
        # pandas' separator inference handles ordinary comma/tab-delimited exports.
        try:
            df = pd.read_csv(StringIO(raw.decode("utf-8", errors="ignore")), sep=None, engine="python")
        except Exception:
            df = pd.read_csv(StringIO(raw.decode("utf-8", errors="ignore")))
        return _coerce_grain_columns(df), "uploaded delimited table"
    raise ValueError("Supported reconstruction uploads are .csv, .tsv, .txt, .ang and .ctf.")


def _segment_points_to_grains(
    points: pd.DataFrame,
    child_sym: Iterable[np.ndarray],
    misorientation_threshold_deg: float,
    min_points: int,
    max_points: int = 40000,
    neighbor_radius_factor: float = 1.35,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Auditable segmentation for modest raw EBSD maps.

    Pixels are connected only when they are spatial nearest-neighbors, belong to the same
    indexed phase, and differ by no more than the user-supplied disorientation threshold.
    The mean grain orientation is symmetry-aligned before averaging.  This is deliberately
    conservative: for production maps, vendor software/MTEX is still preferred, while this
    path is intended to be transparent and reproducible for small maps and cross-checks.
    """
    required = {"x", "y", "phi1_deg", "Phi_deg", "phi2_deg"}
    if not required.issubset(points.columns):
        raise ValueError("Raw-point segmentation requires x, y, phi1_deg, Phi_deg and phi2_deg columns.")
    df = points.dropna(subset=list(required)).reset_index(drop=True).copy()
    if "phase" in df.columns:
        df = df[pd.to_numeric(df["phase"], errors="coerce").fillna(0) > 0].reset_index(drop=True)
    if len(df) < 2:
        raise ValueError("At least two indexed EBSD points are required after phase/quality filtering.")
    if len(df) > max_points:
        raise ValueError(
            f"The auditable in-app segmenter is capped at {max_points:,} points; this selection has {len(df):,}. "
            "Segment the full map in a dedicated EBSD package and upload the grain table + adjacency instead."
        )
    if not (0.5 <= float(neighbor_radius_factor) <= 2.5):
        raise ValueError("neighbor_radius_factor must be between 0.5 and 2.5.")

    xy = df[["x", "y"]].to_numpy(float)
    if not np.isfinite(xy).all():
        raise ValueError("x/y coordinates must be finite.")
    tree = cKDTree(xy)
    d, _ = tree.query(xy, k=2)
    nn = d[:, 1]
    nn = nn[np.isfinite(nn) & (nn > 0)]
    if not len(nn):
        raise ValueError("Could not infer EBSD point spacing: coordinates appear duplicated or degenerate.")
    step = float(np.median(nn))
    # A small tolerance above the nominal nearest-neighbor spacing accepts square/hex
    # vendor grids without bridging the diagonal of a square grid by default.
    pairs = sorted(tree.query_pairs(r=float(neighbor_radius_factor) * step))
    if not pairs:
        raise ValueError("No spatial neighbor pairs were found. Check x/y units or increase the neighbor-radius factor slightly.")

    ori = [
        bunge_euler_to_matrix(float(a), float(b), float(c))
        for a, b, c in df[["phi1_deg", "Phi_deg", "phi2_deg"]].to_numpy()
    ]
    phase = df["phase"].to_numpy() if "phase" in df.columns else np.ones(len(df), int)
    syms = tuple(child_sym)

    parent = np.arange(len(df), dtype=int)
    size = np.ones(len(df), dtype=int)

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return int(x)

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra == rb:
            return
        if size[ra] < size[rb]:
            ra, rb = rb, ra
        parent[rb] = ra
        size[ra] += size[rb]

    for i, j in pairs:
        if phase[i] != phase[j]:
            continue
        if misorientation_deg(ori[i], ori[j], syms) <= float(misorientation_threshold_deg):
            union(i, j)

    root = np.array([find(i) for i in range(len(df))], int)
    counts = pd.Series(root).value_counts()
    keep_roots = set(int(r) for r, n in counts.items() if int(n) >= int(min_points))
    if not keep_roots:
        raise ValueError("No grains survived the minimum-point threshold. Lower the minimum grain size or inspect indexing/segmentation settings.")

    roots_sorted = sorted(keep_roots)
    rid_to_gid = {r: k + 1 for k, r in enumerate(roots_sorted)}
    point_gid = np.array([rid_to_gid.get(int(r), 0) for r in root], int)

    def symmetry_aligned_mean(indices: np.ndarray) -> tuple[np.ndarray, float, float]:
        mats = [ori[int(i)] for i in indices]
        ref = mats[0]
        aligned = []
        syma = np.asarray(syms, float)
        for g in mats:
            # Choose the right-acting crystal-symmetry equivalent of g with the largest
            # relative-rotation trace to the reference (equivalently the smallest angle).
            candidates = np.einsum("ij,kjl->kil", g, syma)
            rel = np.einsum("kij,lj->kil", candidates, ref)
            traces = np.trace(rel, axis1=1, axis2=2)
            aligned.append(candidates[int(np.argmax(traces))])
        gmean = Rotation.from_matrix(np.stack(aligned)).mean().as_matrix()
        spread = np.array([misorientation_deg(g, gmean, syms) for g in mats], float)
        return gmean, float(np.sqrt(np.mean(spread ** 2))), float(np.max(spread))

    grain_rows = []
    for r in roots_sorted:
        idx = np.flatnonzero(root == r)
        gmean, gos_rms, gos_max = symmetry_aligned_mean(idx)
        p1, P, p2 = matrix_to_bunge_euler(gmean)
        row = {
            "grain_id": rid_to_gid[r],
            "x": float(df.loc[idx, "x"].mean()),
            "y": float(df.loc[idx, "y"].mean()),
            "phi1_deg": p1,
            "Phi_deg": P,
            "phi2_deg": p2,
            "point_count": int(len(idx)),
            "area_est_mapunit2": float(len(idx) * step * step),
            "grain_orientation_spread_rms_deg": gos_rms,
            "grain_orientation_spread_max_deg": gos_max,
        }
        if "phase" in df.columns:
            row["phase"] = pd.Series(df.loc[idx, "phase"]).mode().iloc[0]
        grain_rows.append(row)
    grains = pd.DataFrame(grain_rows)

    contact_counts: dict[tuple[int, int], int] = {}
    for i, j in pairs:
        gi, gj = int(point_gid[i]), int(point_gid[j])
        if gi > 0 and gj > 0 and gi != gj:
            key = (min(gi, gj), max(gi, gj))
            contact_counts[key] = contact_counts.get(key, 0) + 1
    adjacency = pd.DataFrame([
        {
            "grain_id_1": a,
            "grain_id_2": b,
            "boundary_contact_count": n,
            "boundary_length_est_mapunit": float(n * step),
        }
        for (a, b), n in sorted(contact_counts.items())
    ])
    if adjacency.empty:
        adjacency = pd.DataFrame(columns=["grain_id_1", "grain_id_2", "boundary_contact_count", "boundary_length_est_mapunit"])
    grains.attrs.update({
        "estimated_point_spacing": step,
        "segmentation_threshold_deg": float(misorientation_threshold_deg),
        "minimum_points_per_grain": int(min_points),
        "neighbor_radius_factor": float(neighbor_radius_factor),
    })
    return grains, adjacency


def _parent_summary(result: ReconstructionResult, source_df: pd.DataFrame) -> pd.DataFrame:
    """One academically useful row per reconstructed parent grain."""
    merged = source_df.merge(result.table, on="grain_id", how="right", suffixes=("", "_recon"))
    rows = []
    area_col = next((c for c in ["area_um2", "area_est_mapunit2", "area_est_um2", "area"] if c in merged.columns), None)
    total_area = float(pd.to_numeric(merged[area_col], errors="coerce").fillna(0).sum()) if area_col else np.nan
    for pid, grp in merged.groupby("reconstructed_parent_id", sort=True):
        first = grp.iloc[0]
        fit = pd.to_numeric(grp["fit_deg"], errors="coerce")
        conf = pd.to_numeric(grp["confidence"], errors="coerce")
        gap = pd.to_numeric(grp.get("candidate_separation_deg", np.nan), errors="coerce")
        area = float(pd.to_numeric(grp[area_col], errors="coerce").fillna(0).sum()) if area_col else np.nan
        variants = pd.to_numeric(grp["variant_candidate_id"], errors="coerce").dropna().astype(int)
        vc = variants.value_counts()
        mean_fit = float(fit.mean())
        mean_conf = float(conf.mean())
        row = {
            "parent_id": int(pid),
            "supporting_daughter_grains": int(len(grp)),
            "distinct_selected_variants": int(variants.nunique()),
            "dominant_variant_id": int(vc.index[0]) if len(vc) else np.nan,
            "dominant_variant_fraction_pct": float(100.0 * vc.iloc[0] / len(variants)) if len(variants) else np.nan,
            "supporting_area": area if area_col else np.nan,
            "area_fraction_pct": (100.0 * area / total_area) if area_col and total_area > 0 else np.nan,
            "mean_OR_fit_deg": mean_fit,
            "median_OR_fit_deg": float(fit.median()),
            "p95_OR_fit_deg": float(fit.quantile(0.95)),
            "max_OR_fit_deg": float(fit.max()),
            "mean_candidate_separation_deg": float(gap.mean()) if hasattr(gap, "mean") else np.nan,
            "mean_support_score": mean_conf,
            "minimum_support_score": float(conf.min()),
            "weak_or_ambiguous_children": int(((fit > 5.0) | (conf < 0.35)).sum()),
            "parent_phi1_deg": float(first["parent_phi1_deg"]),
            "parent_Phi_deg": float(first["parent_Phi_deg"]),
            "parent_phi2_deg": float(first["parent_phi2_deg"]),
            "parent_qw": float(first["parent_qw"]),
            "parent_qx": float(first["parent_qx"]),
            "parent_qy": float(first["parent_qy"]),
            "parent_qz": float(first["parent_qz"]),
            "quality_flag": _quality_label(mean_fit, mean_conf),
        }
        if {"x", "y"}.issubset(grp.columns):
            x = pd.to_numeric(grp["x"], errors="coerce")
            y = pd.to_numeric(grp["y"], errors="coerce")
            if area_col:
                w = pd.to_numeric(grp[area_col], errors="coerce").fillna(0).clip(lower=0).to_numpy(float)
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


def _assignment_table(result: ReconstructionResult, source_df: pd.DataFrame) -> pd.DataFrame:
    """Per-daughter assignment with candidate ambiguity exposed instead of hidden."""
    show = source_df.merge(result.table, on="grain_id", how="right")
    show = show.rename(columns={
        "reconstructed_parent_id": "parent_id",
        "variant_candidate_id": "candidate_variant_id",
        "fit_deg": "best_OR_fit_deg",
        "second_best_candidate_fit_deg": "second_best_candidate_fit_deg",
        "candidate_separation_deg": "candidate_separation_deg",
        "candidate_count": "candidate_count",
        "absolute_fit_support": "absolute_fit_support_0_to_1",
        "separation_support": "candidate_separation_support_0_to_1",
        "confidence": "support_score_0_to_1",
    })
    show["quality_flag"] = [
        _quality_label(float(f), float(c))
        for f, c in zip(show["best_OR_fit_deg"], show["support_score_0_to_1"])
    ]
    preferred = [
        "grain_id", "parent_id", "candidate_variant_id", "best_OR_fit_deg",
        "second_best_candidate_fit_deg", "candidate_separation_deg", "candidate_count",
        "absolute_fit_support_0_to_1", "candidate_separation_support_0_to_1",
        "support_score_0_to_1", "quality_flag", "x", "y", "area_um2",
        "area_est_mapunit2", "area_est_um2", "point_count",
        "grain_orientation_spread_rms_deg", "grain_orientation_spread_max_deg",
        "phase", "parent_phi1_deg", "parent_Phi_deg", "parent_phi2_deg",
        "parent_qw", "parent_qx", "parent_qy", "parent_qz",
    ]
    cols = [c for c in preferred if c in show.columns] + [c for c in show.columns if c not in preferred]
    return show[cols]


def _comparison_table(
    results: dict[str, ReconstructionResult],
    runtimes: dict[str, float],
    edges: Iterable[tuple[int, int]] | None = None,
    truth: np.ndarray | None = None,
) -> pd.DataFrame:
    """Comparison table without collapsing different evidence into one opaque score."""
    edges = list(edges or [])
    rows = []
    for name, r in results.items():
        fit = np.asarray(r.fit_deg, float)
        conf = np.asarray(r.confidence, float)
        pids = np.asarray(r.parent_ids, int)
        _, counts = np.unique(pids, return_counts=True)
        singleton_parent_fraction = float(np.mean(counts == 1)) if len(counts) else np.nan
        boundary_count = int(sum(pids[int(i)] != pids[int(j)] for i, j in edges)) if edges else 0
        row = {
            "method": name,
            "reconstructed_parents": int(r.diagnostics["n_reconstructed_parents"]),
            "singleton_parent_fraction_pct": 100.0 * singleton_parent_fraction if np.isfinite(singleton_parent_fraction) else np.nan,
            "mean_OR_fit_deg": float(np.mean(fit)),
            "median_OR_fit_deg": float(np.median(fit)),
            "p95_OR_fit_deg": float(np.quantile(fit, 0.95)),
            "fraction_OR_fit_gt5deg_pct": float(100.0 * np.mean(fit > 5.0)),
            "mean_support_score": float(np.mean(conf)),
            "low_support_fraction_pct": float(100.0 * np.mean(conf < 0.35)),
            "prior_parent_boundary_edges": boundary_count,
            "prior_parent_boundary_fraction_pct": float(100.0 * boundary_count / len(edges)) if edges else np.nan,
            "runtime_s": float(runtimes.get(name, np.nan)),
        }
        if truth is not None:
            tv = known_truth_validation_metrics(r, truth, edges)
            row.update({
                "truth_status": tv["validation_status"],
                "expected_parents": tv["expected_parent_count"],
                "truth_ARI": tv["truth_ARI"],
                "truth_completeness": tv["truth_completeness"],
                "truth_homogeneity": tv["truth_homogeneity"],
                "truth_V_measure": tv["truth_V_measure"],
                "truth_boundary_precision": tv["truth_boundary_precision"],
                "truth_boundary_recall": tv["truth_boundary_recall"],
                "truth_boundary_F1": tv["truth_boundary_F1"],
                "false_parent_boundary_rate": tv["false_parent_boundary_rate"],
                "fragments_per_true_parent_mean": tv["mean_reconstructed_fragments_per_true_parent"],
            })
        rows.append(row)
    return pd.DataFrame(rows)



def _comparison_display_table(comp: pd.DataFrame) -> pd.DataFrame:
    """Human-readable UI labels; canonical machine columns remain in downloads."""
    rename = {
        "method": "Method",
        "truth_status": "Known-truth verdict",
        "expected_parents": "True parents",
        "reconstructed_parents": "Reconstructed parents",
        "singleton_parent_fraction_pct": "Singleton parents (%)",
        "mean_OR_fit_deg": "Mean OR residual (°)",
        "median_OR_fit_deg": "Median OR residual (°)",
        "p95_OR_fit_deg": "P95 OR residual (°)",
        "fraction_OR_fit_gt5deg_pct": "OR residual >5° (%)",
        "mean_support_score": "Mean heuristic support",
        "low_support_fraction_pct": "Low-support daughters (%)",
        "prior_parent_boundary_edges": "Called parent-boundary edges",
        "prior_parent_boundary_fraction_pct": "Called boundary edges (%)",
        "truth_ARI": "Truth ARI",
        "truth_completeness": "Truth completeness",
        "truth_homogeneity": "Truth homogeneity",
        "truth_V_measure": "Truth V-measure",
        "truth_boundary_precision": "Boundary precision",
        "truth_boundary_recall": "Boundary recall",
        "truth_boundary_F1": "Boundary F1",
        "false_parent_boundary_rate": "False-boundary rate",
        "fragments_per_true_parent_mean": "Fragments / true parent",
        "runtime_s": "Runtime (s)",
    }
    return comp.rename(columns={k:v for k,v in rename.items() if k in comp.columns})

def _agreement_matrix(results: dict[str, ReconstructionResult]) -> pd.DataFrame:
    names = list(results)
    arr = np.eye(len(names), dtype=float)
    for i, a in enumerate(names):
        for j in range(i + 1, len(names)):
            b = names[j]
            score = adjusted_rand_score(results[a].parent_ids, results[b].parent_ids)
            arr[i, j] = arr[j, i] = float(score)
    return pd.DataFrame(arr, index=names, columns=names)


def _dataset_signature(ds: PreparedDataset, r_or: np.ndarray, methods: list[str], controls: dict[str, float]) -> str:
    h = hashlib.sha256()
    h.update(ds.name.encode())
    h.update(pd.util.hash_pandas_object(ds.grains, index=True).values.tobytes())
    h.update(pd.util.hash_pandas_object(ds.adjacency, index=True).values.tobytes())
    h.update(np.asarray(r_or, float).round(10).tobytes())
    h.update("|".join(methods).encode())
    h.update(repr(sorted(controls.items())).encode())
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Streamlit UI
# ---------------------------------------------------------------------------

def _render_or_builder(key_prefix: str = "recon") -> tuple[str, ORPreset, np.ndarray, tuple[np.ndarray, ...], tuple[np.ndarray, ...]]:
    presets = orientation_relationship_presets()
    st.subheader("1 · Define the parent → daughter crystallography")
    st.caption(
        "OR = orientation relationship. Internally the engine stores a daughter→parent crystal-frame rotation matrix R_cp. "
        "You may use a standard literature OR, the NiTi CT/Otsuka–Ren correspondence route, the metric-aware NiTi natural/AQ relation, define parallel crystallographic plane/direction pairs, or supply the matrix directly."
    )
    mode = st.radio(
        "How will the orientation relationship (OR) be defined?",
        [
            "Literature preset — KS / NW / Bain / Pitsch / Burgers",
            "NiTi B2 ↔ B19′ — CT + Otsuka–Ren correspondence (model-derived initial OR)",
            "NiTi B2 ↔ B19′ — natural/AQ OR (metric-aware)",
            "Custom plane + direction parallelism",
            "Custom 3×3 OR matrix",
        ],
        horizontal=False,
        key=f"{key_prefix}_or_mode",
    )
    if mode == "Literature preset — KS / NW / Bain / Pitsch / Burgers":
        name = st.selectbox("Orientation relationship", list(presets), key=f"{key_prefix}_preset")
        p = presets[name]
        r = p.matrix_child_to_parent
        st.info(f"{p.name}: {p.plane_relation}; {p.direction_relation}. {p.note}")
        st.caption("Use these standard ORs only for the transformation class they describe. For NiTi B2↔B19′, use one of the two dedicated NiTi routes below rather than KS/NW/Bain/Pitsch.")
    elif mode == "NiTi B2 ↔ B19′ — CT + Otsuka–Ren correspondence (model-derived initial OR)":
        st.info(
            "NiTi-specific Correspondence Theory (CT) route using the Otsuka–Ren B2→B19′ correspondence matrix together with the measured B2/B19′ metrics. "
            "Important: a lattice correspondence is not itself a unique experimental OR. For orientation reconstruction, the app converts the correspondence deformation to a reproducible **model-derived initial OR** by polar decomposition; use/refine an experimental OR when available."
        )
        with st.expander("NiTi lattice parameters for the CT / Otsuka–Ren model", expanded=True):
            c1, c2, c3, c4, c5 = st.columns(5)
            a_b2 = c1.number_input("a_B2 — Å", min_value=0.1, max_value=20.0, value=3.010, step=0.001, format="%.5f", key=f"{key_prefix}_ct_ab2")
            a_b19 = c2.number_input("a_B19′ — Å", min_value=0.1, max_value=20.0, value=2.898, step=0.001, format="%.5f", key=f"{key_prefix}_ct_am")
            b_b19 = c3.number_input("b_B19′ — Å", min_value=0.1, max_value=20.0, value=4.108, step=0.001, format="%.5f", key=f"{key_prefix}_ct_bm")
            c_b19 = c4.number_input("c_B19′ — Å", min_value=0.1, max_value=20.0, value=4.646, step=0.001, format="%.5f", key=f"{key_prefix}_ct_cm")
            beta = c5.number_input("β_B19′ — °", min_value=60.0, max_value=120.0, value=97.78, step=0.01, format="%.3f", key=f"{key_prefix}_ct_beta")
            st.caption("Input bounds are broad software safety ranges, not physical NiTi limits. For academic work, use lattice parameters measured for the same alloy and temperature as the EBSD data whenever possible.")
        lat = NiTiAQLattice(float(a_b2), float(a_b19), float(b_b19), float(c_b19), float(beta))
        r = niti_ct_otsuka_ren_orientation_relationship(lat)
        diag = niti_ct_otsuka_ren_diagnostics(lat, r)
        name = "NiTi B2→B19′ CT + Otsuka–Ren correspondence"
        p = ORPreset(
            name, "B2 NiTi", "B19′ NiTi", "cubic", "monoclinic", r,
            "Otsuka–Ren B2↔B19′ lattice correspondence + measured metrics",
            "model-derived polar rotation used only as the orientation-reconstruction starting OR",
            "CT/Otsuka–Ren correspondence route. Exact C matrix, lattice metrics and polar bridge are exposed for reproducibility.",
        )
        q1, q2, q3 = st.columns(3)
        q1.metric("det(R_cp)", f"{diag['det_R_cp']:.8f}")
        q2.metric("Polar OR ↔ natural/AQ difference", f"{diag['misorientation_to_natural_AQ_OR_deg']:.3f}°")
        q3.metric("Model volume ratio det(F)", f"{diag['det_F']:.6f}")
        st.caption("det(R_cp)=1 confirms a proper rotation. The CT-polar vs natural/AQ angle is a comparison of two different NiTi starting-OR constructions, not a pass/fail criterion. det(F) is the correspondence-deformation volume ratio.")
        with st.expander("CT / Otsuka–Ren mathematics and exact matrices", expanded=False):
            st.latex(r"C^{A\to M}=\begin{pmatrix}0&1&-1\\0&1&1\\1&0&0\end{pmatrix},\qquad C^{M\to A}=(C^{A\to M})^{-1}")
            st.latex(r"F_{M\leftarrow A}=B_M C^{M\to A}B_A^{-1}=R_{A\to M}U,\qquad R_{cp}=R_{A\to M}^{T}")
            st.markdown(
                """
- **A = B2 parent (austenite)**; **M = B19′ daughter (martensite)**.
- **C** = crystallographic correspondence matrix. It tells which lattice coordinates correspond; it is not by itself an experimental OR.
- **B_A, B_M** = Cartesian direct-lattice basis matrices built from the entered lattice parameters.
- **F** = metric-aware correspondence deformation from the B2 Cartesian basis to the B19′ Cartesian basis.
- **R** = rotational factor of the polar decomposition; **U** = symmetric stretch factor.
- **R_cp** = proper daughter→parent rotation actually consumed by the reconstruction algorithms.
                """
            )
            st.write("C^{A→M} used by the app:")
            st.dataframe(pd.DataFrame(np.asarray(C_A_TO_M, float)), use_container_width=True, hide_index=True)
            st.write("C^{M→A} used by the app:")
            st.dataframe(pd.DataFrame(np.asarray(C_M_TO_A, float)), use_container_width=True, hide_index=True)
            st.write("Correspondence deformation F_{M←A}:")
            st.dataframe(pd.DataFrame(np.asarray(diag['F_child_from_parent'], float)), use_container_width=True, hide_index=True)
            st.write("Principal stretches of F (sorted λ₁≤λ₂≤λ₃):", [float(x) for x in diag["principal_stretches_sorted"]])
            st.caption("Source-derived part: the NiTi B2/B19′ metrics and Otsuka–Ren correspondence matrix. Software bridge: polar decomposition is used here to create an explicit initial OR for an orientation-reconstruction engine. Report this distinction in a paper.")
    elif mode == "NiTi B2 ↔ B19′ — natural/AQ OR (metric-aware)":
        st.info(
            "One-click NiTi setup using the reported natural/AQ relation: (010)B19′ ∥ (110)B2 and [101]B19′ ∥ [−1 1 1]B2. "
            "Unlike a naive Miller-index rotation, the B19′ direct/reciprocal lattice metric is used before constructing R_cp."
        )
        with st.expander("NiTi lattice parameters used to construct this OR", expanded=True):
            c1, c2, c3, c4, c5 = st.columns(5)
            a_b2 = c1.number_input("a_B2 — Å", min_value=0.1, max_value=20.0, value=3.010, step=0.001, format="%.5f", key=f"{key_prefix}_aq_ab2")
            a_b19 = c2.number_input("a_B19′ — Å", min_value=0.1, max_value=20.0, value=2.898, step=0.001, format="%.5f", key=f"{key_prefix}_aq_am")
            b_b19 = c3.number_input("b_B19′ — Å", min_value=0.1, max_value=20.0, value=4.108, step=0.001, format="%.5f", key=f"{key_prefix}_aq_bm")
            c_b19 = c4.number_input("c_B19′ — Å", min_value=0.1, max_value=20.0, value=4.646, step=0.001, format="%.5f", key=f"{key_prefix}_aq_cm")
            beta = c5.number_input("β_B19′ — °", min_value=60.0, max_value=120.0, value=97.78, step=0.01, format="%.3f", key=f"{key_prefix}_aq_beta")
            st.caption(
                "Allowed UI ranges are intentionally broad: lattice lengths 0.1–20 Å and monoclinic β 60–120°. "
                "Use measured values from the same alloy/temperature whenever possible. These numbers affect the metric-aware direction [101]B19′ and therefore the exact OR matrix."
            )
        lat = NiTiAQLattice(float(a_b2), float(a_b19), float(b_b19), float(c_b19), float(beta))
        r = niti_aq_orientation_relationship(lat)
        residuals = niti_aq_parallelism_residuals(lat, r)
        name = "NiTi B2→B19′ natural/AQ OR"
        p = ORPreset(
            name,
            "B2 NiTi",
            "B19′ NiTi",
            "cubic",
            "monoclinic",
            r,
            "(110)B2 ∥ (010)B19′",
            "[−1 1 1]B2 ∥ [101]B19′",
            "Metric-aware natural/AQ relation; exact input lattice parameters are part of the reproducibility record.",
        )
        q1, q2 = st.columns(2)
        q1.metric("Plane-parallelism residual", f"{residuals['plane_parallelism_residual_deg']:.3e}°")
        q2.metric("Direction-parallelism residual", f"{residuals['direction_parallelism_residual_deg']:.3e}°")
        st.caption("These residuals verify the OR construction numerically; 0° is exact within floating-point precision. Literature basis: Crystals 10 (2020) 562, DOI 10.3390/cryst10070562.")
    else:
        c1, c2 = st.columns(2)
        parent_sym_name = c1.selectbox("Parent crystal symmetry", SYMMETRIES, key=f"{key_prefix}_psym")
        child_sym_name = c2.selectbox("Daughter crystal symmetry", SYMMETRIES, key=f"{key_prefix}_csym")
        if mode == "Custom plane + direction parallelism":
            p1, p2 = st.columns(2)
            with p1:
                pp = st.text_input("Parent plane normal [h k l]", "1 1 1", key=f"{key_prefix}_pp")
                pdirection = st.text_input("Parent in-plane direction [u v w]", "1 -1 0", key=f"{key_prefix}_pd")
            with p2:
                cp = st.text_input("Daughter plane normal [h k l]", "1 1 0", key=f"{key_prefix}_cp")
                cdirection = st.text_input("Daughter in-plane direction [u v w]", "1 -1 1", key=f"{key_prefix}_cd")
            r = rotation_from_parallelisms(_parse_vec(pp, "parent plane"), _parse_vec(pdirection, "parent direction"), _parse_vec(cp, "daughter plane"), _parse_vec(cdirection, "daughter direction"))
            plane_rel = f"parent ({pp}) ∥ daughter ({cp})"
            dir_rel = f"parent [{pdirection}] ∥ daughter [{cdirection}]"
            st.warning(
                "Generic custom parallelism treats the supplied vectors in an orthonormal crystal-coordinate frame. For non-orthogonal daughter lattices such as monoclinic B19′, use the dedicated metric-aware NiTi AQ option or supply a validated 3×3 OR matrix."
            )
        else:
            txt = st.text_area(
                "Daughter → parent OR matrix R_cp",
                "1 0 0; 0 1 0; 0 0 1",
                help="Three rows separated by semicolons. The matrix is projected to the nearest proper rotation to remove small numerical non-orthogonality.",
                key=f"{key_prefix}_rmat",
            )
            r = _parse_matrix(txt)
            plane_rel = "custom matrix; no plane parallelism supplied"
            dir_rel = "custom matrix; no direction parallelism supplied"
        name = "Custom OR"
        p = ORPreset(name, "custom parent", "custom daughter", parent_sym_name, child_sym_name, r, plane_rel, dir_rel, "User-defined OR.")
        st.warning("Custom OR: report the exact matrix or the defining parallelisms in any publication; the app cannot infer their physical provenance.")
    parent_sym = symmetry_group(p.parent_symmetry)
    child_sym = symmetry_group(p.child_symmetry)
    nvar = len(unique_child_variants(np.eye(3), r, parent_sym, child_sym))
    a, b, c = st.columns(3)
    a.metric("Parent symmetry", p.parent_symmetry)
    b.metric("Daughter symmetry", p.child_symmetry)
    c.metric("Symmetry-distinct daughter variants", nvar)
    st.caption("Variant count = number of daughter orientations generated by applying parent symmetry and deduplicating under daughter symmetry for the selected OR.")
    with st.expander("Exact OR matrix used by the calculation", expanded=False):
        st.dataframe(pd.DataFrame(np.asarray(r, float), columns=["parent-axis 1", "parent-axis 2", "parent-axis 3"]), use_container_width=True, hide_index=True)
        st.caption(f"R_cp is a proper daughter→parent crystal-frame rotation; det(R_cp) = {np.linalg.det(r):.8f}. Save this matrix with the methods record so the calculation is reproducible.")
    return name, p, r, parent_sym, child_sym

def _manual_grain_editor(key: str) -> pd.DataFrame:
    default = pd.DataFrame([
        {"grain_id": 1, "x": 0.0, "y": 0.0, "phi1_deg": 0.0, "Phi_deg": 0.0, "phi2_deg": 0.0, "area_um2": 1.0, "phase": "daughter"},
        {"grain_id": 2, "x": 1.0, "y": 0.0, "phi1_deg": 0.0, "Phi_deg": 0.0, "phi2_deg": 0.0, "area_um2": 1.0, "phase": "daughter"},
        {"grain_id": 3, "x": 0.5, "y": 1.0, "phi1_deg": 0.0, "Phi_deg": 0.0, "phi2_deg": 0.0, "area_um2": 1.0, "phase": "daughter"},
    ])
    st.caption("Manual grain table: one row = one already-segmented daughter grain. φ₁, Φ, φ₂ are Bunge Euler angles in degrees; x/y are grain-centroid coordinates; area_um2 is grain area in µm² when known.")
    return st.data_editor(default, num_rows="dynamic", use_container_width=True, key=key)


def _manual_adjacency_editor(grain_ids: list[int], key: str) -> pd.DataFrame:
    if len(grain_ids) >= 2:
        default = pd.DataFrame({"grain_id_1": grain_ids[:-1], "grain_id_2": grain_ids[1:]})
    else:
        default = pd.DataFrame(columns=["grain_id_1", "grain_id_2"])
    st.caption("Adjacency table: each row says that two daughter grains physically share a boundary. This is more reliable than centroid k-NN when true boundary topology is available.")
    return st.data_editor(default, num_rows="dynamic", use_container_width=True, key=key)


def _prepare_dataset(
    or_name: str,
    r_child_to_parent: np.ndarray,
    parent_sym: tuple[np.ndarray, ...],
    child_sym: tuple[np.ndarray, ...],
) -> PreparedDataset | None:
    st.subheader("2 · Load and audit the measured daughter microstructure")
    st.caption(
        "The reconstruction engines work on daughter **grains**, not raw pixels. Choose a route that matches what you actually have. "
        "Nothing is silently filtered: phase, quality, segmentation and adjacency choices are shown and become part of the report."
    )
    mode = st.radio(
        "Input route",
        [
            "Built-in validation dataset",
            "Manual grain table",
            "Upload pre-segmented daughter-grain table",
            "Upload raw EBSD map / point table and segment here",
        ],
        horizontal=False,
    )
    convention = "crystal_to_specimen"
    length_unit = "µm"

    if mode == "Built-in validation dataset":
        c1, c2 = st.columns(2)
        n = c1.number_input(
            "Daughter grains per known parent — count",
            min_value=4, max_value=30, value=6, step=1,
            help="Synthetic validation only. Larger values give more daughter evidence per known parent."
        )
        noise = c2.number_input(
            "Synthetic orientation noise — degrees",
            min_value=0.0, max_value=3.0, value=0.35, step=0.05,
            help="Random angular perturbation applied to the synthetic daughter orientations."
        )
        df, edges, _ = synthetic_parent_reconstruction_demo_from_or(
            r_child_to_parent, parent_sym, child_sym, int(n), float(noise), seed=11
        )
        st.success(
            "Known truth: exactly two prior-parent grains. The synthetic daughters are generated with the SAME selected OR and parent/daughter symmetries shown in Step 1, so this is a transformation-matched pipeline validation."
        )
        st.caption(f"Synthetic crystallography used to generate truth: {or_name}. This prevents an FCC→BCC validation dataset from being accidentally judged with a NiTi B2→B19′ model (or vice versa).")
        return PreparedDataset(
            "built-in validation", df, edges, convention,
            f"synthetic two-parent validation generated from selected OR: {or_name}",
            "known/shared-boundary synthetic graph", "arbitrary",
            "true_parent_id is validation truth; true_parent_boundary marks the deliberately inserted inter-parent edge"
        )

    if mode == "Manual grain table":
        df = _manual_grain_editor("manual_grains")
        provenance = "manual pre-segmented daughter-grain table"
        length_unit = st.selectbox("x/y coordinate unit", ["µm", "nm", "mm", "arbitrary/map unit"], index=0, key="manual_length_unit")
    else:
        allowed = ["csv", "tsv", "txt"] if mode == "Upload pre-segmented daughter-grain table" else ["ang", "ctf", "csv", "tsv", "txt"]
        uploaded = st.file_uploader(
            "Daughter orientation file",
            type=allowed,
            accept_multiple_files=False,
            help=(
                "Pre-segmented route: one row per daughter grain with grain_id + Euler/quaternion orientation. "
                "Raw route: .ang, .ctf, or a delimited point table with x/y + Euler angles."
            ),
        )
        if uploaded is None:
            if mode == "Upload pre-segmented daughter-grain table":
                template = pd.DataFrame([
                    {"grain_id": 1, "x": 0.0, "y": 0.0, "phi1_deg": 10.0, "Phi_deg": 20.0, "phi2_deg": 30.0, "area_um2": 1.0, "phase": "daughter"}
                ])
                st.download_button("Download daughter-grain CSV template", template.to_csv(index=False), "daughter_grains_template.csv", "text/csv")
            st.info("Upload a file to continue.")
            return None
        try:
            df, provenance = _read_uploaded_orientation_file(uploaded)
        except Exception as exc:
            st.error(f"Import rejected: {exc}")
            st.caption("The importer fails closed when required orientation/position fields cannot be identified; it does not guess missing crystallography.")
            return None

        raw_route = mode == "Upload raw EBSD map / point table and segment here"
        if raw_route:
            if "grain_id" in df.columns:
                st.info("A grain_id column exists, but you selected the raw-point route. The app will segment from point positions/orientations and will not trust those grain IDs automatically.")
            if not {"x", "y", "phi1_deg", "Phi_deg", "phi2_deg"}.issubset(df.columns):
                st.error("Raw EBSD segmentation requires x, y, phi1_deg, Phi_deg and phi2_deg after import.")
                return None

            st.markdown("**Import audit — check this before segmentation**")
            st.dataframe(_raw_ebsd_audit(df), use_container_width=True, hide_index=True)
            if df.attrs.get("angle_note"):
                st.caption(str(df.attrs["angle_note"]))

            # Explicit daughter-phase selection.  Phase 0 is conventionally not indexed.
            phase_numeric = pd.to_numeric(df.get("phase", pd.Series(np.ones(len(df)))), errors="coerce")
            phases = sorted(float(x) for x in phase_numeric.dropna().unique() if float(x) > 0)
            if not phases:
                st.error("No indexed phase (>0) was found. Check the phase column or export settings.")
                return None
            phase_counts = pd.DataFrame({
                "phase_id": phases,
                "indexed_points": [int(np.sum(phase_numeric == x)) for x in phases],
            })
            phase_counts["fraction_of_indexed_pct"] = 100.0 * phase_counts["indexed_points"] / max(1, int(np.sum(phase_numeric > 0)))
            st.dataframe(phase_counts, use_container_width=True, hide_index=True)
            selected_phase = st.selectbox(
                "Daughter phase ID to reconstruct",
                phases,
                index=0,
                help="Only this indexed phase is segmented as the daughter phase. Phase 0/non-indexed points and other phases are excluded."
            )
            work = df[phase_numeric == selected_phase].copy().reset_index(drop=True)

            with st.expander("Optional EBSD quality filters — OFF by default", expanded=False):
                st.caption("Quality filters are acquisition-dependent. They are never applied silently. If enabled, report the field and threshold.")
                filter_notes = []
                if "CI" in work.columns:
                    use_ci = st.checkbox("Filter by minimum EDAX confidence index (CI)", value=False)
                    ci_min = st.number_input("Minimum CI", min_value=-1.0, max_value=1.0, value=0.05, step=0.01, disabled=not use_ci)
                    if use_ci:
                        before = len(work)
                        work = work[pd.to_numeric(work["CI"], errors="coerce") >= float(ci_min)].copy()
                        filter_notes.append(f"CI >= {ci_min:g} ({before}->{len(work)} points)")
                if "MAD" in work.columns:
                    use_mad = st.checkbox("Filter by maximum Oxford mean angular deviation (MAD)", value=False)
                    mad_max = st.number_input("Maximum MAD", min_value=0.0, max_value=10.0, value=1.5, step=0.1, disabled=not use_mad)
                    if use_mad:
                        before = len(work)
                        work = work[pd.to_numeric(work["MAD"], errors="coerce") <= float(mad_max)].copy()
                        filter_notes.append(f"MAD <= {mad_max:g} ({before}->{len(work)} points)")
                if "fit" in work.columns:
                    use_fit = st.checkbox("Filter by maximum ANG fit field", value=False)
                    fit_max = st.number_input("Maximum ANG fit", min_value=0.0, max_value=10.0, value=2.0, step=0.1, disabled=not use_fit)
                    if use_fit:
                        before = len(work)
                        work = work[pd.to_numeric(work["fit"], errors="coerce") <= float(fit_max)].copy()
                        filter_notes.append(f"fit <= {fit_max:g} ({before}->{len(work)} points)")
                if filter_notes:
                    st.write("Applied filters:", "; ".join(filter_notes))
                else:
                    st.write("Applied filters: none")

            if len(work) < 2:
                st.error("Fewer than two daughter-phase points remain after selection/filtering.")
                return None

            cunit = st.selectbox(
                "Map coordinate unit",
                ["µm", "nm", "mm", "arbitrary/map unit"],
                index=0,
                help="ANG/CTF coordinates carry a scale but the lightweight parser cannot prove the physical unit from every vendor export. Confirm it from acquisition metadata."
            )
            length_unit = cunit
            st.warning(
                "Reference-frame check required: vendor EBSD formats can encode map/Euler reference-frame conventions differently. "
                "Verify the specimen axes and Euler reference frame against the acquisition software or a trusted EBSD package before publishing orientation-dependent results."
            )
            frame_ok = st.checkbox("I verified the Euler-angle and specimen/map reference-frame convention for this export", value=False)
            if not frame_ok:
                st.info("Confirm the reference frame above to enable segmentation. This prevents a plausible-looking reconstruction from an unverified orientation convention.")
                return None

            st.markdown("**Daughter-grain segmentation**")
            s1, s2, s3 = st.columns(3)
            seg_deg = s1.number_input(
                "Grain-boundary disorientation threshold — degrees",
                min_value=0.5, max_value=15.0, value=3.0, step=0.25,
                help="Common EBSD workflows often start near a few degrees, but this is material/data dependent. Report the exact value."
            )
            min_pts = s2.number_input(
                "Minimum indexed points per retained daughter grain",
                min_value=2, max_value=1000, value=5, step=1,
                help="Small objects below this point count are discarded from reconstruction. Report the value."
            )
            radius_factor = s3.number_input(
                "Spatial neighbor radius / median point spacing",
                min_value=1.05, max_value=1.70, value=1.35, step=0.05,
                help="1.35 captures nearest neighbors on common square/hex grids without normally bridging square-grid diagonals. Change only when the map geometry requires it."
            )
            st.caption(f"Selected daughter points ready for segmentation: {len(work):,}. Coordinate unit: {length_unit}.")

            raw_hash = hashlib.sha256(uploaded.getvalue()).hexdigest()
            seg_key = f"{raw_hash}|phase={selected_phase}|n={len(work)}|thr={seg_deg}|min={min_pts}|rf={radius_factor}"
            if st.button("Segment daughter grains", type="primary", use_container_width=True):
                try:
                    grains, adj = _segment_points_to_grains(
                        work, child_sym, float(seg_deg), int(min_pts), neighbor_radius_factor=float(radius_factor)
                    )
                    st.session_state["seg_grains"] = grains
                    st.session_state["seg_adj"] = adj
                    st.session_state["seg_key"] = seg_key
                    st.session_state["seg_notes"] = {
                        "phase": selected_phase,
                        "threshold_deg": float(seg_deg),
                        "min_points": int(min_pts),
                        "neighbor_radius_factor": float(radius_factor),
                        "length_unit": length_unit,
                    }
                except Exception as exc:
                    st.error(f"Segmentation failed: {exc}")
            if st.session_state.get("seg_key") != seg_key:
                st.info("Run segmentation with the current phase/filter/threshold settings.")
                return None

            grains = st.session_state["seg_grains"].copy()
            adj = st.session_state["seg_adj"].copy()
            # Relabel generic map-unit areas if user confirmed µm.
            if "area_est_mapunit2" in grains.columns and length_unit == "µm":
                grains["area_um2"] = grains["area_est_mapunit2"]
            a1, a2, a3, a4 = st.columns(4)
            a1.metric("Retained daughter grains", len(grains))
            a2.metric("Shared-boundary grain pairs", len(adj))
            a3.metric("Median points / grain", f"{float(grains['point_count'].median()):.1f}")
            a4.metric("Median grain orientation spread", f"{float(grains['grain_orientation_spread_rms_deg'].median()):.3f}°")
            st.dataframe(grains.head(40), use_container_width=True, hide_index=True)
            if len(adj) == 0:
                st.error("Segmentation produced no inter-grain adjacency. Parent reconstruction methods need a grain-neighbor graph.")
                return None
            return PreparedDataset(
                uploaded.name, grains, adj, convention,
                f"{provenance}; explicit phase selection + in-app spatial/disorientation segmentation",
                "segmented shared-boundary graph", length_unit,
                json.dumps(st.session_state.get("seg_notes", {}), sort_keys=True),
            )

        # Pre-segmented uploaded grain table continues below.
        if "grain_id" not in df.columns:
            st.error("The pre-segmented route requires one row per daughter grain and a unique grain_id column. Choose the raw EBSD route for pixel/point data.")
            return None
        length_unit = st.selectbox("x/y coordinate unit", ["µm", "nm", "mm", "arbitrary/map unit"], index=0, key="upload_length_unit")

    # Manual or pre-segmented daughter-grain data.
    df = _coerce_grain_columns(df)
    st.markdown("**Orientation representation and convention**")
    has_euler = {"phi1_deg", "Phi_deg", "phi2_deg"}.issubset(df.columns)
    has_quat = {"qw", "qx", "qy", "qz"}.issubset(df.columns)
    st.write("Detected orientation columns:", "Bunge Euler angles" if has_euler else "quaternion (w,x,y,z)" if has_quat else "not recognized")
    conv = st.radio(
        "Orientation mapping convention",
        ["Crystal → specimen", "Specimen → crystal"],
        horizontal=True,
        help="Internal reconstruction uses crystal→specimen proper rotations. Choose specimen→crystal only when your export explicitly uses that inverse convention."
    )
    convention = "crystal_to_specimen" if conv == "Crystal → specimen" else "specimen_to_crystal"
    try:
        grain_ids, _ = orientations_from_dataframe(df, convention=convention)
    except Exception as exc:
        st.error(f"Daughter-grain orientation table rejected: {exc}")
        st.caption("Required: unique grain_id plus either phi1_deg/Phi_deg/phi2_deg in degrees, or qw/qx/qy/qz. Optional: x, y, area_um2, phase, true_parent_id.")
        return None

    if has_euler:
        e = df[["phi1_deg", "Phi_deg", "phi2_deg"]].apply(pd.to_numeric, errors="coerce")
        if e.isna().any().any():
            st.error("Euler-angle columns contain missing or non-numeric values.")
            return None
        if ((e["Phi_deg"] < -1e-6) | (e["Phi_deg"] > 180.0 + 1e-6)).any():
            st.warning("Some Bunge Φ values are outside the common 0–180° range. This can be valid under alternate parameterizations, but verify the export convention.")

    st.markdown("**Physical daughter-grain adjacency**")
    adj_mode_options = ["Upload measured/shared-boundary adjacency", "Enter adjacency manually"]
    if {"x", "y"}.issubset(df.columns):
        adj_mode_options.append("Approximate from centroids (k-nearest neighbors)")
    adj_mode = st.radio("How are daughter-grain neighbors defined?", adj_mode_options, horizontal=False)
    if adj_mode.startswith("Upload"):
        af = st.file_uploader(
            "Adjacency CSV/TSV",
            type=["csv", "tsv", "txt"],
            help="Required columns: grain_id_1, grain_id_2. Optional extra columns such as boundary length are preserved."
        )
        if af is None:
            template = pd.DataFrame([{"grain_id_1": int(grain_ids[0]), "grain_id_2": int(grain_ids[min(1, len(grain_ids)-1)])}])
            st.download_button("Download adjacency template", template.to_csv(index=False), "daughter_adjacency_template.csv", "text/csv")
            st.info("Upload adjacency to continue, or choose another adjacency route.")
            return None
        adj = pd.read_csv(af, sep=None, engine="python")
        graph_kind = "measured/shared-boundary"
    elif adj_mode.startswith("Enter"):
        adj = _manual_adjacency_editor([int(x) for x in grain_ids], "manual_adj")
        graph_kind = "manual/shared-boundary claim"
    else:
        k = st.number_input(
            "Centroid neighbors per daughter grain — k",
            min_value=1, max_value=min(20, max(1, len(df)-1)), value=min(4, max(1, len(df)-1)), step=1,
            help="Approximation only. Larger k adds non-boundary neighbors and can bias graph reconstruction."
        )
        st.warning("Centroid k-NN is not measured grain-boundary topology. Use it only when true adjacency is unavailable and report that limitation explicitly.")
        idx_edges = approximate_knn_edges(df, int(k))
        adj = pd.DataFrame({
            "grain_id_1": [int(grain_ids[i]) for i, _ in idx_edges],
            "grain_id_2": [int(grain_ids[j]) for _, j in idx_edges],
        })
        graph_kind = f"approximate centroid k-NN (k={int(k)})"

    try:
        idx_edges = edges_from_dataframe(adj, grain_ids)
    except Exception as exc:
        st.error(f"Adjacency table rejected: {exc}")
        return None

    audit = audit_reconstruction_input(df, idx_edges, graph_kind)
    a1, a2, a3, a4 = st.columns(4)
    a1.metric("Daughter grains", audit.n_grains)
    a2.metric("Neighbor pairs", audit.n_edges)
    a3.metric("Grains with ≥1 neighbor", f"{100*audit.connected_fraction:.1f}%")
    a4.metric("Isolated grains", audit.isolated_grains)
    for warning in audit.warnings:
        st.warning(warning)
    with st.expander("Method readiness for this input", expanded=False):
        st.dataframe(recommend_methods(audit), use_container_width=True, hide_index=True)
    return PreparedDataset(
        "manual" if mode == "Manual grain table" else uploaded.name,
        df, adj, convention, provenance, graph_kind, length_unit,
        "Pre-segmented grain orientations supplied by user; no grain segmentation performed in this route.",
    )


def _render_method_controls(selected: list[str]) -> dict[str, float]:
    st.subheader("4 · Method-specific controls")
    st.caption(
        "Only controls actually used by the selected algorithms are shown. Angular values are degrees. "
        "The ranges are safety/UI ranges, not material constants; the exact chosen values belong in the Methods section."
    )
    with st.expander("What each reconstruction method needs and returns", expanded=False):
        st.dataframe(method_requirements_table(), use_container_width=True, hide_index=True)

    controls: dict[str, float] = {}
    if any(m in selected for m in ["Neighbor voting", "Grain graph + Markov clustering", "Variant graph probability propagation"]):
        controls["sigma_deg"] = float(st.number_input(
            "Local candidate/edge angular weighting width σ — degrees",
            min_value=0.2, max_value=10.0, value=2.5, step=0.1,
            help="Smaller σ penalizes angular disagreement more strongly; larger σ makes the method more permissive. This is an algorithmic width, not a measured material property."
        ))
    else:
        controls["sigma_deg"] = 2.5
    if any(m in selected for m in ["Grain graph + Markov clustering", "Variant graph probability propagation"]):
        controls["inflation"] = float(st.number_input(
            "Markov clustering inflation — dimensionless",
            min_value=1.05, max_value=3.0, value=1.6, step=0.05,
            help="Higher inflation generally produces more strongly separated/smaller graph clusters. Report it because it can change prior-parent partitioning."
        ))
    else:
        controls["inflation"] = 1.6
    if "Neighbor voting" in selected:
        controls["merge_deg"] = float(st.number_input(
            "Neighbor-voting parent merge tolerance — degrees",
            min_value=0.5, max_value=15.0, value=5.0, step=0.25,
            help="Selected local parent candidates are merged across daughter adjacency only when their parent disorientation is within this value."
        ))
    else:
        controls["merge_deg"] = 5.0
    if "Nucleation + growth" in selected:
        a, b = st.columns(2)
        controls["nucleation_deg"] = float(a.number_input(
            "Nucleation tolerance — degrees",
            min_value=0.5, max_value=10.0, value=3.0, step=0.25,
            help="Strict criterion for forming parent seeds. Lower = fewer but cleaner seeds."
        ))
        controls["growth_deg"] = float(b.number_input(
            "Growth tolerance — degrees",
            min_value=1.0, max_value=20.0, value=10.0, step=0.5,
            help="Looser criterion for adding neighboring daughter grains after a seed exists. Must be interpreted together with the nucleation tolerance."
        ))
        if controls["growth_deg"] < controls["nucleation_deg"]:
            st.warning("Growth tolerance is normally at least as permissive as nucleation tolerance. Your current growth value is smaller.")
    else:
        controls["nucleation_deg"] = 3.0
        controls["growth_deg"] = 10.0
    if "Operator / groupoid consistency" in selected:
        a, b = st.columns(2)
        controls["operator_tol_deg"] = float(a.number_input(
            "Daughter–daughter operator residual tolerance — degrees",
            min_value=0.5, max_value=15.0, value=5.0, step=0.25,
            help="A measured daughter-neighbor misorientation must lie this close to one theoretical transformation operator before the edge is accepted for operator-consistent growth."
        ))
        controls["parent_consistency_deg"] = float(b.number_input(
            "Common-parent consistency tolerance — degrees",
            min_value=0.5, max_value=15.0, value=5.0, step=0.25,
            help="After an operator match, candidate parent orientations must agree within this angular tolerance."
        ))
        st.caption("This route follows the operator/groupoid reconstruction idea and now exports operator IDs/residual statistics, but it is a transparent research implementation rather than a binary clone of ARPGE/GenOVa.")
    else:
        controls["operator_tol_deg"] = 5.0
        controls["parent_consistency_deg"] = 5.0
    return controls


def _run_one(
    name: str,
    grain_ids: np.ndarray,
    child_orientations: list[np.ndarray],
    edges: list[tuple[int, int]],
    r_or: np.ndarray,
    child_sym: tuple[np.ndarray, ...],
    parent_sym: tuple[np.ndarray, ...],
    c: dict[str, float],
) -> ReconstructionResult:
    if name == "Neighbor voting":
        return neighbor_voting_reconstruction(grain_ids, child_orientations, edges, r_or, child_sym, parent_sym, sigma_deg=c["sigma_deg"], merge_deg=c["merge_deg"])
    if name == "Grain graph + Markov clustering":
        return grain_graph_reconstruction(grain_ids, child_orientations, edges, r_or, child_sym, parent_sym, sigma_deg=c["sigma_deg"], inflation=c["inflation"])
    if name == "Variant graph probability propagation":
        return variant_graph_reconstruction(grain_ids, child_orientations, edges, r_or, child_sym, parent_sym, sigma_deg=max(c["sigma_deg"], 0.5), inflation=max(c["inflation"], 1.0))
    if name == "Nucleation + growth":
        return nucleation_growth_reconstruction(grain_ids, child_orientations, edges, r_or, child_sym, parent_sym, nucleation_deg=c["nucleation_deg"], growth_deg=c["growth_deg"])
    return operator_groupoid_reconstruction(grain_ids, child_orientations, edges, r_or, child_sym, parent_sym, operator_tol_deg=c["operator_tol_deg"], parent_consistency_deg=c["parent_consistency_deg"])



def _parent_summary_display(df: pd.DataFrame) -> pd.DataFrame:
    rename = {
        "parent_id": "Parent ID",
        "supporting_daughter_grains": "Supporting daughters",
        "distinct_selected_variants": "Distinct selected variants",
        "dominant_variant_id": "Dominant variant ID",
        "dominant_variant_fraction_pct": "Dominant variant (%)",
        "supporting_area": "Supporting area",
        "area_fraction_pct": "Area fraction (%)",
        "mean_OR_fit_deg": "Mean OR residual (°)",
        "median_OR_fit_deg": "Median OR residual (°)",
        "p95_OR_fit_deg": "P95 OR residual (°)",
        "max_OR_fit_deg": "Max OR residual (°)",
        "mean_candidate_separation_deg": "Mean candidate separation (°)",
        "mean_support_score": "Mean heuristic support",
        "minimum_support_score": "Minimum heuristic support",
        "ambiguous_daughter_grains": "Ambiguous daughters",
        "parent_phi1_deg": "Parent φ₁ (°)",
        "parent_Phi_deg": "Parent Φ (°)",
        "parent_phi2_deg": "Parent φ₂ (°)",
    }
    out = df.rename(columns={k:v for k,v in rename.items() if k in df.columns}).copy()
    return out


def _daughter_assignment_display(df: pd.DataFrame) -> pd.DataFrame:
    rename = {
        "grain_id": "Daughter grain ID",
        "parent_id": "Reconstructed parent ID",
        "candidate_variant_id": "Candidate variant ID",
        "best_OR_fit_deg": "Best OR residual (°)",
        "second_best_candidate_fit_deg": "2nd-best candidate residual (°)",
        "candidate_separation_deg": "Candidate separation (°)",
        "candidate_count": "Candidate count",
        "absolute_fit_support_0_to_1": "Absolute-fit support (0–1)",
        "candidate_separation_support_0_to_1": "Separation support (0–1)",
        "support_score_0_to_1": "Heuristic support (0–1)",
        "quality_flag": "Review flag",
        "x": "x", "y": "y",
    }
    return df.rename(columns={k:v for k,v in rename.items() if k in df.columns}).copy()

def _render_result(
    name: str,
    result: ReconstructionResult,
    source_df: pd.DataFrame,
    comparison: pd.DataFrame,
    *,
    child_orientations: list[np.ndarray] | None = None,
    edges: list[tuple[int, int]] | None = None,
    grain_ids: np.ndarray | None = None,
    r_or: np.ndarray | None = None,
    parent_sym: tuple[np.ndarray, ...] | None = None,
    child_sym: tuple[np.ndarray, ...] | None = None,
    controls: dict[str, float] | None = None,
) -> None:
    st.subheader(f"Detailed evidence · {name}")
    parent_summary = _parent_summary(result, source_df)
    assignments = _assignment_table(result, source_df)
    q = reconstruction_quality_summary(result)

    m1, m2, m3, m4, m5, m6 = st.columns(6)
    m1.metric("Prior parents", int(result.diagnostics["n_reconstructed_parents"]))
    m2.metric("Singleton parents", f"{100*float(q['singleton_parent_fraction']):.1f}%")
    m3.metric("Mean OR residual", f"{float(np.mean(result.fit_deg)):.3f}°")
    m4.metric("P95 OR residual", f"{float(np.quantile(result.fit_deg, 0.95)):.3f}°")
    m5.metric("Fit ≤2.5°", f"{100*float(q['fraction_fit_le_2_5deg']):.1f}%")
    m6.metric("Mean support", f"{float(np.mean(result.confidence)):.3f}")
    if float(q["singleton_parent_fraction"]) > 0.50:
        st.error("Structural warning: more than half of reconstructed parents contain only one daughter grain. Near-zero OR residuals can then be trivial self-fits and must not be interpreted as a successful parent reconstruction.")
    st.caption(
        "OR residual = angular disagreement between the reconstructed parent orientation and the best crystallographically allowed parent candidate for a daughter grain. "
        "Support = a 0–1 heuristic combining absolute fit and separation from the next candidate; it is not a posterior probability. "
        f"Overall interpretation: {q['interpretation']}."
    )

    t1, t2, t3, t4, t5 = st.tabs([
        "Parent summary", "Daughter evidence", "Variants / operators", "Maps & distributions", "How to interpret"
    ])
    with t1:
        st.markdown("**One row = one reconstructed prior-parent grain.** This is normally the table to cite/report first.")
        st.dataframe(_parent_summary_display(parent_summary), use_container_width=True, hide_index=True)
        st.caption(
            "Key fields: supporting_daughter_grains = independent daughter-grain evidence count; distinct_selected_variants = selected crystallographic candidate diversity; "
            "mean/median/P95/max OR fit = residual distribution in degrees; candidate separation = ambiguity against the next candidate; quality_flag is descriptive, not a theorem."
        )
        with st.expander("Column meanings for the parent table"):
            dictionary = {
                "parent_id": "Arbitrary reconstructed cluster label; the number itself has no crystallographic meaning.",
                "supporting_daughter_grains": "Number of measured daughter grains assigned to this reconstructed parent.",
                "distinct_selected_variants": "How many different selected candidate-variant IDs occur among those daughter grains.",
                "dominant_variant_fraction_pct": "Fraction of the parent's daughter grains using the most frequent selected candidate variant.",
                "supporting_area / area_fraction_pct": "Summed daughter area and its share of the analyzed daughter area when area metadata exist.",
                "mean_OR_fit_deg / median_OR_fit_deg / p95_OR_fit_deg / max_OR_fit_deg": "Angular reconstruction residual statistics; smaller indicates closer agreement with the chosen OR.",
                "mean_candidate_separation_deg": "Average gap between best and second-best crystallographic parent candidate; larger means less candidate ambiguity.",
                "mean_support_score / minimum_support_score": "Heuristic 0–1 support, not a calibrated probability.",
                "parent_phi1/Phi/phi2 and quaternion": "The reconstructed parent orientation in Bunge Euler and unit-quaternion forms.",
            }
            st.dataframe(pd.DataFrame([{"field": k, "meaning": v} for k, v in dictionary.items()]), use_container_width=True, hide_index=True)
        st.download_button("Download parent summary CSV", parent_summary.to_csv(index=False), f"{name.replace(' ','_')}_parent_summary.csv", "text/csv")

    with t2:
        st.markdown("**One row = one measured daughter grain and the evidence behind its parent assignment.**")
        qfilter = st.multiselect(
            "Quality bands to display",
            ["Very strong", "Strong", "Review", "Weak / ambiguous"],
            default=["Very strong", "Strong", "Review", "Weak / ambiguous"],
            key=f"quality_filter_{name}",
        )
        view = assignments[assignments["quality_flag"].isin(qfilter)]
        st.dataframe(_daughter_assignment_display(view), use_container_width=True, hide_index=True)
        st.caption(
            "best_OR_fit_deg = best candidate-to-parent angular residual. second_best_candidate_fit_deg and candidate_separation_deg expose ambiguity. "
            "candidate_variant_id is an internal symmetry-generated candidate index unless a transformation-specific canonical packet/Bain/variant mapping is supplied."
        )
        with st.expander("Exact heuristic support formula"):
            st.latex(r"s_{abs}=\exp[-\tfrac12(\theta_{best}/3^\circ)^2]")
            st.latex(r"s_{sep}=1-\exp[-(\theta_{2nd}-\theta_{best})/3^\circ]")
            st.latex(r"s=0.55\,s_{abs}+0.45\,s_{sep}")
            st.caption("θ_best = best candidate fit; θ_2nd = second-best candidate fit. This deterministic score is for ranking/flagging ambiguity only; it has not been calibrated as a statistical probability.")
        st.download_button("Download daughter evidence CSV", assignments.to_csv(index=False), f"{name.replace(' ','_')}_daughter_assignments.csv", "text/csv")

    with t3:
        variants = variant_frequency_table(result, source_df)
        st.markdown("**Selected candidate-variant frequency**")
        st.dataframe(variants, use_container_width=True, hide_index=True)
        st.caption("Variant-frequency tables are useful for variant-selection studies. Internal candidate IDs are reproducible for this calculation but should not be called canonical KS/Burgers packet/Bain labels unless that mapping is explicitly defined.")
        st.download_button("Download variant-frequency CSV", variants.to_csv(index=False), f"{name.replace(' ','_')}_variant_frequency.csv", "text/csv")

        if name == "Operator / groupoid consistency" and all(x is not None for x in [child_orientations, edges, grain_ids, r_or, parent_sym, child_sym, controls]):
            optol = float((controls or {}).get("operator_tol_deg", result.diagnostics.get("operator_tol_deg", 5.0)))
            op_edges = operator_edge_table(
                child_orientations or [], edges or [], grain_ids if grain_ids is not None else [],
                np.asarray(r_or), parent_sym or tuple(), child_sym or tuple(), optol, result=result,
            )
            op_freq = operator_frequency_table(op_edges, accepted_only=True)
            st.markdown("**Operator/groupoid edge classification**")
            st.dataframe(op_edges, use_container_width=True, hide_index=True)
            st.caption("For each neighboring daughter-grain pair: nearest_theoretical_operator_id = closest OR-induced daughter–daughter operator; operator_residual_deg = angular distance to that operator modulo daughter symmetry; within_operator_tolerance = whether the edge passes the user threshold.")
            st.markdown("**Accepted operator frequency**")
            st.dataframe(op_freq, use_container_width=True, hide_index=True)
            st.caption("This gives ARPGE-style operator statistics while keeping the residual and acceptance threshold visible. Operator IDs here follow this app's deterministic generated-operator ordering; do not equate them to an external package's canonical O1/O2/... numbering unless an explicit mapping is supplied.")
            c1, c2 = st.columns(2)
            c1.download_button("Download operator-edge CSV", op_edges.to_csv(index=False), "operator_groupoid_edges.csv", "text/csv", use_container_width=True)
            c2.download_button("Download operator-frequency CSV", op_freq.to_csv(index=False), "operator_groupoid_frequency.csv", "text/csv", use_container_width=True)

    with t4:
        if {"x", "y"}.issubset(assignments.columns):
            fig = px.scatter(
                assignments,
                x="x", y="y", color="parent_id", symbol="candidate_variant_id",
                hover_data=[c for c in ["grain_id", "best_OR_fit_deg", "candidate_separation_deg", "support_score_0_to_1", "quality_flag"] if c in assignments.columns],
                title="Reconstructed prior-parent map from daughter-grain centroids",
                labels={"x": "Map x coordinate", "y": "Map y coordinate", "parent_id": "Prior-parent ID", "candidate_variant_id": "Candidate variant ID"},
            )
            fig.update_yaxes(scaleanchor="x", scaleratio=1)
            st.plotly_chart(fig, use_container_width=True)
        c1, c2 = st.columns(2)
        with c1:
            st.plotly_chart(px.histogram(assignments, x="best_OR_fit_deg", nbins=30, title="OR residual distribution", labels={"best_OR_fit_deg": "Best OR residual (degrees)"}), use_container_width=True)
        with c2:
            st.plotly_chart(px.histogram(assignments, x="candidate_separation_deg", nbins=30, title="Candidate-separation distribution", labels={"candidate_separation_deg": "Second-best minus best candidate fit (degrees)"}), use_container_width=True)
        c3, c4 = st.columns(2)
        with c3:
            st.plotly_chart(px.histogram(assignments, x="support_score_0_to_1", nbins=20, title="Heuristic support distribution", labels={"support_score_0_to_1": "Support score (0–1, not probability)"}), use_container_width=True)
        with c4:
            if "grain_orientation_spread_rms_deg" in assignments.columns:
                st.plotly_chart(px.histogram(assignments, x="grain_orientation_spread_rms_deg", nbins=25, title="Daughter grain orientation spread", labels={"grain_orientation_spread_rms_deg": "RMS intragranular spread (degrees)"}), use_container_width=True)

    with t5:
        st.markdown(
            "**What you can claim from this result:** the algorithm produced a specific parent partition and parent orientations that are internally consistent with the supplied daughter orientations, adjacency, symmetries, OR and thresholds to the residuals shown above."
        )
        st.markdown(
            "**What you should not claim from this table alone:** that the reconstruction is uniquely true, that the heuristic support is a probability, or that a low OR residual proves a correct prior-parent boundary. Use cross-method boundary/cluster agreement, retained-parent validation, known truth, morphology and sensitivity to thresholds when available."
        )
        st.write("Automatic internal consistency summary:", q["interpretation"])
        st.dataframe(comparison[comparison["method"] == name], use_container_width=True, hide_index=True)

    st.download_button("Download method-comparison CSV", comparison.to_csv(index=False), "reconstruction_method_comparison.csv", "text/csv")


def _reference_parent_validation(result: ReconstructionResult, ref_df: pd.DataFrame, parent_sym: tuple[np.ndarray, ...], convention: str) -> pd.DataFrame:
    ids, refs = orientations_from_dataframe(ref_df, convention=convention)
    rows = []
    for pid, gp in result.parent_orientations.items():
        vals = [(misorientation_deg(gp, gref, parent_sym), int(rid)) for rid, gref in zip(ids, refs)]
        fit, rid = min(vals)
        rows.append({"reconstructed_parent_id": int(pid), "nearest_reference_parent_id": rid, "misorientation_to_reference_deg": float(fit)})
    return pd.DataFrame(rows)


def render_reconstruction_workbench() -> None:
    st.header("Parent ↔ daughter orientation reconstruction")
    st.markdown(
        "**Academic workflow, without hiding the crystallography:** start with the measured daughter map, define the transformation crystallography, reconstruct the prior parent, then inspect residuals and cross-method agreement before making a scientific interpretation."
    )
    st.info("Minimum reconstruction input = daughter orientations + parent/daughter symmetries + an orientation relationship (OR) + daughter-grain adjacency. For the dedicated NiTi CT/Otsuka–Ren route, measured B2/B19′ lattice parameters are additionally used to build the metric/correspondence model-derived starting OR.")

    with st.expander("Start here · which workspace do I need?", expanded=True):
        st.dataframe(pd.DataFrame([
            {"goal": "I measured daughter EBSD and want the prior parent", "use": "1 · Reconstruct parent", "main academic output": "parent map + parent orientations + daughter variant assignments + OR residuals + method agreement"},
            {"goal": "I have B19′, reconstructed B2, and want a B19′ re-transformation check", "use": "2 · NiTi cycle", "main academic output": "B19′→B2→B19′ closure + regenerated branches + later-cycle comparison"},
            {"goal": "I know parent orientation(s) and want allowed daughter orientations", "use": "3 · Forward & batch", "main academic output": "symmetry-distinct daughter orientation library"},
            {"goal": "I need definitions, assumptions, equations or literature mapping", "use": "4 · Academic guide", "main academic output": "notation + methods scope + reporting checklist + references"},
        ]), use_container_width=True, hide_index=True)
        st.caption("Recommended paper workflow: reconstruct parent → inspect the parent/daughter tables → compare methods → validate with retained parent or independent data when available. The app does not convert internal agreement into a universal 'correct/incorrect' verdict.")

    tab_recon, tab_cycle, tab_tools, tab_academic = st.tabs([
        "1 · Reconstruct parent", "2 · NiTi B19′↔B2 cycle", "3 · Forward & batch", "4 · Academic guide"
    ])
    # Keep every existing workflow; two pairs are simply grouped into cleaner top-level tabs.
    tab_forward = tab_tools
    tab_batch = tab_tools
    tab_dictionary = tab_academic
    tab_sources = tab_academic

    with tab_recon:
        # OR is needed before raw EBSD segmentation because segmentation uses daughter symmetry.
        or_name, preset, initial_or, parent_sym, child_sym = _render_or_builder("main")
        ds = _prepare_dataset(or_name, initial_or, parent_sym, child_sym)
        if ds is None:
            return

        st.subheader("3 · Validate data and optionally refine the OR")
        try:
            grain_ids, child_orientations = orientations_from_dataframe(ds.grains, convention=ds.convention)
            edges = edges_from_dataframe(ds.adjacency, grain_ids)
        except Exception as exc:
            st.error(f"Input validation failed: {exc}")
            return
        v1, v2, v3, v4 = st.columns(4)
        v1.metric("Daughter grains", len(grain_ids))
        v2.metric("Adjacency pairs", len(edges))
        v3.metric("Orientation convention", "crystal→specimen" if ds.convention == "crystal_to_specimen" else "specimen→crystal")
        v4.metric("Input provenance", ds.provenance)
        if len(grain_ids) > 600:
            st.warning("Candidate-graph methods can become expensive above several hundred grains in this transparent Python implementation. For production-scale EBSD maps, use a representative region or an optimized reconstruction package and use this app for auditable cross-checks.")

        refine = st.checkbox("Refine the selected OR near its starting value before reconstruction", value=False)
        r_or = initial_or
        if refine:
            bound = st.number_input("Maximum allowed OR correction — degrees", 0.1, 10.0, 4.0, 0.1)
            if st.button("Run bounded OR refinement"):
                with st.spinner("Refining OR against neighboring daughter-grain consistency…"):
                    rr, diag = refine_orientation_relationship(child_orientations, edges, initial_or, child_sym, parent_sym, max_rotation_deg=float(bound))
                st.session_state["wb_refined_or"] = rr
                st.session_state["wb_refined_or_diag"] = diag
                st.session_state["wb_refined_or_key"] = or_name
            if st.session_state.get("wb_refined_or_key") == or_name:
                diag = st.session_state.get("wb_refined_or_diag", {})
                r_or = st.session_state["wb_refined_or"]
                a, b, c = st.columns(3)
                a.metric("OR objective before", f"{diag.get('objective_before', np.nan):.5g}")
                b.metric("OR objective after", f"{diag.get('objective_after', np.nan):.5g}")
                c.metric("OR correction", f"{diag.get('correction_angle_deg', np.nan):.3f}°")
                st.caption("The refinement searches only a bounded neighborhood of the supplied physical OR. It is not an unconstrained OR-discovery algorithm.")

        st.subheader("4 · Run one or several reconstruction methods at the same time")
        selected = st.multiselect(
            "Reconstruction methods",
            METHODS,
            default=["Variant graph probability propagation", "Grain graph + Markov clustering", "Operator / groupoid consistency"],
            help="Select one method for a focused analysis or several/all methods for cross-method comparison on exactly the same input and OR.",
        )
        if not selected:
            st.warning("Select at least one reconstruction method.")
            return
        for m in selected:
            label, description, caveat = _method_explanation(m)
            with st.expander(f"{m} · {label}", expanded=False):
                st.write(description)
                st.caption(caveat)
        controls = _render_method_controls(selected)
        sig = _dataset_signature(ds, r_or, selected, controls)

        if st.button("Run selected reconstruction methods", type="primary", use_container_width=True):
            results: dict[str, ReconstructionResult] = {}
            runtimes: dict[str, float] = {}
            with st.spinner("Reconstructing parent orientations…"):
                for name in selected:
                    t0 = time.perf_counter()
                    results[name] = _run_one(name, grain_ids, child_orientations, edges, r_or, child_sym, parent_sym, controls)
                    runtimes[name] = time.perf_counter() - t0
            truth = ds.grains["true_parent_id"].to_numpy() if "true_parent_id" in ds.grains.columns else None
            comp = _comparison_table(results, runtimes, edges=edges, truth=truth)
            st.session_state["wb_results"] = results
            st.session_state["wb_comparison"] = comp
            st.session_state["wb_signature"] = sig
            st.session_state["wb_source"] = ds.grains.copy()
            st.session_state["wb_cycle_context"] = {
                "dataset_name": ds.name,
                "dataset_provenance": ds.provenance,
                "source_df": ds.grains.copy(),
                "grain_ids": np.asarray(grain_ids, int),
                "child_orientations": [np.asarray(g, float) for g in child_orientations],
                "r_child_to_parent": np.asarray(r_or, float),
                "orientation_relationship_name": or_name,
                "parent_phase": preset.parent_phase,
                "daughter_phase": preset.child_phase,
                "parent_symmetry": preset.parent_symmetry,
                "daughter_symmetry": preset.child_symmetry,
                "orientation_convention": ds.convention,
            }

        results = st.session_state.get("wb_results")
        comp = st.session_state.get("wb_comparison")
        if results and st.session_state.get("wb_signature") != sig:
            st.warning("Inputs, OR, selected methods or thresholds changed after the last run. Run reconstruction again before interpreting the old result.")
            results = None
            comp = None
        if results and comp is not None:
            st.subheader("5 · Validate first, then compare methods")
            truth = ds.grains["true_parent_id"].to_numpy() if "true_parent_id" in ds.grains.columns else None
            truth_validation_df = None
            if truth is not None:
                st.markdown("### A · Known-truth validation — this outranks cross-method agreement")
                truth_rows = []
                for method_name, rr in results.items():
                    tv = known_truth_validation_metrics(rr, truth, edges)
                    truth_rows.append({"method": method_name, **tv})
                truth_validation_df = pd.DataFrame(truth_rows)
                vcols = [
                    "method", "validation_status", "expected_parent_count", "reconstructed_parent_count",
                    "truth_ARI", "truth_homogeneity", "truth_completeness", "truth_V_measure",
                    "truth_boundary_precision", "truth_boundary_recall", "truth_boundary_F1",
                    "false_parent_boundary_rate", "singleton_parent_fraction",
                    "mean_reconstructed_fragments_per_true_parent",
                ]
                st.dataframe(truth_validation_df[vcols], use_container_width=True, hide_index=True)
                st.caption(
                    "ARI is chance-adjusted partition recovery. Homogeneity penalizes merging different true parents. Completeness penalizes fragmenting one true parent into many reconstructed parents. "
                    "V-measure balances homogeneity/completeness. Boundary precision asks whether called parent boundaries are truly boundaries; recall asks whether true boundaries were found; F1 balances both. "
                    "The PASS/REVIEW/FAIL label is a transparent synthetic-software checklist (ARI≥0.95, homogeneity/completeness≥0.95, boundary F1≥0.90, exact parent count, ≤10% singleton parents), not a universal experimental acceptance law."
                )
                if any(str(x).startswith("FAIL") for x in truth_validation_df["validation_status"]):
                    st.error("At least one method does NOT recover the known synthetic truth. Do not interpret near-zero OR residuals, high heuristic support, or method-to-method agreement as proof of a correct parent reconstruction.")
                elif all(str(x).startswith("PASS") for x in truth_validation_df["validation_status"]):
                    st.success("All selected methods recover the known synthetic truth under the currently selected crystallography and noise level.")

            st.markdown("### B · Method evidence table")
            st.dataframe(_comparison_display_table(comp), use_container_width=True, hide_index=True)
            st.caption(
                "Read across the table. Low OR residual is only an internal fit diagnostic: a singleton parent can have an almost-zero residual trivially. "
                "A defensible reconstruction should also avoid singleton splitting, recover sensible boundaries, remain stable across methods/thresholds, and agree with independent truth or retained-parent evidence when available. Runtime is not scientific quality."
            )

            agreement = None
            consensus = None
            if len(results) > 1:
                st.markdown("### C · Cross-method agreement — agreement is not accuracy")
                if truth is not None:
                    st.info("When known truth exists, use the truth table above to judge accuracy. ARI/NMI below only ask whether algorithms agree with EACH OTHER; several methods can agree perfectly on the same wrong partition.")
                agreement = build_agreement_summary(results, edges, parent_sym)
                consensus = boundary_consensus_table(results, edges, grain_ids)
                ctab1, ctab2, ctab3, ctab4 = st.tabs([
                    "Cluster agreement", "Boundary agreement", "Parent orientation agreement", "Boundary consensus"
                ])
                with ctab1:
                    st.markdown("**Adjusted Rand Index (ARI)**")
                    st.dataframe(agreement.ari.style.format("{:.3f}"), use_container_width=True)
                    st.caption("ARI = 1: identical parent partition up to arbitrary IDs. ARI ≈ 0: agreement near chance after correction. Negative values are possible for worse-than-chance partitions.")
                    st.markdown("**Normalized Mutual Information (NMI)**")
                    st.dataframe(agreement.nmi.style.format("{:.3f}"), use_container_width=True)
                    st.caption("NMI = 1: identical information about the partition; NMI = 0: no mutual clustering information. Unlike ARI, NMI is not chance-adjusted.")
                with ctab2:
                    st.markdown("**Prior-parent boundary Jaccard agreement**")
                    st.dataframe(agreement.boundary_jaccard.style.format("{:.3f}"), use_container_width=True)
                    st.caption("For the same daughter-grain adjacency graph, each method marks an edge as same-parent or prior-parent boundary. Jaccard = intersection/union of the boundary-edge sets; 1 means identical boundary calls.")
                with ctab3:
                    st.markdown("**Matched reconstructed-parent orientation disagreement — degrees**")
                    st.dataframe(agreement.matched_parent_orientation_deg.style.format("{:.3f}"), use_container_width=True)
                    st.caption("Parent clusters are matched by maximum shared daughter-grain overlap before comparing orientations modulo parent symmetry. Smaller values mean methods not only cluster similarly, but recover similar parent orientations.")
                with ctab4:
                    st.dataframe(consensus, use_container_width=True, hide_index=True)
                    st.caption("boundary_consensus_fraction = fraction of selected methods calling that daughter-grain adjacency edge a prior-parent boundary. Unanimous edges are the strongest cross-method boundary evidence; disagreements deserve map inspection.")

            st.markdown("### 6 · Inspect one method and its academic outputs")
            inspect = st.selectbox("Method to inspect", list(results))
            _render_result(
                inspect, results[inspect], ds.grains, comp,
                child_orientations=child_orientations, edges=edges, grain_ids=grain_ids,
                r_or=r_or, parent_sym=parent_sym, child_sym=child_sym, controls=controls,
            )

            with st.expander("Independent retained/reference parent validation", expanded=False):
                st.caption("If independently measured parent orientations are available, compare them here. This evidence is independent only because the uploaded anchors do not constrain the reconstruction itself.")
                rf = st.file_uploader("Reference parent orientation CSV", type=["csv"], key="ref_parent")
                if rf is not None:
                    rdf = _coerce_grain_columns(pd.read_csv(rf))
                    if "grain_id" not in rdf.columns:
                        rdf.insert(0, "grain_id", np.arange(1, len(rdf) + 1))
                    try:
                        rv = _reference_parent_validation(results[inspect], rdf, parent_sym, "crystal_to_specimen")
                        st.dataframe(rv, use_container_width=True, hide_index=True)
                        st.caption("misorientation_to_reference_deg is the nearest symmetry-reduced parent-orientation distance. Without a known spatial/ID correspondence it is a nearest-reference diagnostic, not proof that the IDs are the same physical grain.")
                    except Exception as exc:
                        st.error(str(exc))

            st.markdown("### 7 · Download a reproducible academic evidence bundle")
            parent_tables = {m: _parent_summary(r, ds.grains) for m, r in results.items()}
            daughter_tables = {m: _assignment_table(r, ds.grains) for m, r in results.items()}
            variant_tables = {m: variant_frequency_table(r, ds.grains) for m, r in results.items()}
            operator_tables = {}
            if "Operator / groupoid consistency" in results:
                operator_tables["Operator / groupoid consistency"] = operator_edge_table(
                    child_orientations, edges, grain_ids, r_or, parent_sym, child_sym,
                    float(controls.get("operator_tol_deg", 5.0)), result=results["Operator / groupoid consistency"]
                )
            metadata = {
                "dataset_name": ds.name,
                "input_provenance": ds.provenance,
                "graph_kind": ds.graph_kind,
                "map_length_unit": ds.length_unit,
                "input_notes": ds.input_notes,
                "orientation_convention": ds.convention,
                "parent_symmetry": preset.parent_symmetry,
                "daughter_symmetry": preset.child_symmetry,
                "orientation_relationship_name": or_name,
                "orientation_relationship_matrix_daughter_to_parent": np.asarray(r_or, float).tolist(),
                "orientation_relationship_was_refined": bool(refine and st.session_state.get("wb_refined_or_key") == or_name),
                "selected_methods": list(results.keys()),
                "controls": controls,
                "daughter_grains": int(len(grain_ids)),
                "adjacency_pairs": int(len(edges)),
                "known_truth_available": bool(truth is not None),
                "known_truth_primary_metrics": "ARI, homogeneity, completeness, V-measure, boundary precision/recall/F1" if truth is not None else None,
                "cross_method_agreement_warning": "ARI/NMI between methods measure agreement, not truth accuracy",
                "support_score_status": "deterministic heuristic, not calibrated probability",
                "variant_id_status": "internal symmetry-generated candidate ID unless an external canonical mapping is supplied",
            }
            methods_text = (
                "PARENT/DAUGHTER RECONSTRUCTION EVIDENCE BUNDLE\n\n"
                f"Dataset: {ds.name}\nInput provenance: {ds.provenance}\nAdjacency: {ds.graph_kind}\n"
                f"Orientation convention: {ds.convention}\nParent symmetry: {preset.parent_symmetry}\nDaughter symmetry: {preset.child_symmetry}\n"
                f"Orientation relationship: {or_name}\nMethods: {', '.join(results.keys())}\nControls: {json.dumps(controls, sort_keys=True)}\n\n"
                "Interpretation rules: OR residuals are angular internal-consistency diagnostics; the 0–1 support value is heuristic and is not a probability. "
                "Variant IDs are internal candidate IDs unless a canonical transformation-specific mapping is supplied. When known truth exists, truth ARI/homogeneity/completeness/V-measure and boundary precision/recall/F1 are the primary validation evidence. Cross-method ARI/NMI measure agreement, not accuracy; boundary Jaccard compares prior-parent boundary calls, and matched-parent orientation disagreement compares recovered parent orientations after overlap-based matching.\n"
            )
            bundle = academic_export_zip(
                source_df=ds.grains, adjacency_df=ds.adjacency, results=results, comparison=comp,
                agreement=agreement, parent_tables=parent_tables, daughter_tables=daughter_tables,
                variant_tables=variant_tables, operator_tables=operator_tables or None,
                metadata=metadata, methods_text=methods_text, boundary_consensus=consensus,
                truth_validation=truth_validation_df,
            )
            st.download_button(
                "Download complete reconstruction evidence ZIP", bundle,
                "parent_daughter_reconstruction_evidence.zip", "application/zip", use_container_width=True
            )
            st.caption("The ZIP contains the analyzed input grains/adjacency, exact OR matrix, all method controls, comparison matrices, parent summaries, daughter assignments, variant/operator statistics, and machine-readable metadata. When synthetic/known truth exists it also contains validation/known_truth_metrics.csv with ARI, homogeneity, completeness, V-measure and boundary precision/recall/F1.")

    with tab_cycle:
        st.subheader("Round-trip crystallography: measured B19′ → reconstructed B2 → regenerated B19′")
        st.markdown(
            "This workflow takes a completed parent reconstruction and immediately transforms each reconstructed parent orientation forward again. "
            "For NiTi it answers two different questions: **(1)** does the measured B19′ close consistently through the reconstructed B2 under the same OR, and **(2)** what B19′ orientation branches are crystallographically allowed if that B2 parent transforms again?"
        )
        st.warning(
            "Scientific scope guard: round-trip closure is an internal consistency check, not independent validation, because the same OR is used in both directions. "
            "The regenerated branch library does not predict which future B19′ variant nucleates. A separately measured later-cycle EBSD map provides independent experimental evidence and can be matched below."
        )

        cycle_results = st.session_state.get("wb_results")
        ctx = st.session_state.get("wb_cycle_context")
        if not cycle_results or not ctx:
            st.info("Run a Daughter → parent reconstruction first. The cycle tab will then reuse the exact reconstructed parents, exact OR matrix, symmetries and measured daughter orientations automatically—no retyping.")
        else:
            method = st.selectbox("Parent reconstruction to use for the cycle", list(cycle_results), key="cycle_reconstruction_method")
            result = cycle_results[method]
            rcp = np.asarray(ctx["r_child_to_parent"], float)
            psg = symmetry_group(str(ctx["parent_symmetry"]))
            csg = symmetry_group(str(ctx["daughter_symmetry"]))
            gids = np.asarray(ctx["grain_ids"], int)
            measured = [np.asarray(g, float) for g in ctx["child_orientations"]]
            source_df = ctx["source_df"].copy()

            is_niti = str(ctx.get("parent_phase", "")).startswith("B2") and "B19" in str(ctx.get("daughter_phase", ""))
            phase_label = "B19′ → B2 → B19′" if is_niti else "daughter → parent → daughter"
            a, b, c, d = st.columns(4)
            a.metric("Cycle", phase_label)
            b.metric("Reconstruction method", method)
            c.metric("Reconstructed parents", len(result.parent_orientations))
            d.metric("Measured daughter grains", len(gids))
            st.caption(
                f"OR reused exactly: {ctx['orientation_relationship_name']} · parent symmetry: {ctx['parent_symmetry']} · daughter symmetry: {ctx['daughter_symmetry']}. "
                "Reusing the identical R_cp prevents a hidden convention change between reconstruction and regeneration."
            )

            st.markdown("### Math used by the round-trip")
            st.latex(r"g_D^{(k)} = g_P\,S_P^{(k)}\,R_{cp}")
            st.caption("g_P = reconstructed parent orientation; S_P^(k) = k-th proper parent-symmetry operation; R_cp = daughter→parent crystal-frame OR rotation; g_D^(k) = regenerated daughter orientation branch k. All are proper rotations except k, which is a discrete branch index.")
            st.latex(r"\delta_i = \min_k d_{G_D}\!\left(g_{D,i}^{\mathrm{meas}},\,g_{D,P_i}^{(k)}\right)")
            st.caption("δ_i = round-trip cycle-closure misorientation for measured daughter grain i, in degrees; d_GD = minimum disorientation modulo daughter crystal symmetry G_D; P_i = reconstructed parent assigned to grain i. Lower δ_i means stronger internal OR/parent consistency.")
            st.latex(r"k_i^* = \operatorname*{arg\,min}_k d_{G_D}\!\left(g_{D,i}^{\mathrm{meas}},\,g_{D,P_i}^{(k)}\right)")
            st.caption("k_i* = regenerated branch that best reproduces the measured daughter orientation. It is an internal branch ID for this exact OR/symmetry enumeration, not automatically a canonical literature variant number.")
            st.latex(r"\Delta\theta_{j\rightarrow k}=d_{G_D}\!\left(g_D^{(j)},g_D^{(k)}\right)")
            st.caption("Δθ_j→k = symmetry-reduced orientation change between two allowed regenerated daughter branches, in degrees. It describes geometric reorientation only; it is not a transition probability or energy barrier.")

            st.markdown("### Interpretation thresholds")
            t1, t2, t3 = st.columns(3)
            strong = t1.number_input("Strong closure ≤ — °", min_value=0.05, max_value=5.0, value=1.0, step=0.05, key="cycle_strong")
            acceptable = t2.number_input("Acceptable closure ≤ — °", min_value=0.10, max_value=10.0, value=2.5, step=0.1, key="cycle_acceptable")
            review = t3.number_input("Review threshold ≤ — °", min_value=0.50, max_value=20.0, value=5.0, step=0.25, key="cycle_review")
            st.caption("These are interpretation thresholds, not material constants. Enforced UI ranges are 0.05–5°, 0.10–10° and 0.50–20°. Report the chosen values in Methods and do not compare studies using different thresholds without stating them.")
            if not (float(strong) <= float(acceptable) <= float(review)):
                st.error("Thresholds must satisfy: strong ≤ acceptable ≤ review.")
            else:
                library, _ = regenerated_variant_library(result.parent_orientations, rcp, psg, csg)
                closure = cycle_closure_table(
                    gids, measured, result, rcp, psg, csg,
                    strong_deg=float(strong), acceptable_deg=float(acceptable), review_deg=float(review),
                )
                psummary = parent_cycle_summary(closure, library)
                occupancy = observed_branch_occupancy(closure, source_df)
                switch = branch_switch_catalog(result.parent_orientations, rcp, psg, csg)

                q1, q2, q3, q4 = st.columns(4)
                q1.metric("Allowed regenerated branches / parent", int(library.groupby("parent_id").size().median()))
                q2.metric("Mean cycle closure", f"{closure['cycle_closure_misorientation_deg'].mean():.3f}°")
                q3.metric("P95 cycle closure", f"{closure['cycle_closure_misorientation_deg'].quantile(0.95):.3f}°")
                q4.metric("Strong-closure grains", f"{100.0*np.mean(closure['cycle_closure_misorientation_deg'] <= float(strong)):.1f}%")
                st.caption("For NiTi B2 parent symmetry and B19′ monoclinic daughter symmetry, the natural/AQ OR produces 12 symmetry-distinct regenerated B19′ orientation branches per reconstructed B2 parent.")

                ctab1, ctab2, ctab3, ctab4, ctab5 = st.tabs([
                    "Round-trip closure", "Regenerated B19′ library", "Observed branch use", "Possible branch switching", "Measure a later B19′ cycle"
                ])
                new_cycle_matches = None
                with ctab1:
                    st.markdown("**Parent-level round-trip summary**")
                    st.dataframe(psummary, use_container_width=True, hide_index=True)
                    st.caption("The strongest academic readout is the distribution of δ_i, not a single opaque score. Mean/median show typical closure; P95 and maximum expose bad tails; branch coverage tells how much of the allowed orientation library is actually represented by the measured daughters.")
                    st.markdown("**Daughter-grain round-trip evidence**")
                    show = source_df.merge(closure, left_on="grain_id", right_on="daughter_grain_id", how="left") if "grain_id" in source_df.columns else closure
                    st.dataframe(show, use_container_width=True, hide_index=True)
                    if {"x", "y"}.issubset(show.columns):
                        fig = px.scatter(
                            show, x="x", y="y", color="cycle_closure_misorientation_deg",
                            hover_data=["daughter_grain_id", "reconstructed_parent_id", "observed_regenerated_branch_id", "cycle_closure_quality"],
                            title="Round-trip B19′ → B2 → B19′ closure map",
                            labels={"cycle_closure_misorientation_deg": "cycle closure δ (°)", "x": "map x", "y": "map y"},
                        )
                        fig.update_yaxes(scaleanchor="x", scaleratio=1)
                        st.plotly_chart(fig, use_container_width=True)
                    st.download_button("Download round-trip closure CSV", closure.to_csv(index=False), "b19p_b2_b19p_cycle_closure.csv", "text/csv")

                with ctab2:
                    st.dataframe(library, use_container_width=True, hide_index=True)
                    st.caption("One row = one symmetry-distinct regenerated daughter orientation for one reconstructed parent. Euler angles, quaternion and all 3×3 orientation-matrix entries are exported so the library can be reproduced or compared in MTEX/other crystallographic software without losing orientation precision.")
                    st.download_button("Download regenerated daughter library", library.to_csv(index=False), "regenerated_b19p_variant_library.csv", "text/csv")

                with ctab3:
                    st.dataframe(occupancy, use_container_width=True, hide_index=True)
                    st.caption("Observed branch occupancy assigns each currently measured daughter grain to the regenerated branch with the smallest symmetry-reduced angular residual. Grain fractions are always shown; area fractions appear when a reliable grain-area column exists.")
                    st.download_button("Download observed branch occupancy", occupancy.to_csv(index=False), "observed_regenerated_branch_occupancy.csv", "text/csv")

                with ctab4:
                    pid = st.selectbox("Reconstructed parent for switch matrix", sorted(result.parent_orientations), key="cycle_switch_parent")
                    psw = switch[switch["parent_id"] == int(pid)].copy()
                    matrix = psw.pivot(index="from_regenerated_branch_id", columns="to_regenerated_branch_id", values="daughter_orientation_change_deg")
                    st.dataframe(matrix.style.format("{:.3f}"), use_container_width=True)
                    st.caption("Entry (j,k) = Δθ_j→k in degrees. The diagonal is 0° by definition; off-diagonal entries are allowed daughter-orientation changes if the same reconstructed parent forms a different branch. This matrix says nothing about which switch is energetically preferred.")
                    st.dataframe(psw, use_container_width=True, hide_index=True)
                    st.download_button("Download branch-switch catalogue", switch.to_csv(index=False), "regenerated_branch_switch_catalog.csv", "text/csv")

                with ctab5:
                    st.markdown("**Optional independent experiment: upload a separately measured later-cycle daughter map**")
                    st.caption("This is stronger than round-trip self-consistency because the later-cycle orientations were not used to reconstruct the original parent. Use a grain-level CSV when already segmented, or upload raw ANG/CTF and segment it here.")
                    nmode = st.radio("Later-cycle input", ["Grain-level CSV", "Raw .ang / .ctf EBSD"], horizontal=True, key="cycle_new_mode")
                    new_grains = None
                    new_cycle_convention = "crystal_to_specimen"
                    if nmode == "Grain-level CSV":
                        nf = st.file_uploader("Later-cycle daughter grain CSV", type=["csv"], key="cycle_new_csv")
                        conv = st.radio("Later-cycle grain orientation convention", ["Crystal → specimen", "Specimen → crystal"], horizontal=True, key="cycle_new_conv")
                        new_cycle_convention = "crystal_to_specimen" if conv == "Crystal → specimen" else "specimen_to_crystal"
                        if nf is not None:
                            new_grains = _coerce_grain_columns(pd.read_csv(nf))
                            if "grain_id" not in new_grains.columns:
                                new_grains.insert(0, "grain_id", np.arange(1, len(new_grains)+1))
                    else:
                        nf = st.file_uploader("Later-cycle raw EBSD", type=["ang", "ctf"], key="cycle_new_raw")
                        if nf is not None:
                            suffix = Path(nf.name).suffix.lower()
                            points = _read_ang_bytes(nf.getvalue()) if suffix == ".ang" else _read_ctf_bytes(nf.getvalue())
                            st.markdown("**Later-cycle EBSD import audit**")
                            st.dataframe(_raw_ebsd_audit(points), use_container_width=True, hide_index=True)
                            st.caption(f"File: {nf.name} · detected format: {points.attrs.get('format')} · Euler-unit handling: {points.attrs.get('euler_input_unit', 'see format')}. The importer fails closed when required columns cannot be identified.")
                            phase_series = pd.to_numeric(points["phase"], errors="coerce") if "phase" in points.columns else pd.Series(np.ones(len(points)), index=points.index, dtype=float)
                            phases = sorted(float(x) for x in phase_series.dropna().unique() if float(x) > 0)
                            if not phases:
                                st.error("No indexed phase (>0) is present in the later-cycle map.")
                            else:
                                phase_counts = pd.DataFrame({
                                    "phase_id": phases,
                                    "indexed_points": [int(np.sum(phase_series == x)) for x in phases],
                                })
                                st.dataframe(phase_counts, use_container_width=True, hide_index=True)
                                phase = st.selectbox("Later-cycle B19′ phase code", phases, key="cycle_new_phase")
                                points = points.loc[phase_series == float(phase)].copy().reset_index(drop=True)
                                with st.expander("Optional later-cycle EBSD quality filters — OFF by default", expanded=False):
                                    st.caption("Quality fields are vendor/acquisition dependent and are never filtered silently. Any enabled threshold must be reported in Methods.")
                                    if "CI" in points.columns:
                                        use_ci = st.checkbox("Minimum EDAX confidence index (CI)", value=False, key="cycle_new_use_ci")
                                        ci_min = st.number_input("Later-cycle minimum CI", min_value=-1.0, max_value=1.0, value=0.05, step=0.01, disabled=not use_ci, key="cycle_new_ci")
                                        if use_ci:
                                            points = points[pd.to_numeric(points["CI"], errors="coerce") >= float(ci_min)].copy()
                                    if "MAD" in points.columns:
                                        use_mad = st.checkbox("Maximum Oxford mean angular deviation (MAD)", value=False, key="cycle_new_use_mad")
                                        mad_max = st.number_input("Later-cycle maximum MAD — °", min_value=0.0, max_value=10.0, value=1.5, step=0.1, disabled=not use_mad, key="cycle_new_mad")
                                        if use_mad:
                                            points = points[pd.to_numeric(points["MAD"], errors="coerce") <= float(mad_max)].copy()
                                    if "fit" in points.columns:
                                        use_fit = st.checkbox("Maximum ANG fit field", value=False, key="cycle_new_use_fit")
                                        fit_max = st.number_input("Later-cycle maximum ANG fit", min_value=0.0, max_value=10.0, value=2.0, step=0.1, disabled=not use_fit, key="cycle_new_fit")
                                        if use_fit:
                                            points = points[pd.to_numeric(points["fit"], errors="coerce") <= float(fit_max)].copy()
                                if len(points) < 2:
                                    st.error("Fewer than two indexed later-cycle B19′ points remain after phase/quality filtering.")
                                else:
                                    s1, s2, s3 = st.columns(3)
                                    seg_deg = s1.number_input("Later-cycle grain-boundary disorientation — °", min_value=0.5, max_value=15.0, value=3.0, step=0.25, key="cycle_new_seg")
                                    min_pts = s2.number_input("Minimum indexed points / grain", min_value=2, max_value=1000, value=5, step=1, key="cycle_new_minpts")
                                    radius_factor = s3.number_input("Spatial neighbor radius / median spacing", min_value=1.05, max_value=1.70, value=1.35, step=0.05, key="cycle_new_radius")
                                    st.caption("Allowed ranges: segmentation 0.5–15°, minimum size 2–1000 points, neighbor-radius factor 1.05–1.70. Defaults match the auditable raw-EBSD route in the parent reconstruction tab; all three must be reported.")
                                    confirmed = st.checkbox("I verified the later-cycle EBSD Euler/specimen reference frame against the original reconstruction", value=False, key="cycle_new_frame_confirm")
                                    if confirmed:
                                        try:
                                            new_grains, _ = _segment_points_to_grains(points, csg, float(seg_deg), int(min_pts), neighbor_radius_factor=float(radius_factor))
                                            st.success(f"Later-cycle EBSD segmented into {len(new_grains)} B19′ daughter grains.")
                                            st.dataframe(new_grains.head(200), use_container_width=True, hide_index=True)
                                        except Exception as exc:
                                            st.error(f"Later-cycle segmentation failed: {exc}")
                    if new_grains is not None and len(new_grains):
                        try:
                            ngids, nori = orientations_from_dataframe(new_grains, convention=new_cycle_convention)
                            known = new_grains["parent_id"].tolist() if "parent_id" in new_grains.columns else None
                            new_cycle_matches = match_new_cycle_daughters(ngids, nori, result.parent_orientations, rcp, psg, csg, known_parent_ids=known)
                            nshow = new_grains.merge(new_cycle_matches, left_on="grain_id", right_on="new_cycle_daughter_grain_id", how="left")
                            r1, r2, r3 = st.columns(3)
                            r1.metric("Later-cycle grains matched", len(new_cycle_matches))
                            r2.metric("Mean later-cycle library fit", f"{new_cycle_matches['new_cycle_OR_library_fit_deg'].mean():.3f}°")
                            r3.metric("P95 later-cycle library fit", f"{new_cycle_matches['new_cycle_OR_library_fit_deg'].quantile(0.95):.3f}°")
                            st.dataframe(nshow, use_container_width=True, hide_index=True)
                            if known is None:
                                st.warning("No parent_id column was supplied, so parent assignment is orientation-only across all reconstructed parents. Use the best-vs-second-parent separation column to detect ambiguity; spatially known parent IDs are stronger evidence.")
                            else:
                                st.success("parent_id was supplied for the later-cycle grains, so each grain was tested against its known reconstructed parent rather than choosing a parent by orientation alone.")
                            st.download_button("Download later-cycle matches", new_cycle_matches.to_csv(index=False), "later_cycle_b19p_matches.csv", "text/csv")
                        except Exception as exc:
                            st.error(f"Could not match later-cycle daughter orientations: {exc}")

                st.markdown("### Academic export")
                metadata = {
                    "workflow": phase_label,
                    "dataset_name": ctx.get("dataset_name"),
                    "parent_reconstruction_method": method,
                    "orientation_relationship_name": ctx.get("orientation_relationship_name"),
                    "orientation_relationship_matrix_daughter_to_parent": rcp.tolist(),
                    "parent_symmetry": ctx.get("parent_symmetry"),
                    "daughter_symmetry": ctx.get("daughter_symmetry"),
                    "strong_closure_threshold_deg": float(strong),
                    "acceptable_closure_threshold_deg": float(acceptable),
                    "review_closure_threshold_deg": float(review),
                    "round_trip_validation_status": "internal consistency; not independent validation",
                    "future_variant_selection_status": "allowed orientation library only; no nucleation probability claimed",
                }
                bundle = cycle_evidence_zip(
                    variant_library=library, closure=closure, parent_summary=psummary,
                    occupancy=occupancy, switch_catalog=switch, metadata=metadata,
                    new_cycle_matches=new_cycle_matches,
                )
                st.download_button("Download complete B19′→B2→B19′ cycle evidence ZIP", bundle, "b19p_b2_b19p_cycle_evidence.zip", "application/zip", use_container_width=True)
                st.caption("The ZIP contains the exact regenerated orientation library, round-trip closure table, parent summary, observed branch occupancy, full branch-switch catalogue, metadata/thresholds and later-cycle matches when supplied.")

    with tab_forward:
        st.subheader("Forward crystallography: known parent orientation(s) → all possible daughter variants")
        fname, fp, rcp, psg, csg = _render_or_builder("forward")
        fmode = st.radio("Parent-orientation input", ["Manual table", "Upload CSV"], horizontal=True)
        if fmode == "Manual table":
            parents = st.data_editor(pd.DataFrame([
                {"grain_id": 1, "phi1_deg": 0.0, "Phi_deg": 0.0, "phi2_deg": 0.0},
            ]), num_rows="dynamic", use_container_width=True, key="forward_parent_editor")
        else:
            pf = st.file_uploader("Parent orientation CSV", type=["csv"], key="forward_parent_file")
            if pf is None:
                st.info("Upload parent orientations to continue.")
                parents = None
            else:
                parents = _coerce_grain_columns(pd.read_csv(pf))
        if parents is not None:
            try:
                ids, pori = orientations_from_dataframe(parents, convention="crystal_to_specimen")
                rows = []
                for pid, gp in zip(ids, pori):
                    variants = unique_child_variants(gp, rcp, psg, csg)
                    for vid, g in enumerate(variants, 1):
                        p1, P, p2 = matrix_to_bunge_euler(g)
                        qw, qx, qy, qz = quaternion_wxyz(g)
                        rows.append({
                            "parent_id": int(pid), "daughter_variant_id": int(vid),
                            "daughter_phi1_deg": p1, "daughter_Phi_deg": P, "daughter_phi2_deg": p2,
                            "daughter_qw": qw, "daughter_qx": qx, "daughter_qy": qy, "daughter_qz": qz,
                        })
                out = pd.DataFrame(rows)
                a, b = st.columns(2)
                a.metric("Parent orientations supplied", len(ids))
                b.metric("Total daughter variants generated", len(out))
                st.dataframe(out, use_container_width=True, hide_index=True)
                st.caption("Forward output gives crystallographically allowed daughter orientations only. It does NOT predict which variant nucleates, its volume fraction, morphology or spatial position; those require additional thermomechanical/microstructural physics.")
                st.download_button("Download predicted daughter variants", out.to_csv(index=False), "predicted_daughter_variants.csv", "text/csv")
            except Exception as exc:
                st.error(str(exc))

    with tab_batch:
        st.divider()
        st.subheader("Batch compare several already-segmented daughter datasets")
        st.caption("One step beyond a single-map workflow: upload several grain-level CSV files at once, apply the same OR and one reconstruction method, and compare map-level diagnostics. Each CSV must contain grain_id, orientation columns and x/y so the same explicitly approximate k-NN adjacency can be built unless you analyze it separately with measured adjacency in the main tab.")
        bname, bp, br, bps, bcs = _render_or_builder("batch")
        files = st.file_uploader("Grain-level CSV datasets", type=["csv"], accept_multiple_files=True, key="batch_files")
        method = st.selectbox("Batch reconstruction method", METHODS, index=2)
        k = st.number_input("Batch centroid k-NN adjacency — k", 1, 12, 4, 1)
        sigma = st.number_input("Batch primary angular width / tolerance — degrees", 0.2, 15.0, 2.5, 0.1)
        if files and st.button("Run batch reconstruction"):
            rows = []
            for f in files:
                try:
                    df = _coerce_grain_columns(pd.read_csv(f))
                    gids, ori = orientations_from_dataframe(df)
                    idx_edges = approximate_knn_edges(df, int(k))
                    controls = {
                        "sigma_deg": float(sigma), "inflation": 1.6, "merge_deg": 8.0,
                        "nucleation_deg": 3.0, "growth_deg": 8.0,
                        "operator_tol_deg": 8.0, "parent_consistency_deg": 5.0,
                    }
                    t0 = time.perf_counter()
                    rr = _run_one(method, gids, ori, idx_edges, br, bcs, bps, controls)
                    rows.append({
                        "dataset": f.name,
                        "daughter_grains": len(gids),
                        "reconstructed_parents": rr.diagnostics["n_reconstructed_parents"],
                        "mean_OR_fit_deg": float(np.mean(rr.fit_deg)),
                        "p95_OR_fit_deg": float(np.quantile(rr.fit_deg, 0.95)),
                        "mean_support_score": float(np.mean(rr.confidence)),
                        "low_support_fraction_pct": float(100*np.mean((rr.confidence < 0.35) | (rr.fit_deg > 5.0))),
                        "runtime_s": time.perf_counter() - t0,
                    })
                except Exception as exc:
                    rows.append({"dataset": f.name, "error": str(exc)})
            bout = pd.DataFrame(rows)
            st.dataframe(bout, use_container_width=True, hide_index=True)
            st.download_button("Download batch reconstruction summary", bout.to_csv(index=False), "batch_reconstruction_summary.csv", "text/csv")
            st.warning("Batch mode uses centroid k-NN adjacency for comparability and convenience. For publication-quality boundary topology, analyze each dataset in the main tab with measured adjacency when available.")

    with tab_dictionary:
        st.subheader("What every reconstruction input and output means")
        dictionary = pd.DataFrame([
            ["grain_id", "Unique daughter-grain identifier", "count / ID", "Required"],
            ["x, y", "Daughter-grain centroid coordinates used for maps or approximate adjacency", "usually µm", "Recommended"],
            ["phi1_deg, Phi_deg, phi2_deg", "Bunge Euler orientation angles of a grain", "degrees", "Required unless quaternion is supplied"],
            ["qw, qx, qy, qz", "Unit quaternion describing the same grain orientation", "dimensionless", "Alternative to Euler angles"],
            ["grain_id_1, grain_id_2", "Two daughter grains that physically share a boundary", "IDs", "Required for graph methods; k-NN fallback possible"],
            ["OR / R_cp", "Orientation relationship; daughter→parent crystal-frame rotation used to generate candidate parents", "rotation / degrees", "Required or supplied by preset"],
            ["parent_id", "Reconstructed prior-parent cluster assigned to a daughter grain", "ID", "Primary output"],
            ["candidate_variant_id", "Index of the crystallographically allowed parent candidate selected for that daughter grain", "ID", "Primary output"],
            ["best_OR_fit_deg", "Angular mismatch between selected candidate and representative reconstructed parent", "degrees", "Primary quality metric; lower is better"],
            ["second_best_candidate_fit_deg", "Angular fit of the next-best crystallographic parent candidate", "degrees", "Used to expose ambiguity"],
            ["candidate_separation_deg", "Second-best fit minus best fit", "degrees", "Larger separation = less candidate ambiguity"],
            ["support_score_0_to_1", "Deterministic heuristic support built from absolute fit and candidate separation", "0–1", "Not a probability"],
            ["parent_phi1_deg, parent_Phi_deg, parent_phi2_deg", "Reconstructed parent orientation in Bunge Euler coordinates", "degrees", "Parent-level output"],
            ["parent_qw, parent_qx, parent_qy, parent_qz", "Same reconstructed parent orientation as a unit quaternion", "dimensionless", "Parent-level output"],
            ["operator_residual_deg", "Distance of a measured daughter-neighbor misorientation to its nearest theoretical transformation operator", "degrees", "Operator/groupoid route; lower is better"],
            ["cycle closure δ_i", "Minimum daughter-symmetry-reduced misorientation between measured daughter grain i and all daughter orientations regenerated from its reconstructed parent", "degrees", "Lower = stronger round-trip internal consistency; not independent validation"],
            ["regenerated branch k", "One symmetry-distinct daughter orientation generated from a reconstructed parent using the exact OR", "discrete ID", "Internal enumeration; not automatically a canonical literature variant number"],
            ["Δθ_j→k", "Daughter-symmetry-reduced orientation change between regenerated branches j and k", "degrees", "Geometry of an allowed branch change; not a switching probability"],
            ["later-cycle library fit", "Minimum misorientation between an independently measured later-cycle daughter and the regenerated library of a reconstructed parent", "degrees", "Independent experimental comparison when the later-cycle EBSD was not used in the original reconstruction"],
            ["ARI", "Adjusted Rand Index comparing parent partitions from two methods", "dimensionless", "1 = identical clustering up to labels; chance-adjusted"],
            ["NMI", "Normalized Mutual Information comparing parent partitions", "dimensionless", "1 = identical clustering information; not chance-adjusted"],
            ["boundary Jaccard", "Intersection/union agreement between prior-parent boundary-edge sets", "dimensionless", "1 = identical boundary calls on the same adjacency graph"],
            ["matched parent orientation disagreement", "Misorientation between overlap-matched reconstructed parent orientations from two methods", "degrees", "Smaller = stronger orientation agreement"],
        ], columns=["field / symbol", "meaning", "unit", "how to use it"])
        st.dataframe(dictionary, use_container_width=True, hide_index=True)
        st.markdown(
            "**Minimum defensible Methods report:** daughter phase and symmetry; parent phase and symmetry; OR definition and whether refined; orientation convention; grain segmentation method/threshold; adjacency source; selected reconstruction algorithm(s); every angular/graph threshold; number of daughter grains; reconstructed parent count; OR-fit statistics; confidence/support definition; and any retained-parent or known-truth validation."
        )
        st.info("These implementations are transparent members of established reconstruction families, not binary-identical reproductions of MTEX, ARPGE, OIM Analysis, AZtecCrystal or other vendor/research packages.")

    with tab_sources:
        st.divider()
        st.subheader("How this workbench relates to established reconstruction tools")
        st.markdown(
            "The app intentionally combines several established **families** on the same daughter-grain dataset. "
            "It does not claim binary identity with external software; instead it exposes the crystallographic inputs, thresholds, residuals and cross-method disagreements needed for an auditable research comparison."
        )
        basis = pd.DataFrame([
            {
                "our route": "Operator / groupoid consistency",
                "closest established idea": "ARPGE / GenOVa operator-groupoid reconstruction",
                "core evidence": "daughter–daughter theoretical operator match + common-parent consistency",
                "extra output here": "per-edge operator ID/residual, operator frequencies, parent/daughter tables",
                "reference": "J. Appl. Cryst. 40 (2007) 1183–1188, DOI 10.1107/S0021889807048777",
            },
            {
                "our route": "Grain graph + Markov clustering",
                "closest established idea": "MTEX parent grain graph",
                "core evidence": "weighted compatibility of neighboring daughter grains followed by graph clustering",
                "extra output here": "cross-method ARI/NMI, boundary Jaccard and orientation disagreement",
                "reference": "J. Appl. Cryst. 55 (2022), DOI 10.1107/S1600576721011560",
            },
            {
                "our route": "Variant graph probability propagation",
                "closest established idea": "variant-graph parent reconstruction",
                "core evidence": "multiple parent candidates retained for each daughter grain before graph resolution",
                "extra output here": "candidate separation + explicit heuristic-support decomposition",
                "reference": "Materialia 22 (2022) 101399, DOI 10.1016/j.mtla.2022.101399",
            },
            {
                "our route": "Neighbor voting",
                "closest established idea": "local neighbor-level voting",
                "core evidence": "local candidate-parent consensus across daughter neighbors",
                "extra output here": "same academic tables and agreement diagnostics as graph methods",
                "reference": "Compared with graph/variant-graph approaches in Materialia 22 (2022) 101399",
            },
            {
                "our route": "Nucleation + growth",
                "closest established idea": "strict seed formation followed by looser spatial growth",
                "core evidence": "high-confidence parent seeds + controlled expansion to neighboring daughter grains",
                "extra output here": "explicit nucleation/growth thresholds and cross-method validation",
                "reference": "Historical parent reconstruction workflow; thresholds must be dataset-specific",
            },
        ])
        st.dataframe(basis, use_container_width=True, hide_index=True)
        st.markdown(
            "**NiTi Correspondence Theory (CT) + Otsuka–Ren route:** the source-derived NiTi model uses the measured B2/B19′ metric tensors and the Otsuka–Ren correspondence matrix `C^{A→M}=[[0,1,-1],[0,1,1],[1,0,0]]`. CT then uses correspondence + metrics + symmetries for correspondence variants, intercorrespondence/operator classes, twins and CMC/SMC compatibility. **This is not itself an EBSD reconstruction algorithm or a unique experimental OR.** For the reconstruction engine, the dedicated CT option transparently uses the polar rotational factor of the metric-aware correspondence deformation only as a model-derived starting OR; an experimental/refined OR remains preferable when available."
        )
        st.markdown(
            "**NiTi natural/AQ orientation route:** the separate metric-aware quick setup uses (010)B19′ ∥ (110)B2 and [101]B19′ ∥ [−1 1 1]B2 and converts the monoclinic direct/reciprocal lattice vectors to Cartesian crystal coordinates before constructing the OR. It remains available beside the CT/Otsuka–Ren route rather than being replaced by it. The cycle module reuses whichever exact OR was selected for reconstruction."
        )
        st.markdown(
            "**External-tool lesson built into this page:** dedicated EBSD packages normally import vendor metadata, verify reference frames, segment grains, fit/refine the OR, reconstruct parents, then transfer parent IDs/variant IDs/fit back to pixels. "
            "This app follows the same separation of concerns and refuses to hide raw-EBSD phase selection, reference-frame confirmation, segmentation, adjacency source or OR refinement."
        )
        st.markdown(
            "**Recommended academic evidence hierarchy:** (1) verify EBSD reference frame and phase indexing; (2) report segmentation and adjacency; (3) define/refine the OR; (4) inspect parent and daughter residual tables; "
            "(5) compare independent reconstruction families; (6) validate against retained parent/known truth when possible; (7) report variant/operator statistics only with their exact ID convention."
        )

