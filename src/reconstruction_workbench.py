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
    reconstruction_accuracy_against_labels,
    refine_orientation_relationship,
    rotation_from_parallelisms,
    symmetry_group,
    synthetic_parent_reconstruction_demo,
    unique_child_variants,
    variant_graph_reconstruction,
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
    """Read the standard numeric part of a TSL/OIM .ang file.

    The first three columns of conventional .ang data are Euler angles in radians;
    x and y follow.  This parser intentionally returns point data and does not invent
    grain IDs.  Grain segmentation is a separate explicit step in the UI.
    """
    text = raw.decode("utf-8", errors="ignore")
    rows: list[list[float]] = []
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.replace(",", " ").split()
        try:
            vals = [float(x) for x in parts]
        except ValueError:
            continue
        if len(vals) >= 5:
            rows.append(vals)
    if not rows:
        raise ValueError("No numeric EBSD rows were found in the .ang file.")
    arr = np.asarray(rows, float)
    out = pd.DataFrame({
        "phi1_deg": np.degrees(arr[:, 0]),
        "Phi_deg": np.degrees(arr[:, 1]),
        "phi2_deg": np.degrees(arr[:, 2]),
        "x": arr[:, 3],
        "y": arr[:, 4],
    })
    if arr.shape[1] > 7:
        out["phase"] = arr[:, 7].astype(int)
    return out


def _read_ctf_bytes(raw: bytes) -> pd.DataFrame:
    """Read common Oxford/HKL .ctf point columns (Euler angles are in degrees)."""
    text = raw.decode("utf-8", errors="ignore")
    lines = text.splitlines()
    header_idx = None
    for i, line in enumerate(lines):
        l = line.strip().lower()
        if l.startswith("phase") and "euler1" in l and "euler2" in l and "euler3" in l:
            header_idx = i
            break
    if header_idx is None:
        raise ValueError("Could not find the CTF data header containing Phase/X/Y/Euler1/Euler2/Euler3.")
    data = pd.read_csv(StringIO("\n".join(lines[header_idx:])), sep=r"\s+|\t+", engine="python")
    lower = {str(c).lower(): c for c in data.columns}
    needed = ["phase", "x", "y", "euler1", "euler2", "euler3"]
    missing = [c for c in needed if c not in lower]
    if missing:
        raise ValueError(f"CTF file is missing expected columns: {', '.join(missing)}")
    return pd.DataFrame({
        "phase": pd.to_numeric(data[lower["phase"]], errors="coerce"),
        "x": pd.to_numeric(data[lower["x"]], errors="coerce"),
        "y": pd.to_numeric(data[lower["y"]], errors="coerce"),
        "phi1_deg": pd.to_numeric(data[lower["euler1"]], errors="coerce"),
        "Phi_deg": pd.to_numeric(data[lower["euler2"]], errors="coerce"),
        "phi2_deg": pd.to_numeric(data[lower["euler3"]], errors="coerce"),
    }).dropna()


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
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Explicit, transparent spatial/misorientation segmentation for modest EBSD maps.

    This is intentionally conservative and is not advertised as a replacement for a
    production EBSD segmentation package.  It is useful for small maps / teaching / an
    auditable first pass when raw point data are supplied.
    """
    required = {"x", "y", "phi1_deg", "Phi_deg", "phi2_deg"}
    if not required.issubset(points.columns):
        raise ValueError("Raw-point segmentation requires x, y, phi1_deg, Phi_deg and phi2_deg columns.")
    df = points.dropna(subset=list(required)).reset_index(drop=True).copy()
    if len(df) < 2:
        raise ValueError("At least two indexed EBSD points are required.")
    if len(df) > max_points:
        raise ValueError(
            f"This transparent in-app segmenter is capped at {max_points:,} points; your file has {len(df):,}. "
            "Segment the EBSD map in MTEX/OIM/AZtec/your normal EBSD workflow and upload the grain table instead."
        )

    xy = df[["x", "y"]].to_numpy(float)
    tree = cKDTree(xy)
    d, _ = tree.query(xy, k=2)
    step = float(np.median(d[:, 1][np.isfinite(d[:, 1]) & (d[:, 1] > 0)]))
    if not np.isfinite(step) or step <= 0:
        raise ValueError("Could not infer a positive EBSD map step size from x/y coordinates.")
    pairs = sorted(tree.query_pairs(r=1.30 * step))

    ori = [bunge_euler_to_matrix(float(a), float(b), float(c)) for a, b, c in df[["phi1_deg", "Phi_deg", "phi2_deg"]].to_numpy()]
    phase = df["phase"].to_numpy() if "phase" in df.columns else np.ones(len(df), int)

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
        if misorientation_deg(ori[i], ori[j], child_sym) <= misorientation_threshold_deg:
            union(i, j)

    root = np.array([find(i) for i in range(len(df))], int)
    counts = pd.Series(root).value_counts()
    keep_roots = set(int(r) for r, n in counts.items() if int(n) >= int(min_points))
    kept = np.array([r in keep_roots for r in root])
    if not np.any(kept):
        raise ValueError("No grains survived the minimum-point threshold. Lower the minimum grain size or inspect indexing/thresholds.")

    roots_sorted = sorted(keep_roots)
    rid_to_gid = {r: k + 1 for k, r in enumerate(roots_sorted)}
    point_gid = np.array([rid_to_gid.get(int(r), 0) for r in root], int)

    grain_rows = []
    for r in roots_sorted:
        idx = np.flatnonzero(root == r)
        mats = np.stack([ori[i] for i in idx])
        # Within a segmented grain orientations are already close by construction;
        # scipy's chordal mean is a transparent representative orientation.
        gmean = Rotation.from_matrix(mats).mean().as_matrix()
        p1, P, p2 = matrix_to_bunge_euler(gmean)
        row = {
            "grain_id": rid_to_gid[r],
            "x": float(df.loc[idx, "x"].mean()),
            "y": float(df.loc[idx, "y"].mean()),
            "phi1_deg": p1,
            "Phi_deg": P,
            "phi2_deg": p2,
            "point_count": int(len(idx)),
            "area_est_um2": float(len(idx) * step * step),
        }
        if "phase" in df.columns:
            row["phase"] = pd.Series(df.loc[idx, "phase"]).mode().iloc[0]
        grain_rows.append(row)
    grains = pd.DataFrame(grain_rows)

    edge_set: set[tuple[int, int]] = set()
    for i, j in pairs:
        gi, gj = int(point_gid[i]), int(point_gid[j])
        if gi > 0 and gj > 0 and gi != gj:
            edge_set.add((min(gi, gj), max(gi, gj)))
    adjacency = pd.DataFrame(sorted(edge_set), columns=["grain_id_1", "grain_id_2"])
    return grains, adjacency


def _parent_summary(result: ReconstructionResult, source_df: pd.DataFrame) -> pd.DataFrame:
    merged = source_df.merge(result.table, on="grain_id", how="right", suffixes=("", "_recon"))
    rows = []
    area_col = next((c for c in ["area_um2", "area_est_um2", "area"] if c in merged.columns), None)
    total_area = float(pd.to_numeric(merged[area_col], errors="coerce").fillna(0).sum()) if area_col else np.nan
    for pid, grp in merged.groupby("reconstructed_parent_id", sort=True):
        first = grp.iloc[0]
        fit = pd.to_numeric(grp["fit_deg"], errors="coerce")
        conf = pd.to_numeric(grp["confidence"], errors="coerce")
        area = float(pd.to_numeric(grp[area_col], errors="coerce").fillna(0).sum()) if area_col else np.nan
        mean_fit = float(fit.mean())
        mean_conf = float(conf.mean())
        rows.append({
            "parent_id": int(pid),
            "supporting_daughter_grains": int(len(grp)),
            "supporting_area_um2": area if area_col else np.nan,
            "area_fraction_pct": (100.0 * area / total_area) if area_col and total_area > 0 else np.nan,
            "mean_OR_fit_deg": mean_fit,
            "median_OR_fit_deg": float(fit.median()),
            "p95_OR_fit_deg": float(fit.quantile(0.95)),
            "max_OR_fit_deg": float(fit.max()),
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
        })
    return pd.DataFrame(rows)


def _assignment_table(result: ReconstructionResult, source_df: pd.DataFrame) -> pd.DataFrame:
    show = source_df.merge(result.table, on="grain_id", how="right")
    show = show.rename(columns={
        "reconstructed_parent_id": "parent_id",
        "variant_candidate_id": "candidate_variant_id",
        "fit_deg": "best_OR_fit_deg",
        "confidence": "support_score_0_to_1",
    })
    show["quality_flag"] = [
        _quality_label(float(f), float(c))
        for f, c in zip(show["best_OR_fit_deg"], show["support_score_0_to_1"])
    ]
    preferred = [
        "grain_id", "parent_id", "candidate_variant_id", "best_OR_fit_deg",
        "support_score_0_to_1", "quality_flag", "x", "y", "area_um2",
        "area_est_um2", "point_count", "phase", "parent_phi1_deg", "parent_Phi_deg",
        "parent_phi2_deg", "parent_qw", "parent_qx", "parent_qy", "parent_qz",
    ]
    cols = [c for c in preferred if c in show.columns] + [c for c in show.columns if c not in preferred]
    return show[cols]


def _comparison_table(results: dict[str, ReconstructionResult], runtimes: dict[str, float], truth: np.ndarray | None = None) -> pd.DataFrame:
    rows = []
    for name, r in results.items():
        fit = np.asarray(r.fit_deg, float)
        conf = np.asarray(r.confidence, float)
        row = {
            "method": name,
            "reconstructed_parents": int(r.diagnostics["n_reconstructed_parents"]),
            "mean_OR_fit_deg": float(np.mean(fit)),
            "median_OR_fit_deg": float(np.median(fit)),
            "p95_OR_fit_deg": float(np.quantile(fit, 0.95)),
            "mean_support_score": float(np.mean(conf)),
            "low_support_fraction_pct": float(100.0 * np.mean((conf < 0.35) | (fit > 5.0))),
            "runtime_s": float(runtimes.get(name, np.nan)),
        }
        if truth is not None:
            row["known_truth_clustering_accuracy"] = float(reconstruction_accuracy_against_labels(r, truth))
        rows.append(row)
    return pd.DataFrame(rows)


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
    st.subheader("2 · Define the parent → daughter crystallography")
    st.caption("OR = orientation relationship. Internally the engine stores a daughter→parent crystal-frame rotation matrix R_cp. You may use a literature OR, define one by parallel crystallographic plane/direction pairs, or supply the matrix directly.")
    mode = st.radio(
        "How will the orientation relationship (OR) be defined?",
        ["Literature preset", "Custom plane + direction parallelism", "Custom 3×3 OR matrix"],
        horizontal=True,
        key=f"{key_prefix}_or_mode",
    )
    if mode == "Literature preset":
        name = st.selectbox("Orientation relationship", list(presets), key=f"{key_prefix}_preset")
        p = presets[name]
        r = p.matrix_child_to_parent
        st.info(f"{p.name}: {p.plane_relation}; {p.direction_relation}. {p.note}")
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


def _prepare_dataset(child_sym: tuple[np.ndarray, ...]) -> PreparedDataset | None:
    st.subheader("1 · Provide the measured daughter microstructure")
    mode = st.radio(
        "Input route",
        [
            "Built-in validation dataset",
            "Manual grain table",
            "Upload pre-segmented grain table",
            "Upload raw EBSD .ang / .ctf and segment here",
        ],
        horizontal=False,
        help="Reconstruction algorithms operate on daughter-grain orientations plus grain adjacency. Raw EBSD pixels therefore require an explicit segmentation step before reconstruction.",
    )
    convention = "crystal_to_specimen"

    if mode == "Built-in validation dataset":
        c1, c2 = st.columns(2)
        n = c1.number_input("Daughter grains per known parent", 4, 30, 6, 1)
        noise = c2.number_input("Synthetic orientation noise — degrees", 0.0, 3.0, 0.35, 0.05)
        preset_name = "Kurdjumov–Sachs (FCC parent → BCC child)"
        df, edges, _ = synthetic_parent_reconstruction_demo(preset_name, int(n), float(noise), seed=11)
        st.success("Validation truth: this synthetic dataset contains exactly two parent grains. The true_parent_id column is retained only so reconstruction accuracy can be checked after the algorithms run.")
        return PreparedDataset("built-in validation", df, edges, convention, "synthetic two-parent validation dataset")

    if mode == "Manual grain table":
        df = _manual_grain_editor("manual_grains")
    else:
        files = st.file_uploader(
            "Daughter orientation data",
            type=["csv", "tsv", "txt", "ang", "ctf"],
            accept_multiple_files=False,
            help="For a pre-segmented grain table use one row per daughter grain. For .ang/.ctf the app reads EBSD points and requires the explicit segmentation controls below.",
        )
        if files is None:
            st.info("Upload a daughter-orientation file to continue.")
            return None
        try:
            df, provenance = _read_uploaded_orientation_file(files)
        except Exception as exc:
            st.error(f"Could not read the uploaded file: {exc}")
            return None
        if mode == "Upload raw EBSD .ang / .ctf and segment here" or "grain_id" not in df.columns:
            st.warning("Raw/pixel data detected. Reconstruction requires daughter grains, so the following segmentation is performed explicitly before any parent reconstruction. For publication-scale EBSD maps, a dedicated EBSD package remains preferable for segmentation.")
            s1, s2 = st.columns(2)
            seg_deg = s1.number_input("Daughter-grain segmentation misorientation threshold — degrees", 0.5, 15.0, 5.0, 0.5)
            min_pts = s2.number_input("Minimum indexed points per retained daughter grain", 2, 200, 5, 1)
            if st.button("Segment raw EBSD points into daughter grains", type="primary"):
                try:
                    grains, adj = _segment_points_to_grains(df, child_sym, float(seg_deg), int(min_pts))
                    st.session_state["seg_grains"] = grains
                    st.session_state["seg_adj"] = adj
                    st.session_state["seg_name"] = files.name
                except Exception as exc:
                    st.error(f"Segmentation failed: {exc}")
            if st.session_state.get("seg_name") != files.name:
                return None
            grains = st.session_state["seg_grains"]
            adj = st.session_state["seg_adj"]
            st.write(f"Segmented daughter grains: **{len(grains)}**; measured inter-grain adjacency pairs: **{len(adj)}**")
            st.dataframe(grains.head(30), use_container_width=True, hide_index=True)
            return PreparedDataset(files.name, grains, adj, convention, f"{provenance}; in-app explicit spatial/misorientation segmentation")

    df = _coerce_grain_columns(df)
    c1, c2 = st.columns(2)
    with c1:
        conv = st.radio("Orientation matrix direction / Euler convention", ["Crystal → specimen", "Specimen → crystal"], horizontal=True)
        convention = "crystal_to_specimen" if conv == "Crystal → specimen" else "specimen_to_crystal"
    with c2:
        st.caption("Internal convention is crystal → specimen. Select specimen → crystal only when your export documentation explicitly says the orientation maps specimen coordinates back into crystal coordinates.")

    # Check grain table before asking for adjacency.
    try:
        grain_ids, _ = orientations_from_dataframe(df, convention=convention)
    except Exception as exc:
        st.error(f"Grain-orientation table is not ready: {exc}")
        st.caption("Required: unique grain_id plus either phi1_deg/Phi_deg/phi2_deg (degrees) or qw/qx/qy/qz. Optional but useful: x, y, area_um2, phase, true_parent_id.")
        return None

    st.markdown("**Grain adjacency / physical neighborhood**")
    adj_mode_options = ["Upload measured adjacency", "Enter adjacency manually"]
    if {"x", "y"}.issubset(df.columns):
        adj_mode_options.append("Approximate from grain centroids (k-NN)")
    adj_mode = st.radio("How are neighboring daughter grains defined?", adj_mode_options, horizontal=True)
    if adj_mode == "Upload measured adjacency":
        af = st.file_uploader("Adjacency CSV", type=["csv", "tsv", "txt"], help="Columns: grain_id_1, grain_id_2. Each row is one measured grain-boundary adjacency pair.")
        if af is None:
            st.info("Upload the adjacency table to continue, or choose manual/k-NN adjacency.")
            return None
        adj = pd.read_csv(af, sep=None, engine="python")
    elif adj_mode == "Enter adjacency manually":
        adj = _manual_adjacency_editor([int(x) for x in grain_ids], "manual_adj")
    else:
        k = st.number_input("Centroid neighbors per grain — k", 1, 12, 4, 1)
        st.warning("k-NN adjacency is an approximation based only on grain centroids. It must not be described as measured grain-boundary topology in a paper.")
        try:
            idx_edges = approximate_knn_edges(df, int(k))
            adj = pd.DataFrame({
                "grain_id_1": [int(grain_ids[i]) for i, _ in idx_edges],
                "grain_id_2": [int(grain_ids[j]) for _, j in idx_edges],
            })
        except Exception as exc:
            st.error(str(exc))
            return None

    try:
        edges_from_dataframe(adj, grain_ids)
    except Exception as exc:
        st.error(f"Adjacency table is not valid: {exc}")
        return None
    return PreparedDataset("manual" if mode == "Manual grain table" else files.name, df, adj, convention, mode)


def _render_method_controls(selected: list[str]) -> dict[str, float]:
    st.subheader("4 · Choose reconstruction algorithms and reportable thresholds")
    st.caption("Only controls used by the selected algorithms are shown. Every angular tolerance below is in degrees and should be reported in a paper or supplementary methods.")
    controls: dict[str, float] = {}
    if any(m in selected for m in ["Neighbor voting", "Grain graph + Markov clustering", "Variant graph probability propagation"]):
        controls["sigma_deg"] = float(st.number_input("Candidate/edge angular width σ — degrees", 0.2, 15.0, 2.5, 0.1, help="Controls how rapidly support decays as parent-candidate angular disagreement increases."))
    else:
        controls["sigma_deg"] = 2.5
    if any(m in selected for m in ["Grain graph + Markov clustering", "Variant graph probability propagation"]):
        controls["inflation"] = float(st.number_input("Markov graph inflation parameter — dimensionless", 1.01, 4.0, 1.6, 0.05, help="Higher inflation generally favors more strongly separated / smaller graph clusters."))
    else:
        controls["inflation"] = 1.6
    if "Neighbor voting" in selected:
        controls["merge_deg"] = float(st.number_input("Neighbor-voting parent merge tolerance — degrees", 0.5, 20.0, 8.0, 0.5))
    else:
        controls["merge_deg"] = 8.0
    if "Nucleation + growth" in selected:
        a, b = st.columns(2)
        controls["nucleation_deg"] = float(a.number_input("Nucleation tolerance — degrees", 0.2, 10.0, 3.0, 0.2, help="Strict angular criterion used to create high-confidence parent nuclei."))
        controls["growth_deg"] = float(b.number_input("Growth tolerance — degrees", 0.5, 20.0, 8.0, 0.5, help="Looser angular criterion used to add neighboring daughter grains to established nuclei."))
    else:
        controls["nucleation_deg"] = 3.0
        controls["growth_deg"] = 8.0
    if "Operator / groupoid consistency" in selected:
        a, b = st.columns(2)
        controls["operator_tol_deg"] = float(a.number_input("Transformation-operator tolerance — degrees", 0.2, 20.0, 8.0, 0.5))
        controls["parent_consistency_deg"] = float(b.number_input("Common-parent consistency tolerance — degrees", 0.2, 20.0, 5.0, 0.5))
    else:
        controls["operator_tol_deg"] = 8.0
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


def _render_result(name: str, result: ReconstructionResult, source_df: pd.DataFrame, comparison: pd.DataFrame) -> None:
    st.subheader(f"Detailed result · {name}")
    parent_summary = _parent_summary(result, source_df)
    assignments = _assignment_table(result, source_df)

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Reconstructed parent grains", int(result.diagnostics["n_reconstructed_parents"]))
    m2.metric("Mean OR fit", f"{float(np.mean(result.fit_deg)):.3f}°")
    m3.metric("95th-percentile OR fit", f"{float(np.quantile(result.fit_deg, 0.95)):.3f}°")
    m4.metric("Mean heuristic support", f"{float(np.mean(result.confidence)):.3f}")
    st.caption("OR fit = angular mismatch between a daughter grain's selected parent candidate and the representative reconstructed parent orientation. Lower is better. Heuristic support = 0–1 algorithmic support from fit/separation; it is NOT a calibrated probability.")

    t1, t2, t3 = st.tabs(["Parent summary", "Daughter assignments", "Maps & diagnostics"])
    with t1:
        st.dataframe(parent_summary, use_container_width=True, hide_index=True)
        st.caption("Parent summary: one row per reconstructed parent. supporting_daughter_grains = evidence count; mean/median/p95/max_OR_fit_deg = angular residual statistics in degrees; mean/minimum_support_score = heuristic 0–1 support; weak_or_ambiguous_children counts daughters with fit >5° or support <0.35; Euler/quaternion columns describe the reconstructed parent orientation.")
        st.download_button("Download parent summary CSV", parent_summary.to_csv(index=False), f"{name.replace(' ','_')}_parent_summary.csv", "text/csv")
    with t2:
        qfilter = st.multiselect("Show quality bands", ["Very strong", "Strong", "Review", "Weak / ambiguous"], default=["Very strong", "Strong", "Review", "Weak / ambiguous"])
        view = assignments[assignments["quality_flag"].isin(qfilter)]
        st.dataframe(view, use_container_width=True, hide_index=True)
        st.caption("Daughter assignment: parent_id = reconstructed prior-parent cluster; candidate_variant_id = which crystallographically allowed parent candidate was selected for this daughter grain; best_OR_fit_deg = angular fit in degrees; support_score_0_to_1 = heuristic algorithmic support, not probability. Variant IDs are internal candidate IDs and are not automatically packet/Bain/KS labels unless a dedicated mapping is supplied.")
        st.download_button("Download daughter assignment CSV", assignments.to_csv(index=False), f"{name.replace(' ','_')}_daughter_assignments.csv", "text/csv")
    with t3:
        if {"x", "y"}.issubset(assignments.columns):
            fig = px.scatter(
                assignments,
                x="x", y="y", color="parent_id", symbol="candidate_variant_id",
                hover_data=["grain_id", "best_OR_fit_deg", "support_score_0_to_1", "quality_flag"],
                title="Reconstructed parent map from daughter-grain centroids",
                labels={"x": "Map x coordinate", "y": "Map y coordinate", "parent_id": "Reconstructed parent ID", "candidate_variant_id": "Candidate variant ID"},
            )
            fig.update_yaxes(scaleanchor="x", scaleratio=1)
            st.plotly_chart(fig, use_container_width=True)
        c1, c2 = st.columns(2)
        with c1:
            figf = px.histogram(assignments, x="best_OR_fit_deg", nbins=30, title="Distribution of OR angular residuals", labels={"best_OR_fit_deg": "Best OR fit (degrees)"})
            st.plotly_chart(figf, use_container_width=True)
        with c2:
            figc = px.histogram(assignments, x="support_score_0_to_1", nbins=20, title="Distribution of heuristic support", labels={"support_score_0_to_1": "Support score (0–1; not probability)"})
            st.plotly_chart(figc, use_container_width=True)

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
        "**Purpose:** infer prior-parent orientations/grains from a measured daughter microstructure, or generate all crystallographically possible daughter variants from known parent orientations. "
        "The page separates data preparation, OR definition, OR refinement, reconstruction method and result confidence so every published number has a clear origin."
    )
    st.info("Core reconstruction input = daughter-grain orientations + parent/daughter symmetries + orientation relationship (OR) + grain adjacency. Lattice parameters are not required for orientation-only reconstruction unless a separate metric/compatibility calculation needs them.")

    tab_recon, tab_forward, tab_batch, tab_dictionary = st.tabs([
        "Daughter → parent", "Parent → daughter variants", "Compare many datasets", "Input/output dictionary"
    ])

    with tab_recon:
        # OR is needed before raw EBSD segmentation because segmentation uses daughter symmetry.
        or_name, preset, initial_or, parent_sym, child_sym = _render_or_builder("main")
        ds = _prepare_dataset(child_sym)
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
            comp = _comparison_table(results, runtimes, truth=truth)
            st.session_state["wb_results"] = results
            st.session_state["wb_comparison"] = comp
            st.session_state["wb_signature"] = sig
            st.session_state["wb_source"] = ds.grains.copy()

        results = st.session_state.get("wb_results")
        comp = st.session_state.get("wb_comparison")
        if results and st.session_state.get("wb_signature") != sig:
            st.warning("Inputs, OR, selected methods or thresholds changed after the last run. Run reconstruction again before interpreting the old result.")
            results = None
            comp = None
        if results and comp is not None:
            st.subheader("5 · Cross-method evidence before choosing what to report")
            st.dataframe(comp, use_container_width=True, hide_index=True)
            st.caption("Do not choose a method only because it gives the smallest parent count or smallest fit. Compare OR residuals, low-support fraction, physical map coherence, known truth if available, and agreement with independent reconstruction methods.")
            if len(results) > 1:
                st.markdown("**Pairwise clustering agreement — Adjusted Rand Index (ARI)**")
                ari = _agreement_matrix(results)
                st.dataframe(ari.style.format("{:.3f}"), use_container_width=True)
                st.caption("ARI = 1 means two methods partition daughter grains into identical parent clusters up to arbitrary cluster labels; ARI near 0 means agreement is around chance. It compares clustering only, not reconstructed orientation accuracy.")
            inspect = st.selectbox("Inspect one method in detail", list(results))
            _render_result(inspect, results[inspect], ds.grains, comp)

            with st.expander("Optional retained/reference parent validation"):
                st.caption("If independent parent orientations are available (retained parent phase, serial experiment, simulation truth), upload them here. The app reports the nearest crystallographically equivalent reference-parent misorientation. Without a known spatial/ID correspondence this is a nearest-reference check, not proof of identity.")
                rf = st.file_uploader("Reference parent orientation CSV", type=["csv"], key="ref_parent")
                if rf is not None:
                    rdf = _coerce_grain_columns(pd.read_csv(rf))
                    if "grain_id" not in rdf.columns:
                        rdf.insert(0, "grain_id", np.arange(1, len(rdf) + 1))
                    try:
                        rv = _reference_parent_validation(results[inspect], rdf, parent_sym, "crystal_to_specimen")
                        st.dataframe(rv, use_container_width=True, hide_index=True)
                    except Exception as exc:
                        st.error(str(exc))

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
            ["support_score_0_to_1", "Heuristic support built from absolute fit and separation from alternatives", "0–1", "Not a probability"],
            ["parent_phi1_deg, parent_Phi_deg, parent_phi2_deg", "Reconstructed parent orientation in Bunge Euler coordinates", "degrees", "Parent-level output"],
            ["parent_qw, parent_qx, parent_qy, parent_qz", "Same reconstructed parent orientation as a unit quaternion", "dimensionless", "Parent-level output"],
            ["ARI", "Adjusted Rand Index comparing the parent clustering produced by two methods", "dimensionless", "1 = identical clustering up to labels"],
        ], columns=["field / symbol", "meaning", "unit", "how to use it"])
        st.dataframe(dictionary, use_container_width=True, hide_index=True)
        st.markdown(
            "**Minimum defensible Methods report:** daughter phase and symmetry; parent phase and symmetry; OR definition and whether refined; orientation convention; grain segmentation method/threshold; adjacency source; selected reconstruction algorithm(s); every angular/graph threshold; number of daughter grains; reconstructed parent count; OR-fit statistics; confidence/support definition; and any retained-parent or known-truth validation."
        )
        st.info("These implementations are transparent members of established reconstruction families, not binary-identical reproductions of MTEX, ARPGE, OIM Analysis, AZtecCrystal or other vendor/research packages.")
