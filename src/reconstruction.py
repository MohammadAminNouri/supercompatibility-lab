from __future__ import annotations

"""Parent/daughter orientation reconstruction utilities.

The module is intentionally independent from EBSD file formats. It works on grain-level
orientation matrices and an adjacency list. Orientations are represented by proper 3x3
rotation matrices mapping crystal coordinates into the specimen frame.

The orientation relationship (OR) is represented as a child->parent crystal-frame rotation
R_cp. Thus a reference child orientation generated from a parent orientation Gp is
Gc = Gp @ R_cp. Parent symmetry operations generate transformation variants, while
child symmetry operations generate the candidate parent orientations compatible with a
measured child orientation.
"""

from dataclasses import dataclass
from itertools import permutations, product
from typing import Iterable, Literal

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation

Array = np.ndarray


@dataclass(frozen=True)
class ORPreset:
    name: str
    parent_phase: str
    child_phase: str
    parent_symmetry: str
    child_symmetry: str
    matrix_child_to_parent: Array
    plane_relation: str
    direction_relation: str
    note: str = ""


@dataclass(frozen=True)
class ReconstructionResult:
    method: str
    table: pd.DataFrame
    parent_orientations: dict[int, Array]
    parent_ids: Array
    variant_ids: Array
    fit_deg: Array
    confidence: Array
    diagnostics: dict[str, float | int | str]


# ---------- rotation and orientation helpers ----------

def _project_so3(m: Array) -> Array:
    u, _, vt = np.linalg.svd(np.asarray(m, float))
    r = u @ vt
    if np.linalg.det(r) < 0:
        u[:, -1] *= -1
        r = u @ vt
    return r


def rotation_angle_deg(r: Array) -> float:
    r = _project_so3(r)
    c = np.clip((np.trace(r) - 1.0) / 2.0, -1.0, 1.0)
    return float(np.degrees(np.arccos(c)))


def bunge_euler_to_matrix(phi1_deg: float, Phi_deg: float, phi2_deg: float) -> Array:
    """Bunge ZXZ Euler angles in degrees -> crystal-to-specimen orientation matrix."""
    # scipy intrinsic ZXZ is consistent with the common Bunge active rotation matrix
    # after using the crystal->specimen convention adopted here.
    return Rotation.from_euler("ZXZ", [phi1_deg, Phi_deg, phi2_deg], degrees=True).as_matrix()


def matrix_to_bunge_euler(g: Array) -> tuple[float, float, float]:
    ang = Rotation.from_matrix(_project_so3(g)).as_euler("ZXZ", degrees=True)
    ang = np.mod(ang, 360.0)
    # conventionally Phi is restricted to [0,180]
    if ang[1] > 180:
        ang[1] = 360.0 - ang[1]
    return float(ang[0]), float(ang[1]), float(ang[2])


def quaternion_wxyz(g: Array) -> tuple[float, float, float, float]:
    qx, qy, qz, qw = Rotation.from_matrix(_project_so3(g)).as_quat()
    if qw < 0:
        qw, qx, qy, qz = -qw, -qx, -qy, -qz
    return float(qw), float(qx), float(qy), float(qz)


def matrix_from_quaternion_wxyz(qw: float, qx: float, qy: float, qz: float) -> Array:
    q = np.array([qx, qy, qz, qw], float)
    n = np.linalg.norm(q)
    if n <= 1e-15:
        raise ValueError("Quaternion cannot be zero.")
    return Rotation.from_quat(q / n).as_matrix()


# ---------- proper rotation symmetry groups ----------

def cubic_rotations() -> tuple[Array, ...]:
    mats: dict[tuple[int, ...], Array] = {}
    for perm in permutations(range(3)):
        p = np.zeros((3, 3))
        p[np.arange(3), perm] = 1.0
        for signs in product((-1.0, 1.0), repeat=3):
            g = np.diag(signs) @ p
            if np.linalg.det(g) > 0.5:
                key = tuple(np.rint(g).astype(int).ravel())
                mats[key] = g
    return tuple(mats[k] for k in sorted(mats))


def orthorhombic_rotations() -> tuple[Array, ...]:
    return (
        np.eye(3),
        np.diag([1.0, -1.0, -1.0]),
        np.diag([-1.0, 1.0, -1.0]),
        np.diag([-1.0, -1.0, 1.0]),
    )


def monoclinic_rotations() -> tuple[Array, ...]:
    return (np.eye(3), np.diag([-1.0, 1.0, -1.0]))


def hexagonal_rotations() -> tuple[Array, ...]:
    mats: list[Array] = []
    for k in range(6):
        rz = Rotation.from_rotvec(np.array([0.0, 0.0, np.deg2rad(60.0 * k)])).as_matrix()
        mats.append(rz)
    rx = Rotation.from_rotvec(np.array([np.pi, 0.0, 0.0])).as_matrix()
    for k in range(6):
        rz = Rotation.from_rotvec(np.array([0.0, 0.0, np.deg2rad(60.0 * k)])).as_matrix()
        mats.append(rz @ rx)
    # unique / projected for numerical cleanliness
    out: list[Array] = []
    for m in mats:
        m = _project_so3(m)
        if not any(np.allclose(m, q, atol=1e-10) for q in out):
            out.append(m)
    return tuple(out)


def symmetry_group(name: str) -> tuple[Array, ...]:
    key = name.strip().lower()
    if key in {"cubic", "m-3m", "432", "fcc", "bcc", "b2"}:
        return cubic_rotations()
    if key in {"hexagonal", "6/mmm", "622", "hcp"}:
        return hexagonal_rotations()
    if key in {"orthorhombic", "mmm", "222", "b19"}:
        return orthorhombic_rotations()
    if key in {"monoclinic", "2/m", "2", "b19'", "b19p"}:
        return monoclinic_rotations()
    if key in {"triclinic", "1"}:
        return (np.eye(3),)
    raise ValueError(f"Unsupported proper-rotation symmetry group: {name}")


def misorientation_deg(g1: Array, g2: Array, symmetry: Iterable[Array]) -> float:
    """Minimum disorientation angle modulo right-acting crystal symmetry.

    This hot-path implementation is vectorized over the symmetry group.
    """
    g1 = np.asarray(g1, float)
    g2 = np.asarray(g2, float)
    syms = np.asarray(tuple(symmetry), float)
    # trace(g1 S g2^T) = trace((g2^T g1) S)
    a = g2.T @ g1
    traces = np.einsum("ij,kji->k", a, syms)
    cosines = np.clip((traces - 1.0) / 2.0, -1.0, 1.0)
    return float(np.degrees(np.arccos(np.max(cosines))))


# ---------- orientation relationships ----------

def _orthobasis_from_plane_direction(normal: Iterable[float], direction: Iterable[float]) -> Array:
    n = np.asarray(normal, float)
    d = np.asarray(direction, float)
    if np.linalg.norm(n) <= 1e-14 or np.linalg.norm(d) <= 1e-14:
        raise ValueError("Plane normal and direction must be non-zero.")
    n = n / np.linalg.norm(n)
    # remove accidental component normal to the plane
    d = d - n * float(n @ d)
    if np.linalg.norm(d) <= 1e-12:
        raise ValueError("Direction must lie in the specified plane.")
    d = d / np.linalg.norm(d)
    t = np.cross(n, d)
    t /= np.linalg.norm(t)
    # columns: in-plane direction, second in-plane direction, plane normal
    return np.column_stack([d, t, n])


def rotation_from_parallelisms(
    parent_plane_normal: Iterable[float],
    parent_direction: Iterable[float],
    child_plane_normal: Iterable[float],
    child_direction: Iterable[float],
) -> Array:
    """Return child->parent OR satisfying one plane and one in-plane direction parallelism."""
    bp = _orthobasis_from_plane_direction(parent_plane_normal, parent_direction)
    bc = _orthobasis_from_plane_direction(child_plane_normal, child_direction)
    return _project_so3(bp @ bc.T)


def orientation_relationship_presets() -> dict[str, ORPreset]:
    # Representative members of the standard OR families. Symmetry generates the
    # remaining equivalent variants.
    ks = rotation_from_parallelisms(
        [1, 1, 1], [1, -1, 0],
        [1, 1, 0], [1, -1, 1],
    )
    nw = rotation_from_parallelisms(
        [1, 1, 1], [0, 1, -1],
        [1, 1, 0], [0, 0, 1],
    )
    bain = rotation_from_parallelisms(
        [0, 1, 0], [0, 0, 1],
        [0, 1, 0], [1, 0, 1],
    )
    pitsch = rotation_from_parallelisms(
        [1, 0, 0], [0, 1, 1],
        [1, 1, 0], [1, -1, 1],
    )
    burgers = rotation_from_parallelisms(
        [1, 1, 0], [1, -1, 1],
        [0, 0, 1], [1, 0, 0],  # child HCP Cartesian: basal a direction and c-axis normal
    )
    return {
        "Kurdjumov–Sachs (FCC parent → BCC child)": ORPreset(
            "Kurdjumov–Sachs", "FCC", "BCC", "cubic", "cubic", ks,
            "{111}parent ∥ {110}child", "<110>parent ∥ <111>child",
            "24 orientation variants for ideal cubic symmetries.",
        ),
        "Nishiyama–Wassermann (FCC parent → BCC child)": ORPreset(
            "Nishiyama–Wassermann", "FCC", "BCC", "cubic", "cubic", nw,
            "{111}parent ∥ {110}child", "<011>parent ∥ <001>child",
            "12 orientation variants for ideal cubic symmetries.",
        ),
        "Bain (FCC parent → BCC child)": ORPreset(
            "Bain", "FCC", "BCC", "cubic", "cubic", bain,
            "{010}parent ∥ {010}child", "<001>parent ∥ <101>child",
            "Classical rational FCC/BCC orientation relationship.",
        ),
        "Pitsch (FCC parent → BCC child)": ORPreset(
            "Pitsch", "FCC", "BCC", "cubic", "cubic", pitsch,
            "{100}parent ∥ {110}child", "<110>parent ∥ <111>child",
            "Classical rational FCC/BCC orientation relationship.",
        ),
        "Burgers (BCC parent → HCP child)": ORPreset(
            "Burgers", "BCC β", "HCP α", "cubic", "hexagonal", burgers,
            "{110}parent ∥ (0001)child", "<111>parent ∥ <11-20>child",
            "Standard β→α relation used in titanium alloys.",
        ),
    }


# ---------- variants and candidate parents ----------

def unique_child_variants(
    parent_orientation: Array,
    r_child_to_parent: Array,
    parent_sym: Iterable[Array],
    child_sym: Iterable[Array],
    tol_deg: float = 1e-5,
) -> list[Array]:
    """Transformation variants Gc = Gp S_parent R_cp, unique modulo child symmetry."""
    out: list[Array] = []
    for sp in parent_sym:
        g = _project_so3(parent_orientation @ sp @ r_child_to_parent)
        if not any(misorientation_deg(g, h, child_sym) <= tol_deg for h in out):
            out.append(g)
    return out


def parent_candidates(
    child_orientation: Array,
    r_child_to_parent: Array,
    child_sym: Iterable[Array],
    parent_sym: Iterable[Array],
    tol_deg: float = 1e-5,
) -> list[tuple[Array, int]]:
    """All crystallographically distinct parent candidates for one child grain.

    Measured child orientations are equivalent under right multiplication by child
    symmetry. For every child symmetry S_c we back-transform
        Gp = Gc S_c R_cp^T
    and deduplicate modulo parent symmetry.
    """
    out: list[tuple[Array, int]] = []
    for j, sc in enumerate(child_sym, 1):
        gp = _project_so3(child_orientation @ sc @ r_child_to_parent.T)
        if not any(misorientation_deg(gp, old, parent_sym) <= tol_deg for old, _ in out):
            out.append((gp, j))
    return out


def candidate_sets(
    child_orientations: list[Array],
    r_child_to_parent: Array,
    child_sym: Iterable[Array],
    parent_sym: Iterable[Array],
) -> list[list[tuple[Array, int]]]:
    return [parent_candidates(g, r_child_to_parent, child_sym, parent_sym) for g in child_orientations]


def _best_pair_misfit(c1: list[tuple[Array, int]], c2: list[tuple[Array, int]], parent_sym: Iterable[Array]) -> float:
    return min(misorientation_deg(a, b, parent_sym) for a, _ in c1 for b, _ in c2)


def _gaussian_weight(angle_deg: float, threshold_deg: float, tolerance_deg: float) -> float:
    # Smooth probability-like weight: 0.5 around threshold when tolerance ~= threshold/1.177.
    sigma = max(float(tolerance_deg), 1e-6)
    return float(np.exp(-0.5 * (max(angle_deg, 0.0) / sigma) ** 2))


def _normalize_edges(edges: Iterable[tuple[int, int]], n: int) -> list[tuple[int, int]]:
    out: set[tuple[int, int]] = set()
    for a, b in edges:
        i, j = int(a), int(b)
        if i == j:
            continue
        if not (0 <= i < n and 0 <= j < n):
            raise ValueError("Adjacency contains an out-of-range grain index.")
        out.add((min(i, j), max(i, j)))
    return sorted(out)


class _UF:
    def __init__(self, n: int):
        self.p = list(range(n))
    def find(self, x: int) -> int:
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x
    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[rb] = ra


def _clusters_from_selected(selected: list[Array], edges: list[tuple[int, int]], parent_sym: Iterable[Array], threshold_deg: float) -> Array:
    uf = _UF(len(selected))
    for i, j in edges:
        if misorientation_deg(selected[i], selected[j], parent_sym) <= threshold_deg:
            uf.union(i, j)
    roots: dict[int, int] = {}
    ids = np.zeros(len(selected), dtype=int)
    for i in range(len(selected)):
        r = uf.find(i)
        if r not in roots:
            roots[r] = len(roots) + 1
        ids[i] = roots[r]
    return ids


def _medoid_orientation(orientations: list[Array], parent_sym: Iterable[Array], weights: Array | None = None) -> Array:
    if len(orientations) == 1:
        return orientations[0]
    w = np.ones(len(orientations)) if weights is None else np.asarray(weights, float)
    scores = []
    for i, g in enumerate(orientations):
        scores.append(sum(w[j] * misorientation_deg(g, h, parent_sym) ** 2 for j, h in enumerate(orientations)))
    return orientations[int(np.argmin(scores))]


def _assign_to_parent(cands: list[tuple[Array, int]], gp: Array, parent_sym: Iterable[Array]) -> tuple[Array, int, float, float]:
    vals = sorted((misorientation_deg(c, gp, parent_sym), idx, c) for c, idx in cands)
    best_fit, best_var, best_c = vals[0]
    second = vals[1][0] if len(vals) > 1 else 180.0
    # confidence combines absolute fit and separation to the next candidate
    abs_conf = np.exp(-0.5 * (best_fit / 3.0) ** 2)
    sep_conf = 1.0 - np.exp(-max(second - best_fit, 0.0) / 3.0)
    conf = float(np.clip(0.55 * abs_conf + 0.45 * sep_conf, 0.0, 1.0))
    return best_c, int(best_var), float(best_fit), conf


def _finalize(
    method: str,
    grain_ids: Array,
    child_orientations: list[Array],
    cands: list[list[tuple[Array, int]]],
    parent_ids: Array,
    parent_sym: Iterable[Array],
    extra: dict[str, float | int | str] | None = None,
) -> ReconstructionResult:
    parent_orientations: dict[int, Array] = {}
    variant_ids = np.zeros(len(child_orientations), dtype=int)
    fit = np.zeros(len(child_orientations), dtype=float)
    conf = np.zeros(len(child_orientations), dtype=float)

    for pid in sorted(set(parent_ids.tolist())):
        idx = np.flatnonzero(parent_ids == pid)
        # first medoid from all locally selected best-consensus candidates
        pool: list[Array] = []
        for i in idx:
            # choose candidate most mutually consistent with candidates of other grains in the cluster
            local_scores = []
            for cand, _ in cands[i]:
                sc = 0.0
                for j in idx:
                    if j == i:
                        continue
                    sc += min(misorientation_deg(cand, cj, parent_sym) ** 2 for cj, _ in cands[j])
                local_scores.append(sc)
            pool.append(cands[i][int(np.argmin(local_scores))][0])
        gp = _medoid_orientation(pool, parent_sym)
        # two assignment/refit passes
        for _ in range(2):
            aligned: list[Array] = []
            for i in idx:
                best, _, _, _ = _assign_to_parent(cands[i], gp, parent_sym)
                aligned.append(best)
            gp = _medoid_orientation(aligned, parent_sym)
        parent_orientations[int(pid)] = gp
        for i in idx:
            _, vid, fi, ci = _assign_to_parent(cands[i], gp, parent_sym)
            variant_ids[i] = vid
            fit[i] = fi
            conf[i] = ci

    rows = []
    for i, g in enumerate(child_orientations):
        p = parent_orientations[int(parent_ids[i])]
        ph1, PH, ph2 = matrix_to_bunge_euler(p)
        qw, qx, qy, qz = quaternion_wxyz(p)
        candidate_fits = sorted(float(misorientation_deg(c, p, parent_sym)) for c, _ in cands[i])
        second_fit = candidate_fits[1] if len(candidate_fits) > 1 else 180.0
        candidate_gap = max(second_fit - float(fit[i]), 0.0)
        absolute_fit_support = float(np.exp(-0.5 * (float(fit[i]) / 3.0) ** 2))
        separation_support = float(1.0 - np.exp(-candidate_gap / 3.0))
        rows.append({
            "grain_id": int(grain_ids[i]),
            "reconstructed_parent_id": int(parent_ids[i]),
            "variant_candidate_id": int(variant_ids[i]),
            "fit_deg": float(fit[i]),
            "second_best_candidate_fit_deg": float(second_fit),
            "candidate_separation_deg": float(candidate_gap),
            "candidate_count": int(len(candidate_fits)),
            "absolute_fit_support": absolute_fit_support,
            "separation_support": separation_support,
            "confidence": float(conf[i]),
            "parent_phi1_deg": ph1,
            "parent_Phi_deg": PH,
            "parent_phi2_deg": ph2,
            "parent_qw": qw,
            "parent_qx": qx,
            "parent_qy": qy,
            "parent_qz": qz,
        })
    diag: dict[str, float | int | str] = {
        "n_grains": len(child_orientations),
        "n_reconstructed_parents": len(parent_orientations),
        "mean_fit_deg": float(np.mean(fit)),
        "median_fit_deg": float(np.median(fit)),
        "mean_confidence": float(np.mean(conf)),
    }
    if extra:
        diag.update(extra)
    return ReconstructionResult(method, pd.DataFrame(rows), parent_orientations, np.asarray(parent_ids, int), variant_ids, fit, conf, diag)


# ---------- reconstruction algorithms ----------

def neighbor_voting_reconstruction(
    grain_ids: Array,
    child_orientations: list[Array],
    edges: Iterable[tuple[int, int]],
    r_child_to_parent: Array,
    child_sym: Iterable[Array],
    parent_sym: Iterable[Array],
    sigma_deg: float = 2.5,
    merge_deg: float = 5.0,
) -> ReconstructionResult:
    n = len(child_orientations)
    ed = _normalize_edges(edges, n)
    cands = candidate_sets(child_orientations, r_child_to_parent, child_sym, parent_sym)
    neigh: list[list[int]] = [[] for _ in range(n)]
    for i, j in ed:
        neigh[i].append(j); neigh[j].append(i)
    selected: list[Array] = []
    vote_margin = []
    for i in range(n):
        scores = []
        for cand, _ in cands[i]:
            score = 0.0
            for j in neigh[i]:
                a = min(misorientation_deg(cand, cj, parent_sym) for cj, _ in cands[j])
                score += _gaussian_weight(a, sigma_deg, sigma_deg)
            scores.append(score)
        order = np.argsort(scores)[::-1]
        selected.append(cands[i][int(order[0])][0])
        second = scores[int(order[1])] if len(order) > 1 else 0.0
        vote_margin.append(float(scores[int(order[0])] - second))
    pids = _clusters_from_selected(selected, ed, parent_sym, merge_deg)
    return _finalize("Neighbor voting", grain_ids, child_orientations, cands, pids, parent_sym, {
        "mean_vote_margin": float(np.mean(vote_margin)), "n_edges": len(ed)
    })


def _mcl(matrix: Array, inflation: float = 1.6, expansion: int = 2, max_iter: int = 100, tol: float = 1e-7, prune: float = 1e-5) -> Array:
    m = np.asarray(matrix, float).copy()
    n = m.shape[0]
    m = np.maximum(m, 0.0)
    np.fill_diagonal(m, np.maximum(np.diag(m), 1.0))
    colsum = m.sum(axis=0)
    colsum[colsum == 0] = 1.0
    m /= colsum
    for _ in range(max_iter):
        old = m.copy()
        # expansion
        x = m.copy()
        for _k in range(expansion - 1):
            x = x @ m
        # inflation
        x = np.power(np.maximum(x, 0.0), inflation)
        x[x < prune] = 0.0
        colsum = x.sum(axis=0)
        colsum[colsum == 0] = 1.0
        m = x / colsum
        if np.max(np.abs(m - old)) < tol:
            break
    return m


def _clusters_from_mcl(m: Array, threshold: float = 1e-3) -> Array:
    # connect nodes with persistent mutual probability; then connected components
    a = (m + m.T) > threshold
    np.fill_diagonal(a, True)
    n = len(a)
    uf = _UF(n)
    for i in range(n):
        for j in range(i + 1, n):
            if a[i, j]:
                uf.union(i, j)
    roots: dict[int, int] = {}
    ids = np.zeros(n, dtype=int)
    for i in range(n):
        r = uf.find(i)
        if r not in roots:
            roots[r] = len(roots) + 1
        ids[i] = roots[r]
    return ids


def grain_graph_reconstruction(
    grain_ids: Array,
    child_orientations: list[Array],
    edges: Iterable[tuple[int, int]],
    r_child_to_parent: Array,
    child_sym: Iterable[Array],
    parent_sym: Iterable[Array],
    sigma_deg: float = 2.5,
    inflation: float = 1.6,
) -> ReconstructionResult:
    n = len(child_orientations)
    ed = _normalize_edges(edges, n)
    cands = candidate_sets(child_orientations, r_child_to_parent, child_sym, parent_sym)
    w = np.zeros((n, n), float)
    np.fill_diagonal(w, 1.0)
    edge_misfits = []
    for i, j in ed:
        ang = _best_pair_misfit(cands[i], cands[j], parent_sym)
        edge_misfits.append(ang)
        ww = _gaussian_weight(ang, sigma_deg, sigma_deg)
        w[i, j] = w[j, i] = ww
    m = _mcl(w, inflation=inflation)
    pids = _clusters_from_mcl(m)
    return _finalize("Grain graph + Markov clustering", grain_ids, child_orientations, cands, pids, parent_sym, {
        "n_edges": len(ed), "mean_edge_OR_misfit_deg": float(np.mean(edge_misfits) if edge_misfits else np.nan), "inflation": inflation
    })


def variant_graph_reconstruction(
    grain_ids: Array,
    child_orientations: list[Array],
    edges: Iterable[tuple[int, int]],
    r_child_to_parent: Array,
    child_sym: Iterable[Array],
    parent_sym: Iterable[Array],
    sigma_deg: float = 3.5,
    inflation: float = 1.4,
    max_candidates_per_grain: int | None = None,
) -> ReconstructionResult:
    """Candidate-node variant graph with iterative probability propagation.

    Each possible parent orientation of each child grain is retained. Pairwise
    candidate compatibility across neighboring grains is converted to a probability-like
    weight, and support is propagated iteratively before one parent candidate is selected
    per grain. This retains second/third-best candidate information that a one-node grain
    graph discards.
    """
    n = len(child_orientations)
    ed = _normalize_edges(edges, n)
    cands = candidate_sets(child_orientations, r_child_to_parent, child_sym, parent_sym)
    if max_candidates_per_grain is not None:
        cands = [x[:max_candidates_per_grain] for x in cands]
    neigh: list[list[int]] = [[] for _ in range(n)]
    pair_weights: dict[tuple[int, int], Array] = {}
    for i, j in ed:
        neigh[i].append(j); neigh[j].append(i)
        w = np.zeros((len(cands[i]), len(cands[j])), float)
        for ai, (gi, _) in enumerate(cands[i]):
            for aj, (gj, _) in enumerate(cands[j]):
                ang = misorientation_deg(gi, gj, parent_sym)
                w[ai, aj] = _gaussian_weight(ang, sigma_deg, sigma_deg)
        pair_weights[(i, j)] = w

    probs = [np.ones(len(x), float) / max(len(x), 1) for x in cands]
    # inflation sharpens the probability distribution while preserving all candidates.
    for _ in range(35):
        new_probs: list[Array] = []
        max_change = 0.0
        for i in range(n):
            score = np.full(len(cands[i]), 1e-12)
            for j in neigh[i]:
                if i < j:
                    w = pair_weights[(i, j)]
                    msg = w @ probs[j]
                else:
                    w = pair_weights[(j, i)]
                    msg = w.T @ probs[j]
                score += msg
            score = np.power(np.maximum(score, 1e-15), max(inflation, 1.0))
            score /= np.sum(score)
            max_change = max(max_change, float(np.max(np.abs(score - probs[i]))))
            new_probs.append(score)
        probs = new_probs
        if max_change < 1e-7:
            break

    selected: list[Array] = []
    ambiguity: list[float] = []
    for i in range(n):
        order = np.argsort(probs[i])[::-1]
        selected.append(cands[i][int(order[0])][0])
        second = probs[i][int(order[1])] if len(order) > 1 else 0.0
        ambiguity.append(float(probs[i][int(order[0])] - second))
    pids = _clusters_from_selected(selected, ed, parent_sym, threshold_deg=max(4.0, 1.5 * sigma_deg))
    return _finalize("Variant graph", grain_ids, child_orientations, cands, pids, parent_sym, {
        "n_edges": len(ed),
        "candidate_nodes": int(sum(len(x) for x in cands)),
        "mean_probability_margin": float(np.mean(ambiguity)),
        "inflation": inflation,
    })

def nucleation_growth_reconstruction(
    grain_ids: Array,
    child_orientations: list[Array],
    edges: Iterable[tuple[int, int]],
    r_child_to_parent: Array,
    child_sym: Iterable[Array],
    parent_sym: Iterable[Array],
    nucleation_deg: float = 3.0,
    growth_deg: float = 8.0,
) -> ReconstructionResult:
    n = len(child_orientations)
    ed = _normalize_edges(edges, n)
    cands = candidate_sets(child_orientations, r_child_to_parent, child_sym, parent_sym)
    neigh: list[list[int]] = [[] for _ in range(n)]
    for i, j in ed:
        neigh[i].append(j); neigh[j].append(i)
    # local best candidate and support
    best_cand: list[Array] = []
    support: list[int] = []
    for i in range(n):
        best = (-1, 1e9, cands[i][0][0])
        for cand, _ in cands[i]:
            fits = [min(misorientation_deg(cand, cj, parent_sym) for cj, _ in cands[j]) for j in neigh[i]]
            count = sum(x <= nucleation_deg for x in fits)
            mean = np.mean(fits) if fits else 999.0
            key = (count, -mean)
            if key > (best[0], -best[1]):
                best = (count, float(mean), cand)
        support.append(best[0]); best_cand.append(best[2])
    order = np.argsort(support)[::-1]
    parent_ids = np.zeros(n, dtype=int)
    parents: dict[int, Array] = {}
    pid = 0
    # nucleate from highly supported grains, preventing duplicate seeds
    for i in order:
        if support[i] < max(1, len(neigh[i]) // 2):
            continue
        if parent_ids[i] != 0:
            continue
        existing = next((p for p, g in parents.items() if misorientation_deg(best_cand[i], g, parent_sym) <= growth_deg), None)
        if existing is None:
            pid += 1; parents[pid] = best_cand[i]; existing = pid
        parent_ids[i] = existing
    if not parents:
        # fallback to highest-support grain
        pid = 1; parents[1] = best_cand[int(order[0])]; parent_ids[int(order[0])] = 1
    # grow iteratively from assigned neighbors
    changed = True
    while changed:
        changed = False
        for i in range(n):
            if parent_ids[i] != 0:
                continue
            votes: dict[int, tuple[float, Array]] = {}
            for j in neigh[i]:
                pj = int(parent_ids[j])
                if pj == 0:
                    continue
                gp = parents[pj]
                vals = [(misorientation_deg(c, gp, parent_sym), c) for c, _ in cands[i]]
                fit, cand = min(vals, key=lambda x: x[0])
                if fit <= growth_deg and (pj not in votes or fit < votes[pj][0]):
                    votes[pj] = (fit, cand)
            if votes:
                chosen = min(votes, key=lambda p: votes[p][0])
                parent_ids[i] = chosen; changed = True
    # remaining isolated grains become singleton parents
    for i in range(n):
        if parent_ids[i] == 0:
            pid += 1; parent_ids[i] = pid; parents[pid] = best_cand[i]
    return _finalize("Nucleation + growth", grain_ids, child_orientations, cands, parent_ids, parent_sym, {
        "n_edges": len(ed), "nucleation_deg": nucleation_deg, "growth_deg": growth_deg
    })


def theoretical_child_operators(
    r_child_to_parent: Array,
    parent_sym: Iterable[Array],
    child_sym: Iterable[Array],
    tol_deg: float = 1e-4,
) -> list[Array]:
    """Return ideal child-child operators induced by the parent symmetry group.

    For child variants generated from one parent, relative operators reduce to
    R_cp^T S_parent R_cp. Child-symmetry equivalents are handled during matching.
    """
    ops: list[Array] = []
    for sp in parent_sym:
        op = _project_so3(r_child_to_parent.T @ sp @ r_child_to_parent)
        # inexpensive one-sided dedup is sufficient here; exact two-sided equivalence
        # is accounted for by the matching routine below.
        if not any(rotation_angle_deg(op @ q.T) <= tol_deg for q in ops):
            ops.append(op)
    return ops


def two_sided_operator_distance_deg(measured: Array, target: Array, child_sym: Iterable[Array]) -> float:
    """Distance between child-child operators modulo symmetry of both child grains."""
    m = np.asarray(measured, float)
    t = np.asarray(target, float)
    syms = np.asarray(tuple(child_sym), float)
    best_trace = -1.0
    # vectorize over the right symmetry and loop only over the left symmetry.
    for s1 in syms:
        a = t.T @ s1 @ m
        traces = np.einsum("ij,kji->k", a, syms)
        best_trace = max(best_trace, float(np.max(traces)))
    c = np.clip((best_trace - 1.0) / 2.0, -1.0, 1.0)
    return float(np.degrees(np.arccos(c)))


def operator_groupoid_reconstruction(
    grain_ids: Array,
    child_orientations: list[Array],
    edges: Iterable[tuple[int, int]],
    r_child_to_parent: Array,
    child_sym: Iterable[Array],
    parent_sym: Iterable[Array],
    operator_tol_deg: float = 4.0,
    parent_consistency_deg: float = 5.0,
) -> ReconstructionResult:
    """Operator/groupoid-style reconstruction with consistency-controlled growth.

    Daughter-daughter misorientations are first classified against the theoretical
    operator set induced by the OR and parent symmetry. Because a local operator match
    alone can occur accidentally across a true parent boundary, qualifying edges are
    subsequently resolved by common-parent candidate voting before clusters are grown.
    """
    n = len(child_orientations)
    ed = _normalize_edges(edges, n)
    cands = candidate_sets(child_orientations, r_child_to_parent, child_sym, parent_sym)
    ops = theoretical_child_operators(r_child_to_parent, parent_sym, child_sym)
    qualified: list[tuple[int, int]] = []
    residuals = []
    for i, j in ed:
        measured = _project_so3(child_orientations[i].T @ child_orientations[j])
        r = min(two_sided_operator_distance_deg(measured, op, child_sym) for op in ops)
        residuals.append(r)
        if r <= operator_tol_deg:
            qualified.append((i, j))

    neigh: list[list[int]] = [[] for _ in range(n)]
    for i, j in qualified:
        neigh[i].append(j); neigh[j].append(i)
    selected: list[Array] = []
    for i in range(n):
        if not neigh[i]:
            selected.append(cands[i][0][0])
            continue
        local = []
        for cand, _ in cands[i]:
            fits = [min(misorientation_deg(cand, cj, parent_sym) for cj, _ in cands[j]) for j in neigh[i]]
            count = sum(f <= parent_consistency_deg for f in fits)
            score = (count, -float(np.mean(fits)))
            local.append((score, cand))
        selected.append(max(local, key=lambda x: x[0])[1])
    pids = _clusters_from_selected(selected, qualified, parent_sym, threshold_deg=parent_consistency_deg)
    # grains with no qualifying operator edges remain distinct rather than being forced.
    return _finalize("Operator / groupoid consistency", grain_ids, child_orientations, cands, pids, parent_sym, {
        "n_theoretical_operators": len(ops),
        "n_operator_consistent_edges": len(qualified),
        "mean_operator_residual_deg": float(np.mean(residuals) if residuals else np.nan),
        "operator_tol_deg": operator_tol_deg,
        "parent_consistency_deg": parent_consistency_deg,
    })


# ---------- OR refinement ----------

def refine_orientation_relationship(
    child_orientations: list[Array],
    edges: Iterable[tuple[int, int]],
    initial_r_child_to_parent: Array,
    child_sym: Iterable[Array],
    parent_sym: Iterable[Array],
    max_rotation_deg: float = 5.0,
    max_edges: int = 400,
) -> tuple[Array, dict[str, float]]:
    """Refine an OR with a deterministic bounded rotation-vector pattern search.

    The objective is the robust neighboring-grain candidate-parent mismatch.  A
    small derivative-free coordinate search is used instead of an unconstrained
    optimizer because the candidate-switching objective is non-smooth.  This makes
    the result deterministic, bounded, and fast enough for interactive/CI use.
    """
    n = len(child_orientations)
    ed = _normalize_edges(edges, n)
    if not ed:
        raise ValueError("OR refinement requires grain adjacency.")
    if len(ed) > max_edges:
        ids = np.linspace(0, len(ed) - 1, max_edges).astype(int)
        ed = [ed[k] for k in ids]

    max_rad = np.deg2rad(float(max_rotation_deg))

    def objective(rv: Array) -> float:
        angle = np.linalg.norm(rv)
        if angle > max_rad:
            return 1e4 + angle * 1e3
        r = Rotation.from_rotvec(rv).as_matrix() @ initial_r_child_to_parent
        cs = candidate_sets(child_orientations, r, child_sym, parent_sym)
        vals = np.array([_best_pair_misfit(cs[i], cs[j], parent_sym) for i, j in ed])
        delta = 2.0
        loss = delta * delta * (np.sqrt(1.0 + (vals / delta) ** 2) - 1.0)
        return float(np.mean(loss))

    rv = np.zeros(3, dtype=float)
    before = objective(rv)
    best = before
    evaluations = 1

    # Coarse-to-fine coordinate pattern search.  The final 0.1° step is well below
    # normal EBSD OR-refinement uncertainties while keeping CI/runtime predictable.
    steps_deg = (1.0, 0.5, 0.25, 0.10)
    for step_deg in steps_deg:
        step = np.deg2rad(step_deg)
        for _ in range(3):
            candidate_best = best
            candidate_rv = rv
            for axis in range(3):
                for sign in (-1.0, 1.0):
                    trial = rv.copy()
                    trial[axis] += sign * step
                    if np.linalg.norm(trial) > max_rad:
                        continue
                    value = objective(trial)
                    evaluations += 1
                    if value < candidate_best - 1e-12:
                        candidate_best = value
                        candidate_rv = trial
            if candidate_best < best - 1e-12:
                best = candidate_best
                rv = candidate_rv
            else:
                break

    out = Rotation.from_rotvec(rv).as_matrix() @ initial_r_child_to_parent
    return _project_so3(out), {
        "objective_before": float(before),
        "objective_after": float(best),
        "correction_angle_deg": float(np.degrees(np.linalg.norm(rv))),
        "evaluations": float(evaluations),
        "success": float(best <= before),
    }


# ---------- dataframe / synthetic helpers ----------

def orientations_from_dataframe(
    df: pd.DataFrame,
    convention: Literal["crystal_to_specimen", "specimen_to_crystal"] = "crystal_to_specimen",
) -> tuple[Array, list[Array]]:
    """Read grain orientations with an explicit matrix-direction convention.

    Internally the reconstruction engine always uses proper rotations mapping crystal
    coordinates to specimen coordinates. Some EBSD packages export the transpose
    convention; those matrices are explicitly transposed on import rather than silently
    mixed with the internal convention.
    """
    if "grain_id" not in df.columns:
        raise ValueError("Orientation CSV must contain a grain_id column.")
    if convention not in {"crystal_to_specimen", "specimen_to_crystal"}:
        raise ValueError("Unknown orientation convention.")
    grain_ids = pd.to_numeric(df["grain_id"], errors="raise").astype(int).to_numpy()
    if len(set(grain_ids.tolist())) != len(grain_ids):
        raise ValueError("grain_id values must be unique.")
    if {"phi1_deg", "Phi_deg", "phi2_deg"}.issubset(df.columns):
        ori = [bunge_euler_to_matrix(float(a), float(b), float(c)) for a, b, c in df[["phi1_deg", "Phi_deg", "phi2_deg"]].to_numpy()]
    elif {"qw", "qx", "qy", "qz"}.issubset(df.columns):
        ori = [matrix_from_quaternion_wxyz(float(a), float(b), float(c), float(d)) for a, b, c, d in df[["qw", "qx", "qy", "qz"]].to_numpy()]
    else:
        raise ValueError("Provide either Bunge Euler columns phi1_deg, Phi_deg, phi2_deg or quaternion columns qw,qx,qy,qz.")
    if convention == "specimen_to_crystal":
        ori = [g.T for g in ori]
    return grain_ids, ori


def edges_from_dataframe(edge_df: pd.DataFrame, grain_ids: Array) -> list[tuple[int, int]]:
    cols = set(edge_df.columns)
    if not ({"grain_id_1", "grain_id_2"} <= cols):
        raise ValueError("Adjacency CSV must contain grain_id_1 and grain_id_2.")
    pos = {int(g): i for i, g in enumerate(grain_ids)}
    out = []
    for a, b in edge_df[["grain_id_1", "grain_id_2"]].to_numpy():
        if int(a) not in pos or int(b) not in pos:
            raise ValueError("Adjacency references a grain_id not present in the orientation table.")
        out.append((pos[int(a)], pos[int(b)]))
    return _normalize_edges(out, len(grain_ids))


def approximate_knn_edges(df: pd.DataFrame, k: int = 4) -> list[tuple[int, int]]:
    if not {"x", "y"}.issubset(df.columns):
        raise ValueError("Approximate k-NN adjacency requires x and y centroid columns.")
    xy = df[["x", "y"]].to_numpy(float)
    n = len(xy)
    k = max(1, min(int(k), n - 1))
    out: set[tuple[int, int]] = set()
    for i in range(n):
        d2 = np.sum((xy - xy[i]) ** 2, axis=1)
        ids = np.argsort(d2)[1:k+1]
        for j in ids:
            out.add((min(i, int(j)), max(i, int(j))))
    return sorted(out)


def synthetic_parent_reconstruction_demo(
    preset_name: str = "Kurdjumov–Sachs (FCC parent → BCC child)",
    n_per_parent: int = 12,
    noise_deg: float = 0.35,
    seed: int = 11,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[int, Array]]:
    preset = orientation_relationship_presets()[preset_name]
    ps = symmetry_group(preset.parent_symmetry)
    cs = symmetry_group(preset.child_symmetry)
    rng = np.random.default_rng(seed)
    # two well-separated parent orientations
    p1 = Rotation.from_euler("ZXZ", [12, 28, 41], degrees=True).as_matrix()
    p2 = Rotation.from_euler("ZXZ", [126, 52, 203], degrees=True).as_matrix()
    parents = {1: p1, 2: p2}
    rows = []
    edges = []
    gid = 1
    for pid, gp in parents.items():
        variants = unique_child_variants(gp, preset.matrix_child_to_parent, ps, cs)
        start = gid
        for k in range(n_per_parent):
            g = variants[k % len(variants)]
            axis = rng.normal(size=3); axis /= np.linalg.norm(axis)
            ang = np.deg2rad(rng.normal(0, noise_deg))
            noisy = Rotation.from_rotvec(axis * ang).as_matrix() @ g
            e1, E, e2 = matrix_to_bunge_euler(noisy)
            rows.append({"grain_id": gid, "phi1_deg": e1, "Phi_deg": E, "phi2_deg": e2, "x": float(k % 4 + (pid-1)*6), "y": float(k // 4), "true_parent_id": pid})
            gid += 1
        # chain + short-range edges within each parent
        ids = list(range(start, gid))
        for a, b in zip(ids[:-1], ids[1:]):
            edges.append({"grain_id_1": a, "grain_id_2": b})
        for a, b in zip(ids[:-4], ids[4:]):
            edges.append({"grain_id_1": a, "grain_id_2": b})
    # one boundary edge between parent grains to challenge algorithms
    edges.append({"grain_id_1": n_per_parent, "grain_id_2": n_per_parent + 1})
    return pd.DataFrame(rows), pd.DataFrame(edges), parents


def reconstruction_accuracy_against_labels(result: ReconstructionResult, true_labels: Array) -> float:
    """Permutation-invariant clustering accuracy for validation datasets."""
    pred = np.asarray(result.parent_ids, int)
    true = np.asarray(true_labels, int)
    if len(pred) != len(true):
        raise ValueError("Label length mismatch.")
    # map each predicted cluster to its majority truth label
    correct = 0
    for p in set(pred.tolist()):
        idx = pred == p
        vals, counts = np.unique(true[idx], return_counts=True)
        correct += int(np.max(counts))
    return float(correct / len(true))
