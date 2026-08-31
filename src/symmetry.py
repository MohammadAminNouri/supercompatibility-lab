from __future__ import annotations

from dataclasses import dataclass
from itertools import permutations, product
from math import gcd
from functools import reduce

import numpy as np

from .core import (
    Array,
    C_A_TO_M,
    C_M_TO_A,
    TwinResult,
    type_i_twin,
    type_ii_twin,
)


@dataclass(frozen=True)
class SymmetryElement:
    matrix: Array
    kind: str
    axis_or_normal: Array | None
    notation: str


@dataclass(frozen=True)
class DoubleCoset:
    label: str
    matrices: tuple[Array, ...]
    representative: Array
    polar: bool
    contains_identity: bool
    n_mirrors: int
    n_twofold_rotations: int


@dataclass(frozen=True)
class TwinExplorerEntry:
    double_coset: str
    symmetry_kind: str
    parent_element: str
    twin: TwinResult


def _mat_key(m: Array) -> tuple[int, ...]:
    return tuple(np.rint(np.asarray(m)).astype(int).ravel().tolist())


def cubic_point_group() -> tuple[Array, ...]:
    """Full m-3m cubic point group as 48 signed permutation matrices."""
    mats: dict[tuple[int, ...], Array] = {}
    for perm in permutations(range(3)):
        p = np.zeros((3, 3), dtype=float)
        p[np.arange(3), perm] = 1.0
        for signs in product((-1.0, 1.0), repeat=3):
            g = np.diag(signs) @ p
            mats[_mat_key(g)] = g
    return tuple(mats[k] for k in sorted(mats))


def monoclinic_2m_point_group() -> tuple[Array, ...]:
    """2/m point group in the conventional monoclinic basis with unique b axis."""
    return (
        np.eye(3),
        np.diag([-1.0, 1.0, -1.0]),  # two-fold about b
        np.diag([1.0, -1.0, 1.0]),   # mirror normal to b
        -np.eye(3),                   # inversion
    )


def correspondence_intersection_subgroup() -> tuple[Array, ...]:
    """H_C^A = G_A ∩ C_A→M G_M C_M→A for the built-in B2→B19' model."""
    ga = {_mat_key(g): g for g in cubic_point_group()}
    out: dict[tuple[int, ...], Array] = {}
    for gm in monoclinic_2m_point_group():
        mapped = C_A_TO_M @ gm @ C_M_TO_A
        k = _mat_key(mapped)
        if k in ga and np.allclose(mapped, ga[k], atol=1e-12):
            out[k] = ga[k]
    if len(out) != 4:
        raise RuntimeError(f"Expected a 4-element correspondence subgroup, got {len(out)}.")
    return tuple(out[k] for k in sorted(out))


def symmetry_kind(g: Array) -> str:
    m = np.asarray(g, float)
    det = int(round(np.linalg.det(m)))
    tr = float(np.trace(m))
    if np.allclose(m, np.eye(3), atol=1e-10):
        return "identity"
    if np.allclose(m, -np.eye(3), atol=1e-10):
        return "inversion"
    if np.allclose(m @ m, np.eye(3), atol=1e-10):
        if det == -1 and abs(tr - 1.0) < 1e-8:
            return "mirror"
        if det == 1 and abs(tr + 1.0) < 1e-8:
            return "rotation_180"
    if det == 1:
        return "rotation_other"
    return "improper_other"


def _canonical_integer_direction(v: Array) -> Array:
    x = np.asarray(v, float).copy()
    x[np.abs(x) < 1e-8] = 0.0
    nz = np.abs(x[np.abs(x) > 1e-8])
    if len(nz) == 0:
        return np.zeros(3, dtype=int)
    x = x / np.min(nz)
    xr = np.rint(x).astype(int)
    if not np.allclose(x, xr, atol=1e-5):
        # Fallback for numerical eigenvectors; cubic two-fold elements should
        # still be representable by small integers.
        x = x / max(np.max(np.abs(x)), 1e-12)
        xr = np.rint(8 * x).astype(int)
    vals = [abs(int(a)) for a in xr if int(a) != 0]
    if vals:
        d = reduce(gcd, vals)
        xr //= max(d, 1)
    for a in xr:
        if a != 0:
            if a < 0:
                xr *= -1
            break
    return xr


def axis_or_normal(g: Array) -> Array | None:
    kind = symmetry_kind(g)
    if kind not in {"mirror", "rotation_180"}:
        return None
    target = -1.0 if kind == "mirror" else 1.0
    vals, vecs = np.linalg.eigh(0.5 * (g + g.T))
    idx = int(np.argmin(np.abs(vals - target)))
    return _canonical_integer_direction(vecs[:, idx])


def symmetry_notation(g: Array) -> str:
    kind = symmetry_kind(g)
    v = axis_or_normal(g)
    if kind == "identity":
        return "identity"
    if kind == "inversion":
        return "inversion"
    if v is not None:
        body = " ".join(str(int(x)) for x in v)
        return f"({body}) mirror" if kind == "mirror" else f"[{body}] 180° axis"
    return kind.replace("_", " ")


def double_cosets() -> tuple[DoubleCoset, ...]:
    """Partition G_A into H_C^A g H_C^A double cosets.

    This is the symmetry reduction used to avoid redundant variant-pair/twin
    calculations. For the built-in B2→B19' correspondence it yields 7 classes.
    """
    ga = {_mat_key(g): g for g in cubic_point_group()}
    h = correspondence_intersection_subgroup()
    remaining = set(ga)
    raw: list[set[tuple[int, ...]]] = []
    while remaining:
        k0 = min(remaining)
        g = ga[k0]
        dc: set[tuple[int, ...]] = set()
        for h1 in h:
            for h2 in h:
                dc.add(_mat_key(h1 @ g @ h2))
        raw.append(dc)
        remaining -= dc

    def properties(keys: set[tuple[int, ...]]) -> tuple[bool, bool, int, int, tuple[int, ...]]:
        kinds = [symmetry_kind(ga[k]) for k in keys]
        contains_identity = any(k == _mat_key(np.eye(3)) for k in keys)
        nmir = kinds.count("mirror")
        nrot = kinds.count("rotation_180")
        polar = (nmir + nrot) == 0
        rep = min(keys)
        return contains_identity, polar, nmir, nrot, rep

    # Stable, physically meaningful ordering: identity class first, then polar
    # classes, then remaining ambivalent classes.
    raw.sort(key=lambda keys: (
        0 if properties(keys)[0] else 1,
        0 if properties(keys)[1] else 1,
        properties(keys)[4],
    ))

    out: list[DoubleCoset] = []
    for i, keys in enumerate(raw, 1):
        contains_identity, polar, nmir, nrot, rep_key = properties(keys)
        mats = tuple(ga[k] for k in sorted(keys))
        out.append(DoubleCoset(
            label=f"DC{i}",
            matrices=mats,
            representative=ga[rep_key],
            polar=polar,
            contains_identity=contains_identity,
            n_mirrors=nmir,
            n_twofold_rotations=nrot,
        ))
    return tuple(out)


def full_twin_explorer(m_a: Array, m_m: Array) -> list[TwinExplorerEntry]:
    """All two-fold-generated Type-I/Type-II CT twin candidates by double coset.

    Symmetry-equivalent parent elements are retained and labeled, while the
    double-coset column makes redundancy explicit instead of silently dropping it.
    """
    entries: list[TwinExplorerEntry] = []
    seen: set[tuple[str, str, tuple[int, ...]]] = set()
    for dc in double_cosets():
        for g in dc.matrices:
            kind = symmetry_kind(g)
            v = axis_or_normal(g)
            if v is None:
                continue
            k = (dc.label, kind, tuple(int(x) for x in v))
            if k in seen:
                continue
            seen.add(k)
            notation = symmetry_notation(g)
            try:
                if kind == "mirror":
                    twin = type_i_twin(m_a, m_m, v, f"{dc.label} Type I", notation)
                    if twin.shear_amplitude > 1e-12:
                        entries.append(TwinExplorerEntry(dc.label, kind, notation, twin))
                elif kind == "rotation_180":
                    twin = type_ii_twin(m_a, m_m, v, f"{dc.label} Type II", notation)
                    if twin.shear_amplitude > 1e-12:
                        entries.append(TwinExplorerEntry(dc.label, kind, notation, twin))
            except ValueError:
                continue
    return sorted(entries, key=lambda e: (e.double_coset, e.symmetry_kind, e.parent_element))


def stretch_variants(U: Array, tol: float = 1e-9) -> tuple[Array, ...]:
    """Unique parent-symmetry conjugates g U g^T."""
    variants: list[Array] = []
    for g in cubic_point_group():
        v = g @ U @ g.T
        if not any(np.allclose(v, w, atol=tol, rtol=0.0) for w in variants):
            variants.append(v)
    variants.sort(key=lambda m: tuple(np.round(m, 10).ravel()))
    return tuple(variants)


def commutator_norm(a: Array, b: Array) -> float:
    return float(np.linalg.norm(a @ b - b @ a, ord="fro"))


def commuting_variant_pairs(U: Array, tol: float = 1e-9) -> list[tuple[int, int, float]]:
    variants = stretch_variants(U)
    out: list[tuple[int, int, float]] = []
    for i in range(len(variants)):
        for j in range(i + 1, len(variants)):
            n = commutator_norm(variants[i], variants[j])
            if n <= tol:
                out.append((i + 1, j + 1, n))
    return out

def parent_twofold_axes() -> tuple[Array, ...]:
    """Unique <100> and <110> two-fold axes/normals in the cubic parent."""
    axes: list[Array] = []
    for g in cubic_point_group():
        v = axis_or_normal(g)
        if v is None:
            continue
        vf = np.asarray(v, float)
        vf /= np.linalg.norm(vf)
        if not any(np.allclose(vf, w, atol=1e-10) or np.allclose(vf, -w, atol=1e-10) for w in axes):
            axes.append(vf)
    axes.sort(key=lambda x: tuple(_canonical_integer_direction(x).tolist()))
    return tuple(axes)


def _collinear(v: Array, w: Array, tol: float = 1e-7) -> bool:
    """Return True when non-zero vectors are parallel or antiparallel."""
    a = np.asarray(v, float)
    b = np.asarray(w, float)
    na = float(np.linalg.norm(a))
    nb = float(np.linalg.norm(b))
    if na <= tol or nb <= tol:
        return False
    return abs(float((a @ b) / (na * nb))) >= 1.0 - tol


def compound_twin_pairs(entries: list[TwinExplorerEntry], tol: float = 1e-6) -> list[tuple[int, int]]:
    """Identify Type-I/Type-II entries that represent the same compound twin system.

    A compound twin is simultaneously a Type-I and Type-II solution for the same
    variant pair. Numerically, the two descriptions have the same shear plane
    normal and shear direction (up to sign/scaling) and the same shear amplitude.
    Candidate matching is restricted to the same intercorrespondence/double-coset
    class so unrelated geometrically parallel systems are not merged.
    """
    pairs: list[tuple[int, int]] = []
    for i, left in enumerate(entries):
        if left.twin.twin_type != "Type I" or left.twin.shear_amplitude <= tol:
            continue
        for j, right in enumerate(entries):
            if right.twin.twin_type != "Type II" or right.double_coset != left.double_coset:
                continue
            if abs(left.twin.shear_amplitude - right.twin.shear_amplitude) > tol * max(1.0, left.twin.shear_amplitude):
                continue
            if _collinear(left.twin.twin_shear_vector, right.twin.twin_shear_vector, tol) and _collinear(
                left.twin.twin_plane_normal, right.twin.twin_plane_normal, tol
            ):
                pairs.append((i, j))
    return pairs
