import numpy as np

from src.compatibility_methods import (
    hadamard_rank_one_certificate,
    laminate_volume_fraction_scan,
    pairwise_twin_compatibility,
    single_variant_ips_residual,
    triplet_condition_orthorhombic,
    twin_rank_one_certificate,
)
from src.distances import all_cofactor_systems
from src.presets import PRESETS
from src.ptmc import stretch_from_lattice
from src.symmetry import stretch_variants


def test_hadamard_rank_one_certificate_exact_and_nonrankone():
    A = np.eye(3) + np.outer([1.0, 2.0, 0.0], [0.0, 1.0, 0.0])
    c = hadamard_rank_one_certificate(A, np.eye(3))
    assert c.exact_within_tol
    assert c.normalized_rank_one_residual < 1e-12
    B = np.diag([1.1, 0.9, 1.0])
    c2 = hadamard_rank_one_certificate(B, np.eye(3))
    assert not c2.exact_within_tol


def test_single_variant_ips_matches_middle_stretch():
    inp = PRESETS["Published binary NiTi example"]
    U = stretch_from_lattice(inp).U
    d = single_variant_ips_residual(U)
    assert np.isclose(d["lambda2"], stretch_from_lattice(inp).eigenvalues[1])


def test_pairwise_rank_one_test_finds_compatible_stretch_variants():
    U = stretch_from_lattice(PRESETS["Published binary NiTi example"]).U
    variants = stretch_variants(U)
    df = pairwise_twin_compatibility(variants)
    assert len(variants) == 12
    assert len(df) == 66
    assert int(df.rank_one_compatible.sum()) > 0
    row = df[df.rank_one_compatible].iloc[0]
    cert = twin_rank_one_certificate(variants[int(row.variant_i)-1], variants[int(row.variant_j)-1])
    assert cert.middle_residual < 1e-7


def test_cofactor_passing_system_has_all_volume_fraction_ips():
    inp = PRESETS["C1-compatible teaching example"]
    U = stretch_from_lattice(inp).U
    systems = [s for s in all_cofactor_systems(inp) if s.all_pass]
    assert systems
    for s in systems:
        scan = laminate_volume_fraction_scan(U, s.a, s.n, points=41)
        assert scan.max_middle_stretch_residual < 1e-10


def test_triplet_condition_specialized_formulas_have_exact_nontrivial_solutions():
    alpha = 0.95
    beta = 1.05
    gamma_i = np.sqrt(3 * alpha**2 * beta**2 / (alpha**2 + 2 * beta**2))
    r_i = triplet_condition_orthorhombic(alpha, beta, gamma_i)
    assert abs(r_i.tc_i_raw) < 1e-12
    gamma_ii = np.sqrt((2 * alpha**2 + beta**2) / 3.0)
    r_ii = triplet_condition_orthorhombic(alpha, beta, gamma_ii)
    assert abs(r_ii.tc_ii_raw) < 1e-12
