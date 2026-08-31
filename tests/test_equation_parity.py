from __future__ import annotations

import numpy as np

from src.core import LatticeInput, cmc_degeneracy, cmc_matrix, normalized_metrics
from src.equation_engine import (
    appendix_c_o4,
    cmc_quadratic_explicit,
    cmc_symmetry_group,
    correspondence_left_cosets,
    first_order_families,
    higher_order_degeneracy_families,
    ips_geometry_from_stretch,
    lambdas_from_stretch_shear,
    niti_cmc_analytic,
    niti_cmc_from_input,
    o2_type_i_family,
    o2_type_ii_family,
    paper_defined_distances,
    stretch_shear_from_lambdas,
    verify_analytic_cmc_against_general,
)
from src.provenance import equation_registry
from src.ptmc import stretch_from_lattice
from src.symmetry import correspondence_intersection_subgroup, cubic_point_group


def kudoh() -> LatticeInput:
    return LatticeInput(3.01, 2.898, 4.108, 4.646, 97.78)


def test_analytic_cmc_exactly_matches_general_metric_construction():
    inp = kudoh()
    chk = verify_analytic_cmc_against_general(inp)
    assert chk["matrix_frobenius_residual"] < 1e-14
    assert chk["eigenvalue_max_abs_residual"] < 1e-14


def test_explicit_quadratic_is_four_times_cmc_form():
    inp = kudoh()
    a, b, c = inp.ratios()
    q = niti_cmc_from_input(inp).matrix
    for u in ([1, 0, 0], [1, -2, 0.3], [-0.4, 0.7, 1.2]):
        assert np.isclose(cmc_quadratic_explicit(u, a, b, c, inp.beta_deg), 4*np.asarray(u)@q@np.asarray(u), atol=1e-12)


def test_closed_form_eigenvalues_match_numerical_eigenvalues():
    for args in [
        (0.9628, 1.3648, 1.5435, 97.78),
        (0.9, np.sqrt(2), 1.36201, 97.78),
        (1.05, 1.39, 1.8, 104.0),
    ]:
        r = niti_cmc_analytic(*args)
        assert np.allclose(r.eigenvalues_sorted, np.linalg.eigvalsh(r.matrix), atol=1e-12)


def test_c1_constructed_point_has_first_order_or_higher_degeneracy():
    a, c, beta = 0.9628, 1.5435, 97.78
    fam = {x.name: x for x in first_order_families(a, np.sqrt(2), c, beta)}
    assert fam["C1"].met
    q = niti_cmc_analytic(a, np.sqrt(2), c, beta).matrix
    d = cmc_degeneracy(q, rtol=1e-8)
    assert d.order >= 1
    assert d.habit_planes


def test_c2a_constructed_example_from_source_is_detected():
    # The source figure reports rounded c=2.1451. Build the exact analytical target
    # from its displayed a,b,beta so this test checks the equation rather than
    # the figure's limited decimal precision.
    a, b, beta = 1.0166, 1.3648, 97.78
    s = np.sin(np.deg2rad(beta))
    c = np.sqrt(2.0)*np.sqrt((1.0-a*a)/(1.0-a*a*s*s))
    assert abs(c-2.1451) < 2e-3
    fam = {x.name: x for x in first_order_families(a, b, c, beta, tol=1e-9)}
    assert fam["C2a"].met


def test_c2b_source_discrepancy_is_not_silently_reconciled():
    # The plotted analytical C2b example has b<sqrt(2).
    eq = {x.name: x for x in first_order_families(1.0, 1.3142, 1.5409, 90.0, c2b_interpretation="equation")}
    table = {x.name: x for x in first_order_families(1.0, 1.3142, 1.5409, 90.0, c2b_interpretation="table")}
    assert eq["C2b"].met
    assert not table["C2b"].met
    assert "conflict" in table["C2b"].note.lower()


def test_c3_constructed_example_from_source_is_detected():
    fam = {x.name: x for x in first_order_families(0.98, 1.4242, 1.1767, 97.78, tol=2e-4)}
    assert fam["C3"].met


def test_second_and_third_order_families():
    d1 = {x.name: x for x in higher_order_degeneracy_families(0.9, np.sqrt(2), 1.36201, 97.78, tol=3e-5)}
    assert d1["D1"].met and d1["D1"].order == 2 and d1["D1"].habit_plane is not None
    # Source figure gives x-y-2.290z=0 for this branch.
    p = d1["D1"].habit_plane / d1["D1"].habit_plane[0]
    assert np.allclose(p, [1, -1, -2.290], atol=2e-3)

    d2 = {x.name: x for x in higher_order_degeneracy_families(1.0, 1.2, np.sqrt(2), 90.0)}
    assert d2["D2"].met and d2["D2"].order == 2
    assert np.allclose(d2["D2"].habit_plane, [1, 1, 0])
    assert not d2["E"].met

    e = {x.name: x for x in higher_order_degeneracy_families(1.0, np.sqrt(2), np.sqrt(2), 90.0)}
    assert e["E"].met and e["E"].order == 3
    assert np.allclose(niti_cmc_analytic(1, np.sqrt(2), np.sqrt(2), 90).matrix, 0, atol=1e-14)


def test_table_distance_numbers_and_their_internal_discrepancy_are_reproduced():
    a, b, c = kudoh().ratios()  # the printed table values are generated from the raw lattice constants
    d = paper_defined_distances(a, b, c, 97.78)
    assert np.isclose(d["C1"]["distance_equality"], 0.137364, atol=1e-6)
    assert np.isclose(d["C1"]["distance_inequality_frontier"], 0.068402, atol=1e-6)
    assert np.isclose(d["C2a"]["distance_inequality_frontier"], 0.227385, atol=1e-6)
    assert np.isclose(d["C2b"]["distance_equality"], 0.091359, atol=1e-6)
    assert np.isclose(d["C3"]["distance_inequality_frontier"], 0.210399, atol=1e-6)
    # The source's 0.269706 equality value is reproduced by |c-c_target|, not its parsed/visual squared expression.
    assert np.isclose(d["C3"]["distance_equality_table_value_reproduction"], 0.269706, atol=1e-6)
    assert d["C3"]["distance_equality_printed_expression"] > 0.7


def test_appendix_c_beta_98_closed_form_and_invariant_length():
    a, c = appendix_c_o4(98.0)
    assert np.isclose(a, 0.8825126216, atol=1e-9)
    assert np.isclose(c, 1.6182339082, atol=1e-9)
    beta = np.deg2rad(98.0)
    assert np.isclose(a*a+c*c+2*a*c*np.cos(beta), 3.0, atol=1e-12)
    # C1 is also required and b=sqrt(2).
    fam = {x.name: x for x in first_order_families(a, np.sqrt(2), c, 98.0)}
    assert fam["C1"].met


def test_o2_analytic_families_are_c1_compatible():
    for builder in (o2_type_i_family, o2_type_ii_family):
        a, b, c = builder(98.0)
        fam = {x.name: x for x in first_order_families(a, b, c, 98.0, tol=1e-9)}
        assert fam["C1"].met


def test_stretch_shear_round_trip_and_explicit_ips_rank_one_geometry():
    l1, l3 = 0.93, 1.1195
    tau, delta = stretch_shear_from_lambdas(l1, l3)
    l1b, l3b = lambdas_from_stretch_shear(tau, delta)
    assert np.allclose([l1b, l3b], [l1, l3], atol=1e-12)

    inp = kudoh()
    # Force C1 exactly while keeping other physical lattice values.
    c1 = LatticeInput(inp.a_b2, inp.a_b19p, np.sqrt(2)*inp.a_b2, inp.c_b19p, inp.beta_deg)
    U = stretch_from_lattice(c1).U
    for branch in (-1, 1):
        g = ips_geometry_from_stretch(U, branch)
        assert g.rank_one_residual < 1e-12
        assert g.plane_invariance_residual < 1e-12
        assert abs(g.trace_identity_residual) < 1e-12
        assert abs(g.determinant_identity_residual) < 1e-12


def test_correspondence_left_cosets_partition_parent_group_and_count_12():
    cosets = correspondence_left_cosets()
    H = correspondence_intersection_subgroup()
    G = cubic_point_group()
    assert len(G) == 48 and len(H) == 4 and len(cosets) == 12
    assert all(len(c.elements) == 4 for c in cosets)
    keys = [tuple(np.rint(g).astype(int).ravel()) for c in cosets for g in c.elements]
    assert len(keys) == len(set(keys)) == 48


def test_cmc_symmetry_contains_correspondence_subgroup_for_binary_benchmark():
    ma, mm = normalized_metrics(kudoh())
    entries = cmc_symmetry_group(cmc_matrix(ma, mm))
    assert len(entries) >= 4
    found = [e.symmetry for e in entries]
    for h in correspondence_intersection_subgroup():
        assert any(np.allclose(h, g) for g in found)


def test_equation_registry_covers_new_analytical_engine_and_discrepancies():
    reg = equation_registry()
    required = {
        "PTMC-SC1", "PTMC-SC2", "PTMC-SC3", "PTMC-V", "PTMC-TRACE", "PTMC-DET",
        "HC", "LEFT-COSET", "N-CORR", "DOUBLE-COSET", "CMC", "G-CMC", "SMC",
        "SHEAR-SHEAR", "EPSILON", "NITI-CMC", "NITI-Q", "CT-C1", "CT-C2a",
        "CT-C2b-EQ", "CT-C2b-TABLE", "CT-C3", "CT-D1", "CT-D2", "CT-E", "APP-C-O4", "FTC",
    }
    assert required <= set(reg)
    assert reg["CT-C2b-TABLE"].provenance_class == "source discrepancy"
    assert "ambiguous" in (reg["NITI-Q"].caveat + reg["NITI-Q"].scope).lower()
