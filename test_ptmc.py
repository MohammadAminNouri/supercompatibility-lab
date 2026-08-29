import numpy as np

from src.core import LatticeInput
from src.distances import all_cofactor_systems, compatibility_dashboard
from src.ptmc import stretch_from_lattice


def test_binary_niti_ptmc_lambda2_matches_b_ratio_relation():
    inp = LatticeInput(3.01, 2.898, 4.108, 4.646, 97.78)
    s = stretch_from_lattice(inp)
    expected = (4.108 / 3.01) / np.sqrt(2.0)
    assert np.isclose(s.eigenvalues[1], expected, atol=2e-7)
    assert s.eigenvalues[0] < s.eigenvalues[1] < s.eigenvalues[2]


def test_c1_example_passes_classical_cofactor_for_some_systems():
    a_b2 = 3.01
    inp = LatticeInput(a_b2, 0.9628*a_b2, np.sqrt(2)*a_b2, 1.5435*a_b2, 97.78)
    systems = all_cofactor_systems(inp)
    passing = [s for s in systems if s.all_pass]
    assert passing
    assert any(s.domain_type.startswith("Compound") for s in passing)
    assert all(np.isnan(s.cc2_simplified_residual) for s in passing if s.domain_type.startswith("Compound"))


def test_cross_framework_discrepancy_is_exposed_for_c1_example():
    a_b2 = 3.01
    inp = LatticeInput(a_b2, 0.9628*a_b2, np.sqrt(2)*a_b2, 1.5435*a_b2, 97.78)
    d = compatibility_dashboard(inp)
    assert d.ptmc_all_pass
    assert d.best_epsilon is not None
    assert np.isclose(d.best_epsilon, 0.215, atol=0.003)


def test_domain_enumeration_classifies_compound_pair_without_double_counting():
    from src.ptmc import enumerate_domain_specs
    from src.symmetry import parent_twofold_axes

    inp = LatticeInput(3.01, 2.898, 4.108, 4.646, 97.78)
    U = stretch_from_lattice(inp).U
    specs = enumerate_domain_specs(U, parent_twofold_axes())
    compounds = [s for s in specs if s.domain_type.startswith("Compound")]
    assert len(specs) == 14
    assert len(compounds) == 2
    assert np.isclose(abs(float(compounds[0].axis @ compounds[0].partner_axis)), 0.0, atol=1e-12)


def test_noncompound_full_cc2_agrees_with_simplified_zero_condition_near_c1():
    # For non-compound Type-I/II systems the full and simplified CC2 conditions
    # should identify the same zero; the test avoids the compound axes where the
    # simplified Proposition-6 form is not applicable.
    a_b2 = 3.01
    inp = LatticeInput(a_b2, 0.9628*a_b2, np.sqrt(2)*a_b2, 1.5435*a_b2, 97.78)
    systems = [s for s in all_cofactor_systems(inp) if not s.domain_type.startswith("Compound")]
    assert systems
    # No non-compound system is accidentally reported as CC2-passing in this benchmark.
    assert all(not s.cc2_pass for s in systems)
    assert all(np.isfinite(s.cc2_simplified_residual) for s in systems)
