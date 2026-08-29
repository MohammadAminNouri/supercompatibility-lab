import numpy as np
import pandas as pd

from src.core import C_M_TO_A, LatticeInput
from src.frontier import frontier_diagnostic
from src.literature import load_literature
from src.multistep import GeneralCell, evaluate_stage
from src.ptmc import enumerate_domain_specs, stretch_from_lattice
from src.temperature import temperature_sweep
from src.symmetry import parent_twofold_axes
from src.uncertainty import LatticeUncertainty, monte_carlo_uncertainty


def test_frontier_applicability_guard_for_builtin_correspondence():
    inp = LatticeInput(3.01, 2.898, 4.108, 4.646, 97.78)
    d = frontier_diagnostic(stretch_from_lattice(inp).U)
    assert np.allclose(d.monoclinic_axis_parent_direction, [1,1,0])
    assert not d.built_in_is_monoclinic_ii
    assert d.extreme_target_distance_if_applicable is None
    assert d.n_stretch_variants == 12
    assert len(d.commuting_pairs) == 0


def test_generic_multistep_stage_reproduces_builtin_stretch():
    parent = GeneralCell(3.01,3.01,3.01,90,90,90)
    product = GeneralCell(2.898,4.108,4.646,90,97.78,90)
    r = evaluate_stage(parent, product, C_M_TO_A)
    s = stretch_from_lattice(LatticeInput(3.01,2.898,4.108,4.646,97.78))
    assert np.allclose([r.lambda1,r.lambda2,r.lambda3], s.eigenvalues, atol=1e-10)


def test_temperature_sweep_outputs_expected_columns():
    df = pd.DataFrame([
        {"temperature_K":300,"a_B2_A":3.01,"a_B19p_A":2.898,"b_B19p_A":4.108,"c_B19p_A":4.646,"beta_deg":97.78},
        {"temperature_K":310,"a_B2_A":3.011,"a_B19p_A":2.899,"b_B19p_A":4.110,"c_B19p_A":4.648,"beta_deg":97.76},
    ])
    out = temperature_sweep(df)
    assert len(out) == 2
    assert "abs_lambda2_minus_1" in out.columns
    assert np.isfinite(out["cofactor_cc2_normalized"]).all()


def test_uncertainty_is_deterministic_for_zero_sigmas():
    inp = LatticeInput(3.01, 0.9628*3.01, np.sqrt(2)*3.01, 1.5435*3.01, 97.78)
    stretch = stretch_from_lattice(inp)
    specs = enumerate_domain_specs(stretch.U, parent_twofold_axes())
    compound = next(s for s in specs if s.domain_type == "Compound 1")
    u = monte_carlo_uncertainty(inp, LatticeUncertainty(), compound, n=100, cc1_tol=1e-5, cc2_tol=1e-5, cmc_tol=1e-5)
    assert u.cc1_fraction == 1.0
    assert u.cc2_fraction == 1.0
    assert u.cc3_fraction == 1.0
    assert u.cmc_fraction == 1.0


def test_literature_database_is_curated_and_has_no_forbidden_author_name():
    df = load_literature()
    assert len(df) >= 12
    text = df.astype(str).to_string().lower()
    forbidden = "".join(chr(x) for x in (99, 97, 121, 114, 111, 110))
    assert forbidden not in text
    assert "10.1038/nature12532" in set(df["doi"])
    assert "10.1016/j.joule.2026.102627" in set(df["doi"])


def test_generic_cmc_distance_is_invariant_to_consistent_length_unit_scaling():
    from src.multistep import GeneralCell, evaluate_stage
    from src.core import C_M_TO_A

    parent_A = GeneralCell(3.01, 3.01, 3.01, 90, 90, 90)
    product_A = GeneralCell(2.898, 4.108, 4.646, 90, 97.78, 90)
    r_A = evaluate_stage(parent_A, product_A, C_M_TO_A)

    # Same cells expressed in nm instead of Å.
    parent_nm = GeneralCell(0.301, 0.301, 0.301, 90, 90, 90)
    product_nm = GeneralCell(0.2898, 0.4108, 0.4646, 90, 97.78, 90)
    r_nm = evaluate_stage(parent_nm, product_nm, C_M_TO_A)

    assert np.isclose(r_A.cmc_relative_zero, r_nm.cmc_relative_zero, rtol=1e-10, atol=1e-12)
    assert np.allclose([r_A.lambda1, r_A.lambda2, r_A.lambda3], [r_nm.lambda1, r_nm.lambda2, r_nm.lambda3], rtol=1e-10, atol=1e-12)
