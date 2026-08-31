import numpy as np

from src.design import DesignWeights, inverse_lattice_design, pareto_lattice_scan
from src.presets import PRESETS
from src.ptmc import enumerate_domain_specs, stretch_from_lattice
from src.symmetry import parent_twofold_axes


def test_compound_domain_inverse_and_pareto_design_interfaces():
    inp = PRESETS["Published binary NiTi example"]
    assert inp is not None
    U = stretch_from_lattice(inp).U
    spec = next(s for s in enumerate_domain_specs(U, parent_twofold_axes()) if s.domain_type == "Compound 1")
    a, b, c = inp.ratios()
    bounds = {
        "a": (a - 0.05, a + 0.05),
        "b": (b - 0.05, b + 0.07),
        "c": (c - 0.05, c + 0.05),
        "beta_deg": (96.0, 100.0),
    }

    scan = pareto_lattice_scan(inp.a_b2, spec, bounds, n_samples=100, seed=7)
    assert len(scan) == 100
    assert {"pareto", "abs_lambda2_minus_1", "cc2_normalized", "cc3_penalty", "cmc_relative_zero"} <= set(scan.columns)
    assert scan["pareto"].any()

    res = inverse_lattice_design(
        inp.a_b2,
        spec,
        bounds,
        weights=DesignWeights(),
        reference=(a, b, c, inp.beta_deg),
        seed=7,
        maxiter=10,
    )
    assert res.domain_type == "Compound 1"
    assert res.partner_axis is not None
    assert np.isfinite(res.objective)
