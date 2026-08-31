from __future__ import annotations

"""Deterministic equation/source parity audit used by release CI."""

from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.core import LatticeInput
from src.equation_engine import (
    appendix_c_o4,
    correspondence_left_cosets,
    first_order_families,
    ips_geometry_from_stretch,
    paper_defined_distances,
    verify_analytic_cmc_against_general,
)
from src.ptmc import stretch_from_lattice
from src.provenance import equation_registry


def check(cond: bool, label: str) -> None:
    if not cond:
        raise AssertionError(label)
    print(f"PASS · {label}")


def main() -> None:
    inp = LatticeInput(3.01, 2.898, 4.108, 4.646, 97.78)
    parity = verify_analytic_cmc_against_general(inp)
    check(parity["matrix_frobenius_residual"] < 1e-13, "closed-form CMC == general metric CMC")
    check(parity["eigenvalue_max_abs_residual"] < 1e-13, "closed-form CMC eigenvalues == numerical eigenspectrum")

    a, b, c = inp.ratios()
    dist = paper_defined_distances(a, b, c, inp.beta_deg)
    check(abs(float(dist["C1"]["distance_equality"]) - 0.137364) < 1e-6, "source Table-2 C1 equality checkpoint")
    check(abs(float(dist["C1"]["distance_inequality_frontier"]) - 0.068402) < 1e-6, "source Table-2 C1 frontier checkpoint")
    check(abs(float(dist["C3"]["distance_equality_table_value_reproduction"]) - 0.269706) < 1e-6, "source Table-2 0.269706 numerical-value reproduction")
    check(float(dist["C3"]["distance_equality_printed_expression"]) > 0.7, "source distance expression/value discrepancy remains detectable")

    eq = {x.name: x for x in first_order_families(1.0, 1.3142, 1.5409, 90.0, c2b_interpretation="equation")}
    tab = {x.name: x for x in first_order_families(1.0, 1.3142, 1.5409, 90.0, c2b_interpretation="table")}
    check(eq["C2b"].met and not tab["C2b"].met, "C2b inequality conflict is preserved, not silently reconciled")

    aa, cc = appendix_c_o4(98.0)
    check(abs(aa - 0.8825126216053114) < 1e-10 and abs(cc - 1.6182339082203274) < 1e-10, "Appendix-C beta=98 closed-form checkpoint")
    inv = aa*aa + cc*cc + 2*aa*cc*np.cos(np.deg2rad(98.0))
    check(abs(inv - 3.0) < 1e-12, "Appendix-C invariant-length condition")

    c1 = LatticeInput(inp.a_b2, inp.a_b19p, np.sqrt(2.0)*inp.a_b2, inp.c_b19p, inp.beta_deg)
    for branch in (-1, 1):
        g = ips_geometry_from_stretch(stretch_from_lattice(c1).U, branch)
        check(g.rank_one_residual < 1e-12 and g.plane_invariance_residual < 1e-12, f"explicit IPS geometry branch {branch}")

    cosets = correspondence_left_cosets()
    check(len(cosets) == 12 and sum(len(x.elements) for x in cosets) == 48, "left-coset correspondence-variant partition 48 = 12 x 4")

    reg = equation_registry()
    required = {"CMC", "SMC", "PTMC-SC1", "PTMC-SC2", "PTMC-SC3", "NITI-Q", "CT-C2b-TABLE", "APP-C-O4", "IND-HADAMARD", "IND-TWIN-SPECTRAL"}
    check(required <= set(reg), "equation provenance registry coverage")
    print("Equation parity audit: PASS")


if __name__ == "__main__":
    main()
