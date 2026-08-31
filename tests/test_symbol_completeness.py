from __future__ import annotations

from src.provenance import EQUATIONS, equations_for
from src.symbols import EQUATION_SYMBOLS, INPUT_SYMBOLS, OUTPUT_SYMBOLS, definitions_for_equation


def test_every_registered_equation_has_explicit_symbol_definitions():
    keys = {eq.key for eq in EQUATIONS}
    assert keys == set(EQUATION_SYMBOLS), (keys - set(EQUATION_SYMBOLS), set(EQUATION_SYMBOLS) - keys)
    for eq in EQUATIONS:
        defs = definitions_for_equation(eq.key)
        assert defs, eq.key
        for d in defs:
            assert d.symbol.strip()
            assert d.meaning.strip()
            assert d.unit.strip()


def test_every_primary_lattice_input_has_visible_meaning_and_unit():
    required = {"a_b2", "a_b19p", "b_b19p", "c_b19p", "beta_deg", "cc1_tol", "cc2_tol"}
    assert required <= set(INPUT_SYMBOLS)
    for key in required:
        d = INPUT_SYMBOLS[key]
        assert d.meaning and d.unit


def test_core_outputs_have_notation_entries():
    required = {"lambda2", "cmc_distance", "cc2_residual", "epsilon", "detU", "q", "habit_plane", "dA", "s"}
    assert required <= set(OUTPUT_SYMBOLS)


def test_exported_equation_provenance_carries_symbol_dictionary():
    rows = equations_for([eq.key for eq in EQUATIONS])
    assert len(rows) == len(EQUATIONS)
    for row in rows:
        assert row["symbols"]
        assert "operators" in row
        for symbol in row["symbols"] + row["operators"]:
            assert symbol["symbol"]
            assert symbol["meaning"]
            assert symbol["unit"]
