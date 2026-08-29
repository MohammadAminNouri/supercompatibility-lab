import numpy as np

from src.core import LatticeInput, normalized_metrics
from src.ptmc import stretch_from_lattice
from src.symmetry import (
    compound_twin_pairs,
    correspondence_intersection_subgroup,
    cubic_point_group,
    double_cosets,
    full_twin_explorer,
    parent_twofold_axes,
    stretch_variants,
)


def test_group_and_double_coset_counts():
    assert len(cubic_point_group()) == 48
    assert len(correspondence_intersection_subgroup()) == 4
    dcs = double_cosets()
    assert len(dcs) == 7
    assert sum(len(d.matrices) for d in dcs) == 48
    assert sorted(len(d.matrices) for d in dcs) == [4,4,8,8,8,8,8]


def test_complete_twofold_axis_and_twin_explorer_counts():
    assert len(parent_twofold_axes()) == 9
    inp = LatticeInput(3.01, 2.898, 4.108, 4.646, 97.78)
    ma, mm = normalized_metrics(inp)
    twins = full_twin_explorer(ma, mm)
    assert len(twins) == 16
    assert {t.twin.twin_type for t in twins} == {"Type I", "Type II"}


def test_twelve_stretch_variants_for_builtin_transformation():
    inp = LatticeInput(3.01, 2.898, 4.108, 4.646, 97.78)
    U = stretch_from_lattice(inp).U
    assert len(stretch_variants(U)) == 12


def test_compound_twins_are_identified_in_c1_benchmark():
    inp = LatticeInput(3.01, 0.9628*3.01, np.sqrt(2.0)*3.01, 1.5435*3.01, 97.78)
    ma, mm = normalized_metrics(inp)
    twins = full_twin_explorer(ma, mm)
    pairs = compound_twin_pairs(twins)
    assert len(pairs) == 2
    assert len({twins[i].double_coset for i, _ in pairs}) == 1
    for i, j in pairs:
        assert twins[i].twin.twin_type == "Type I"
        assert twins[j].twin.twin_type == "Type II"
