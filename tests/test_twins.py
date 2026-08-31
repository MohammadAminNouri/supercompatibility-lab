import numpy as np

from src.core import LatticeInput, rank_matches


def c1_input():
    a_b2 = 3.01
    return LatticeInput(
        a_b2,
        0.9628 * a_b2,
        np.sqrt(2) * a_b2,
        1.5435 * a_b2,
        97.78,
    )


def best_for(label: str):
    matches = [m for m in rank_matches(c1_input()) if m.twin.label == label]
    return min(matches, key=lambda m: m.epsilon)


def test_compound_type_i_benchmark():
    m = best_for("Compound / Type-I equivalent")
    assert np.isclose(m.twin.shear_amplitude, 0.27325, atol=8e-5)
    assert np.isclose(m.epsilon, 0.215, atol=0.002)
    assert np.isclose(m.angle_deg, 5.9, atol=0.15)


def test_compound_type_ii_benchmark():
    m = best_for("Compound / Type-II equivalent")
    assert np.isclose(m.twin.shear_amplitude, 0.27325, atol=8e-5)
    assert np.isclose(m.epsilon, 0.809, atol=0.003)
    assert np.isclose(m.angle_deg, 37.0, atol=0.15)


def test_type_i_family_a_benchmark():
    m = best_for("Type I family A")
    assert np.isclose(m.twin.shear_amplitude, 0.30615, atol=8e-5)
    assert np.isclose(m.epsilon, 0.453, atol=0.003)
    assert np.isclose(m.angle_deg, 25.3, atol=0.15)


def test_type_ii_family_a_benchmark():
    m = best_for("Type II family A")
    assert np.isclose(m.twin.shear_amplitude, 0.30615, atol=8e-5)
    assert np.isclose(m.epsilon, 0.566, atol=0.003)
    assert np.isclose(m.angle_deg, 30.0, atol=0.15)


def test_type_i_family_b_benchmark():
    m = best_for("Type I family B")
    assert np.isclose(m.twin.shear_amplitude, 0.25502, atol=8e-5)
    assert np.isclose(m.epsilon, 0.469, atol=0.003)
    assert np.isclose(m.angle_deg, 26.9, atol=0.15)


def test_type_ii_family_b_benchmark():
    m = best_for("Type II family B")
    assert np.isclose(m.twin.shear_amplitude, 0.25502, atol=8e-5)
    assert np.isclose(m.epsilon, 0.906, atol=0.003)
    assert np.isclose(m.angle_deg, 64.8, atol=0.15)


def test_rounded_supercompatible_target():
    inp = LatticeInput(
        3.0,
        0.8825 * 3.0,
        np.sqrt(2) * 3.0,
        1.6182 * 3.0,
        98.0,
    )
    matches = [m for m in rank_matches(inp) if m.twin.label == "Type I family A"]
    best = min(matches, key=lambda m: m.epsilon)
    assert best.epsilon < 5e-4
    assert np.allclose(best.habit_plane, [1, -1, 2], atol=5e-4)
