import numpy as np

from src.core import (
    LatticeInput,
    c1_status,
    cmc_degeneracy,
    cmc_matrix,
    normalized_metrics,
    smc_matrix,
)


def test_published_binary_ratios():
    inp = LatticeInput(3.01, 2.898, 4.108, 4.646, 97.78)
    a, b, c = inp.ratios()
    assert np.isclose(a, 0.9627906977, atol=1e-10)
    assert np.isclose(b, 1.3647840532, atol=1e-10)
    assert np.isclose(c, 1.5435215947, atol=1e-10)


def test_c1_distance_for_binary_example():
    inp = LatticeInput(3.01, 2.898, 4.108, 4.646, 97.78)
    c1 = c1_status(inp)
    assert np.isclose(c1["equality_residual"], -0.137363, atol=2e-5)
    assert c1["inequality_margin"] > 0


def test_c1_habit_planes_and_smc_vectors():
    a_b2 = 3.01
    inp = LatticeInput(
        a_b2,
        0.9628 * a_b2,
        np.sqrt(2) * a_b2,
        1.5435 * a_b2,
        97.78,
    )
    ma, mm = normalized_metrics(inp)
    deg = cmc_degeneracy(cmc_matrix(ma, mm))
    assert deg.order == 1
    assert len(deg.habit_planes) == 2

    # Order-independent comparison against the two benchmark planes.
    targets = [
        np.array([1.0, -1.0, 2.41966]),
        np.array([1.0, -1.0, -0.31568]),
    ]
    for target in targets:
        assert any(np.allclose(p, target, atol=3e-4) for p in deg.habit_planes)

    smc = smc_matrix(ma, mm)
    expected = {
        2.41966: np.array([0.36938, -0.36938, -0.05378]),
        -0.31568: np.array([0.115568, -0.115568, 0.216812]),
    }
    for p in deg.habit_planes:
        key = min(expected, key=lambda k: abs(p[2] - k))
        assert np.allclose(smc @ p, expected[key], atol=8e-5)
