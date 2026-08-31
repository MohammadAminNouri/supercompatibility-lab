from __future__ import annotations

import sys
import types

import numpy as np
import pandas as pd

try:
    import streamlit  # noqa: F401
except ModuleNotFoundError:
    sys.modules["streamlit"] = types.ModuleType("streamlit")

from src.reconstruction import symmetry_group
from src.reconstruction_workbench import (
    _coerce_grain_columns,
    _parse_matrix,
    _parse_vec,
    _quality_label,
    _read_ang_bytes,
    _read_ctf_bytes,
    _segment_points_to_grains,
)


def test_manual_parsers_and_quality_labels():
    assert _quality_label(0.5, 0.9) == "Very strong"
    assert _quality_label(6.0, 0.2) == "Weak / ambiguous"
    np.testing.assert_allclose(_parse_vec("1 -1 0", "v"), [1, -1, 0])
    np.testing.assert_allclose(_parse_matrix("1 0 0;0 1 0;0 0 1"), np.eye(3), atol=1e-12)


def test_common_grain_column_aliases_are_normalized():
    df = _coerce_grain_columns(pd.DataFrame({"grain": [1], "phi1": [10.0], "Phi": [20.0], "phi2": [30.0]}))
    assert {"grain_id", "phi1_deg", "Phi_deg", "phi2_deg"}.issubset(df.columns)


def test_ang_parser_converts_radians_to_degrees():
    raw = b"# header\n0.1 0.2 0.3 1 2 100 0.9 1 0 0\n"
    out = _read_ang_bytes(raw)
    assert np.isclose(out.loc[0, "phi1_deg"], np.degrees(0.1))
    assert np.isclose(out.loc[0, "Phi_deg"], np.degrees(0.2))
    assert out.loc[0, "phase"] == 1


def test_ctf_parser_reads_standard_euler_columns():
    raw = b"Channel Text File\nPhase X Y Bands Error Euler1 Euler2 Euler3 MAD BC BS\n1 0 0 8 0 10 20 30 0.5 100 20\n"
    out = _read_ctf_bytes(raw)
    assert out.loc[0, "phi1_deg"] == 10
    assert out.loc[0, "Phi_deg"] == 20
    assert out.loc[0, "phi2_deg"] == 30


def test_explicit_point_segmentation_separates_misoriented_regions():
    pts = []
    for x, y in [(0, 0), (1, 0), (0, 1), (1, 1)]:
        pts.append({"x": x, "y": y, "phi1_deg": 0.0, "Phi_deg": 20.0, "phi2_deg": 10.0, "phase": 1})
    for x, y in [(3, 0), (4, 0), (3, 1), (4, 1)]:
        pts.append({"x": x, "y": y, "phi1_deg": 30.0, "Phi_deg": 20.0, "phi2_deg": 10.0, "phase": 1})
    grains, _ = _segment_points_to_grains(pd.DataFrame(pts), symmetry_group("cubic"), 5.0, 2)
    assert len(grains) == 2
    assert sorted(grains["point_count"].tolist()) == [4, 4]
