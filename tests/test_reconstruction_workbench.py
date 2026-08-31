from __future__ import annotations

from pathlib import Path

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


def test_ang_parser_preserves_standard_quality_fields_and_nonindexed_phase():
    raw = b"# VERSION: 5\n0.1 0.2 0.3 1 2 100 0.9 1 12 0.5\n0.2 0.3 0.4 2 2 80 -1 0 10 1.0\n"
    out = _read_ang_bytes(raw)
    assert {"IQ", "CI", "phase", "SEM_signal", "fit"}.issubset(out.columns)
    assert out.loc[1, "phase"] == 0
    assert out.attrs["format"] == "ANG"
    assert "radians" in out.attrs["euler_input_unit"]


def test_ang_parser_fails_safe_but_marks_nonstandard_degree_angles():
    raw = b"10 20 30 0 0 100 0.8 1 5 0.5\n"
    out = _read_ang_bytes(raw)
    assert np.isclose(out.loc[0, "phi1_deg"], 10.0)
    assert "inferred" in out.attrs["euler_input_unit"]


def test_ctf_parser_preserves_quality_columns():
    raw = b"Channel Text File\nXCells\t2\nYCells\t1\nPhase\tX\tY\tBands\tError\tEuler1\tEuler2\tEuler3\tMAD\tBC\tBS\n1\t0\t0\t8\t0\t10\t20\t30\t0.5\t100\t20\n"
    out = _read_ctf_bytes(raw)
    assert {"Bands", "Error", "MAD", "BC", "BS"}.issubset(out.columns)
    assert out.attrs["format"] == "CTF"


def test_segmentation_exports_spread_and_boundary_contact_evidence():
    pts = []
    for x, y in [(0, 0), (1, 0), (0, 1), (1, 1)]:
        pts.append({"x": x, "y": y, "phi1_deg": 0.0, "Phi_deg": 20.0, "phi2_deg": 10.0, "phase": 1})
    for x, y in [(2, 0), (3, 0), (2, 1), (3, 1)]:
        pts.append({"x": x, "y": y, "phi1_deg": 30.0, "Phi_deg": 20.0, "phi2_deg": 10.0, "phase": 1})
    grains, adj = _segment_points_to_grains(pd.DataFrame(pts), symmetry_group("cubic"), 5.0, 2, neighbor_radius_factor=1.35)
    assert {"grain_orientation_spread_rms_deg", "grain_orientation_spread_max_deg", "area_est_mapunit2"}.issubset(grains.columns)
    assert {"boundary_contact_count", "boundary_length_est_mapunit"}.issubset(adj.columns)
    assert len(grains) == 2
    assert len(adj) >= 1


def test_workbench_exposes_ct_otsuka_ren_and_clean_four_workspace_navigation():
    source = Path("src/reconstruction_workbench.py").read_text(encoding="utf-8")
    assert "CT + Otsuka–Ren correspondence" in source
    assert "model-derived initial OR" in source
    assert "1 · Reconstruct parent" in source
    assert "2 · NiTi B19′↔B2 cycle" in source
    assert "3 · Forward & batch" in source
    assert "4 · Academic guide" in source
