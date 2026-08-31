from __future__ import annotations

import json

from src.core import LatticeInput
from src.reproducibility import paper_ready_record, record_json, record_markdown, record_sha256


def test_paper_ready_record_is_machine_readable_and_contains_claim_boundaries():
    inp = LatticeInput(3.01, 2.898, 4.108, 4.646, 97.78)
    rec = paper_ready_record(inp)
    text = record_json(rec)
    loaded = json.loads(text)
    assert loaded["schema"] == "supercompatibility-lab-paper-record-v1"
    assert loaded["provenance"]["primary_source"]["doi"] == "10.1016/j.actamat.2026.122399"
    assert len(loaded["equations_used"]) >= 15
    assert len(loaded["claim_boundaries"]) >= 4
    assert len(record_sha256(rec)) == 64
    md = record_markdown(rec)
    assert "Equation provenance" in md and "Claim boundaries" in md
