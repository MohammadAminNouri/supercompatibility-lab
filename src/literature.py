from __future__ import annotations

from pathlib import Path

import pandas as pd


DATA_PATH = Path(__file__).resolve().parents[1] / "data" / "literature.csv"


def load_literature() -> pd.DataFrame:
    df = pd.read_csv(DATA_PATH)
    # Normalize whitespace in hand-curated text fields without converting empty
    # numeric cells into misleading zeros.
    for c in ["record_type", "alloy", "composition_atpct", "parent_phase", "martensite_phase", "key_result", "doi", "title", "evidence_note"]:
        if c in df:
            df[c] = df[c].fillna("").astype(str).str.strip()
    return df


def literature_by_type(record_type: str | None = None) -> pd.DataFrame:
    df = load_literature()
    if record_type and record_type != "All":
        return df[df["record_type"] == record_type].reset_index(drop=True)
    return df
