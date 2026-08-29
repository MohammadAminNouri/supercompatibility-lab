from __future__ import annotations

from .core import LatticeInput


PRESETS: dict[str, LatticeInput | None] = {
    "Published binary NiTi example": LatticeInput(
        a_b2=3.01,
        a_b19p=2.898,
        b_b19p=4.108,
        c_b19p=4.646,
        beta_deg=97.78,
    ),
    "C1-compatible teaching example": LatticeInput(
        a_b2=3.01,
        a_b19p=0.9628 * 3.01,
        b_b19p=(2.0**0.5) * 3.01,
        c_b19p=1.5435 * 3.01,
        beta_deg=97.78,
    ),
    "Rounded supercompatible target": LatticeInput(
        a_b2=3.0,
        a_b19p=0.8825 * 3.0,
        b_b19p=(2.0**0.5) * 3.0,
        c_b19p=1.6182 * 3.0,
        beta_deg=98.0,
    ),
    "Custom": None,
}
