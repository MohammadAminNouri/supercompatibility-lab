from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class FunctionalMetrics:
    thermal_hysteresis_C: float | None = None
    stress_hysteresis_MPa: float | None = None
    deltaTad_K: float | None = None
    cycle_life: float | None = None
    transformation_temp_shift_C: float | None = None
    recoverable_strain_pct: float | None = None
    material_COP: float | None = None
    Ms_C: float | None = None
    Mf_C: float | None = None
    As_C: float | None = None
    Af_C: float | None = None


def derived_metrics(m: FunctionalMetrics) -> dict[str, float | None]:
    midpoint_hysteresis = None
    if None not in (m.Ms_C, m.Mf_C, m.As_C, m.Af_C):
        cooling_mid = (float(m.Ms_C) + float(m.Mf_C)) / 2.0
        heating_mid = (float(m.As_C) + float(m.Af_C)) / 2.0
        midpoint_hysteresis = heating_mid - cooling_mid
    window = None
    if m.Ms_C is not None and m.Af_C is not None:
        window = float(m.Af_C) - float(m.Ms_C)
    return {
        "midpoint_thermal_hysteresis_C": midpoint_hysteresis,
        "Af_minus_Ms_C": window,
    }
