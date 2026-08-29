from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class MicrostructureContext:
    mean_grain_size_nm: float | None = None
    grain_size_distribution: str = "Not specified"
    precipitate_type: str = "Not specified"
    precipitate_size_nm: float | None = None
    precipitate_fraction_pct: float | None = None
    retained_martensite_pct: float | None = None
    dislocation_density_m2: float | None = None
    processing: str = "Not specified"


def evidence_annotations(ctx: MicrostructureContext) -> list[dict[str, str]]:
    """Evidence-linked context notes; deliberately not a fatigue predictor."""
    notes: list[dict[str, str]] = []
    if ctx.mean_grain_size_nm is not None:
        notes.append({
            "factor": "Grain size",
            "note": "Grain size can change transformation pathways, strength, hysteresis and cyclic stability. Recent NiTiCu/NiTiCuCo studies show that both average size and spatial distribution can matter.",
            "scope": "Context only — no universal quantitative law is assumed.",
        })
    if ctx.grain_size_distribution not in {"Not specified", "Uniform"}:
        notes.append({
            "factor": "Grain-size distribution",
            "note": "Bimodal and gradient distributions are an active design variable; simulation work reports lower dissipation for some distributions than for uniform nanograins.",
            "scope": "Mechanism evidence, not a direct experimental-performance prediction.",
        })
    if ctx.precipitate_type not in {"", "Not specified", "None"}:
        notes.append({
            "factor": "Precipitates",
            "note": "Fine/coherent precipitates can strongly affect repeatability. Ti2Cu-containing TiNiCu films are a benchmark example of precipitate-assisted ultralow functional fatigue.",
            "scope": "Effect depends on phase, coherency, size, spacing and processing.",
        })
    if ctx.retained_martensite_pct is not None and ctx.retained_martensite_pct > 0:
        notes.append({
            "factor": "Retained martensite",
            "note": "Retained martensite is a documented contributor to functional/elastocaloric degradation during cycling in superelastic NiTi.",
            "scope": "Measured retained fraction is metadata; this app does not infer damage kinetics from it.",
        })
    if ctx.dislocation_density_m2 is not None and ctx.dislocation_density_m2 > 0:
        notes.append({
            "factor": "Dislocations",
            "note": "Transformation-induced plasticity and dislocation accumulation are documented fatigue mechanisms in NiTi-based systems.",
            "scope": "No universal conversion from dislocation density to fatigue life is imposed.",
        })
    if not notes:
        notes.append({
            "factor": "No microstructure metadata",
            "note": "Crystallographic compatibility is only one layer of functional performance. Add grain/precipitate/fatigue-state metadata when available.",
            "scope": "The geometric calculations remain valid independently of this optional metadata.",
        })
    return notes
