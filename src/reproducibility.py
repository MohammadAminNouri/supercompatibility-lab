from __future__ import annotations

import hashlib
import json
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from .core import LatticeInput, cmc_degeneracy, cmc_matrix, normalized_metrics, smc_matrix
from .distances import compatibility_dashboard
from .equation_engine import (
    appendix_c_o4,
    first_order_families,
    higher_order_degeneracy_families,
    niti_cmc_from_input,
    paper_defined_distances,
    verify_analytic_cmc_against_general,
)
from .provenance import BUILD_ID, equations_for, method_provenance
from .ptmc import stretch_from_lattice


def _jsonable(x: Any) -> Any:
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, np.generic):
        return x.item()
    if hasattr(x, "__dataclass_fields__"):
        return {k: _jsonable(v) for k, v in asdict(x).items()}
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (tuple, list)):
        return [_jsonable(v) for v in x]
    return x


def paper_ready_record(
    inp: LatticeInput,
    *,
    cc1_tol: float = 1e-5,
    cc2_tol: float = 1e-5,
    equation_keys: Iterable[str] | None = None,
    random_seed: int = 20260831,
) -> dict[str, Any]:
    ma, mm = normalized_metrics(inp)
    cmc = cmc_matrix(ma, mm)
    deg = cmc_degeneracy(cmc)
    smc = smc_matrix(ma, mm)
    stretch = stretch_from_lattice(inp)
    dash = compatibility_dashboard(inp, cc1_tol=cc1_tol, cc2_tol=cc2_tol)
    a, b, c = inp.ratios()
    if equation_keys is None:
        equation_keys = [
            "PTMC-SC1", "PTMC-SC2", "PTMC-SC3", "CMC", "CMC-Q", "CMC-DEGEN",
            "SMC", "SMC-D", "SHEAR-SHEAR", "EPSILON", "NITI-CMC", "NITI-Q",
            "CT-C1", "CT-C2a", "CT-C2b-EQ", "CT-C3", "CT-D1", "CT-D2", "CT-E",
        ]
    record: dict[str, Any] = {
        "schema": "supercompatibility-lab-paper-record-v1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "build_id": BUILD_ID,
        "random_seed": int(random_seed),
        "provenance": method_provenance(cc1_tol, cc2_tol),
        "inputs": {
            "a_B2": {"value": inp.a_b2, "unit": "angstrom"},
            "a_B19p": {"value": inp.a_b19p, "unit": "angstrom"},
            "b_B19p": {"value": inp.b_b19p, "unit": "angstrom"},
            "c_B19p": {"value": inp.c_b19p, "unit": "angstrom"},
            "beta": {"value": inp.beta_deg, "unit": "degree"},
            "normalized_ratios": {"a": a, "b": b, "c": c, "unit": "dimensionless"},
        },
        "equations_used": equations_for(equation_keys),
        "matrices": {"M_A": ma, "M_M": mm, "CMC": cmc, "SMC": smc, "U": stretch.U},
        "core_results": {
            "lambda": stretch.eigenvalues,
            "det_U": stretch.determinant,
            "CMC_eigenvalues": deg.eigenvalues,
            "CMC_degeneracy_order": deg.order,
            "CMC_exact_at_tolerance": deg.exact,
            "CMC_habit_planes": deg.habit_planes,
            "dashboard": dash,
        },
        "analytical_cross_checks": {
            "analytic_CMC": niti_cmc_from_input(inp),
            "analytic_vs_general": verify_analytic_cmc_against_general(inp),
            "first_order_families_equation_interpretation": first_order_families(a, b, c, inp.beta_deg, c2b_interpretation="equation"),
            "first_order_families_table_interpretation": first_order_families(a, b, c, inp.beta_deg, c2b_interpretation="table"),
            "higher_order_degeneracy": higher_order_degeneracy_families(a, b, c, inp.beta_deg),
            "paper_defined_distances": paper_defined_distances(a, b, c, inp.beta_deg),
        },
        "claim_boundaries": [
            "Numerical tolerances define pass/fail labels and must be reported with results.",
            "A/M compatibility, classical cofactor conditions and metric shear/shear intercompatibility are reported separately.",
            "A source discrepancy in the C2b inequality is retained rather than silently reconciled.",
            "Literature/microstructure/performance metadata are contextual evidence, not outputs of the crystallographic equations.",
            "Machine-learning predictions are candidate-screening outputs and require independent physical/experimental validation.",
        ],
    }
    return _jsonable(record)


def record_json(record: dict[str, Any]) -> str:
    return json.dumps(record, indent=2, sort_keys=True, allow_nan=False)


def record_sha256(record: dict[str, Any]) -> str:
    return hashlib.sha256(record_json(record).encode("utf-8")).hexdigest()


def record_markdown(record: dict[str, Any]) -> str:
    i = record["inputs"]
    r = record["core_results"]
    p = record["provenance"]
    lines = [
        "# Reproducible crystallographic compatibility record",
        "",
        f"- Build: `{record['build_id']}`",
        f"- Generated UTC: `{record['generated_utc']}`",
        f"- Primary source DOI: `{p['primary_source']['doi']}`",
        f"- Random seed: `{record['random_seed']}`",
        "",
        "## Physical inputs",
        f"- B2 austenite lattice parameter a_B2 = {i['a_B2']['value']} Å",
        f"- B19′ martensite a = {i['a_B19p']['value']} Å",
        f"- B19′ martensite b = {i['b_B19p']['value']} Å",
        f"- B19′ martensite c = {i['c_B19p']['value']} Å",
        f"- Monoclinic β = {i['beta']['value']}°",
        "",
        "## Core numerical results",
        f"- Principal stretches λ = {r['lambda']}",
        f"- det(U) = {r['det_U']}",
        f"- CMC eigenvalues = {r['CMC_eigenvalues']}",
        f"- CMC degeneracy order = {r['CMC_degeneracy_order']}",
        f"- Exact at configured tolerance = {r['CMC_exact_at_tolerance']}",
        "",
        "## Equation provenance",
    ]
    for eq in record["equations_used"]:
        lines.append(f"- **{eq['key']}** — `{eq['source_location']}` — {eq['meaning']}")
    lines += ["", "## Claim boundaries"]
    for x in record["claim_boundaries"]:
        lines.append(f"- {x}")
    lines += ["", f"Record SHA-256: `{record_sha256(record)}`"]
    return "\n".join(lines) + "\n"
