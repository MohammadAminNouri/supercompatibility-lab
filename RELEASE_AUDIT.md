# Final release audit

Release audit date: 2026-08-31

Build identifier: `2026-08-31-final-notation-explicit`

Recommended Streamlit entrypoint: `supercompatibility_final.py`

## Automated scientific/software tests

Command:

```bash
python -m pytest -q
```

Result:

```text
63 passed
```

Coverage includes the metric/CMC/SMC engine, PTMC/cofactor cross-checks, analytical source-equation parity, Type-I/Type-II/Compound domain classification, symmetry/cosets/double cosets, independent compatibility methods, temperature evaluation, uncertainty propagation, inverse design, ML screening, multi-step transformations, parent/daughter reconstruction, OR refinement, reproducibility export and notation/symbol completeness.

## Equation/source parity audit

`python scripts/verify_equation_parity.py` passes the closed-form CMC/general-metric equivalence, analytical/numerical eigenspectrum, source distance checkpoints, preserved source discrepancies, Appendix-C analytical checkpoint, explicit IPS branches, left-coset partition and equation-provenance registry coverage.

## Independent numerical benchmark audit

`python scripts/verify_release.py` reproduces:

- binary normalized ratios `(0.9627907, 1.3647841, 1.5435216)`;
- middle principal stretch `lambda2 = 0.9650480588`;
- 12 symmetry-distinct stretch variants;
- 42/66 pairwise rank-one-compatible variant pairs;
- all-volume-fraction middle-singular-value residual `8.882e-16` for the selected cofactor benchmark;
- two-parent synthetic reconstruction with 100% clustering accuracy and mean parent fit about 0.174 degrees for all five implemented reconstruction families.

## Embedded deployment integrity

`python scripts/verify_embedded_app.py` verifies that all 23 embedded Python modules match `src/` byte-for-byte, all 9 embedded data assets match `data/` byte-for-byte, and `app.py`, `streamlit_app.py`, `supercompatibility_r7.py` and `supercompatibility_final.py` are byte-identical.

## Scientific self-test

`python research_selftest.py` passes using the extracted embedded engine, independently reproducing the main numerical/reconstruction checkpoints.

## Notation completeness

The final build uses a no-symbol-left-behind rule. Every registered equation has a local symbol dictionary plus detected operator definitions. Primary physical inputs show full names, symbols and units directly beside the input. Core outputs and technical tables/plots carry plain-language explanations, and the paper-ready JSON/Markdown/CSV exports preserve the same notation dictionary. `tests/test_symbol_completeness.py` is part of the release gate.

## Compilation and deployment preflight

All source, scripts, tests and four Streamlit entrypoints compile. `python scripts/deployment_preflight.py` passes and confirms the final build identifier and required repository structure.

## Streamlit runtime note

The local packaging environment does not contain the `streamlit` package, so no local browser/server-health claim is made. GitHub Actions installs the declared dependencies, starts `supercompatibility_final.py`, and requires the Streamlit `/_stcore/health` endpoint to return healthy before CI can pass.

## Scientific scope

Passing these checks establishes internal numerical consistency and reproduction of the encoded benchmark/equation cases. It does not prove synthesizability of a proposed composition, generalization of a user-trained ML model, or that crystallographic compatibility alone determines hysteresis, fatigue life or functional performance. Such claims require independent experimental evidence and appropriate validation.
