# Final release audit

Release audit date: 2026-08-31

Build identifier: `2026-08-31-final-truth-valid-reconstruction`

Recommended Streamlit entrypoint: `supercompatibility_final.py`

## Automated scientific/software tests

Collection check:

```bash
python -m pytest --collect-only -q
```

Result: **96 tests collected**. Because the execution harness imposes a time limit on long single commands, the complete suite was executed in deterministic batches. All **96/96 tests passed**; the reconstruction module was split into a fast 8-test batch plus its one heavier multi-preset test.

Coverage includes the metric/CMC/SMC engine, PTMC/cofactor cross-checks, analytical source-equation parity, Type-I/Type-II/Compound domain classification, symmetry/cosets/double cosets, independent compatibility methods, temperature evaluation, uncertainty propagation, inverse design, ML screening, multi-step transformations, five parent/daughter reconstruction families, OR refinement, robust ANG/CTF import, raw-point segmentation, academic reconstruction comparison/export, metric-aware NiTi B2↔B19′ cycle regeneration and later-cycle matching, reproducibility export and notation/symbol completeness.

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

`python scripts/verify_embedded_app.py` verifies that all 27 embedded Python modules match `src/` byte-for-byte, all 9 embedded data assets match `data/` byte-for-byte, and `app.py`, `streamlit_app.py`, `supercompatibility_r7.py` and `supercompatibility_final.py` are byte-identical.

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


## CT/Otsuka–Ren reconstruction addition

The parent/daughter workbench preserves all previous ORs and reconstruction methods and adds a dedicated NiTi CT/Otsuka–Ren route. The source-derived B2/B19′ correspondence and metrics are kept distinct from the software bridge that obtains a model-derived starting OR by polar decomposition. The UI exposes this caveat, the exact correspondence matrices and the resulting rotation/stretches. Top-level navigation is grouped into four academic workspaces without removing any prior workflow.


### Known-truth reconstruction hardening

- Selected-OR/symmetry-matched synthetic generation: PASS.
- CT/Otsuka–Ren two-parent synthetic recovery: PASS (truth ARI/completeness/homogeneity/boundary F1 = 1.000 in dedicated regression tests).
- Singleton over-segmentation counterexample: PASS (truth ARI = 0 and boundary F1 = 0.125; no false 100% academic accuracy).
- Cross-method agreement-vs-accuracy warning: present in UI and exports.
- Academic evidence ZIP includes `validation/known_truth_metrics.csv` when truth labels are available.
