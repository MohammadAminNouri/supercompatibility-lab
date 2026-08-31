# Validation strategy

Supercompatibility Lab is validated by **numerical reproduction, algebraic cross-checks and software invariants**. The tests are designed to catch coordinate-convention, transpose, normalization and symmetry-enumeration errors—not merely Python exceptions.

## Test command

```bash
pip install -r requirements-dev.txt
pytest -q
```

The repository also compiles all source modules in CI.

## 1. B2→B19′ experimental benchmark

A published binary-NiTi lattice dataset is encoded as a preset:

- \(a_{B2}=3.01\) Å
- \(a_{B19'}=2.898\) Å
- \(b_{B19'}=4.108\) Å
- \(c_{B19'}=4.646\) Å
- \(\beta=97.78^\circ\)

The normalized values must reproduce approximately

\[
a=0.9628,\qquad b=1.3648,\qquad c=1.5435.
\]

The test suite checks these ratios and the resulting matrices.

## 2. C1 A/M compatibility benchmark

The C1 teaching case sets

\[
b=\sqrt2
\]

while retaining the benchmark \(a,c,\beta\) ratios.

The CMC quadratic form must degenerate into two habit-plane solutions close to

\[
m_A^+=(1,-1,2.41966),
\]

\[
m_A^-=(1,-1,-0.31568),
\]

up to the canonical plane scaling used by the code.

The associated non-normalized SMC vectors must reproduce approximately

\[
d_A^+=(0.36938,-0.36938,-0.05378)^T,
\]

\[
d_A^-=(0.115568,-0.115568,0.216812)^T.
\]

These checks are sensitive to the exact correspondence/reciprocal-space convention.

## 3. Twin and metric shear/shear benchmark

The compact benchmark twin subset checks published shear amplitudes and mismatch values. For the C1-compatible example, the best included compound system has a metric shear/shear mismatch near

\[
\varepsilon\approx0.215
\]

and angular mismatch near \(5.9^\circ\).

This is intentionally retained as a **non-zero** benchmark; forcing the metric route to return zero simply because a cofactor condition passes would be a conceptual regression.

## 4. Rounded metric-supercompatible benchmark

A rounded target near

\[
a=0.8825,\quad b=\sqrt2,\quad c=1.6182,\quad \beta=98^\circ
\]

with \(a_{B2}=3.0\) Å must produce a very small metric shear/shear residual for the corresponding twin family. Because the published lattice ratios are rounded, the test expects a small finite numerical residual rather than machine-zero equality.

## 5. PTMC/cofactor cross-check

Tests verify that:

- the stretch tensor is symmetric positive definite;
- its ordered eigenvalues reproduce the known \(\lambda_2\) relation for the built-in system;
- the full CC1–CC3 equations are evaluated from Type-I, Type-II and dedicated Compound rank-one vectors;
- a paired set of perpendicular generators is classified as Compound rather than double-counted as Type-I/Type-II;
- the simplified Type-I/Type-II CC2 conditions are used only for non-compound systems;
- trivial parent axes for which \(QUQ^T=U\) are filtered from domain screening.

## 6. Cross-framework discrepancy test

A dedicated test preserves the fact that a C1-compatible lattice can satisfy the classical cofactor conditions for a conventional domain system while the specialized metric shear/shear benchmark remains nonzero.

This test exists because the two formulations have different demonstrated logical status: the metric/shear construction implies the classical conditions, while complete converse equivalence is not assumed by the software.

## 7. Symmetry invariants

The tests require:

- 48 signed-permutation matrices in the cubic parent point group;
- 4 operations in the correspondence intersection subgroup;
- 7 double-coset classes with sizes 4,4,8,8,8,8,8;
- 9 unique cubic two-fold `<100>/<110>` axes up to sign;
- 16 non-trivial symmetry-labeled two-fold-generated Type-I/Type-II explorer entries;
- the known compound-twin equivalences are recovered as two Type-I/Type-II pairings in the relevant intercorrespondence class;
- 14 non-trivial cofactor domain solutions after compound classification removes double-counting and a trivial symmetry relation;
- 12 unique stretch variants under parent-symmetry conjugation for the built-in benchmark.

These are structural tests, not visual UI checks.

## 8. Temperature evaluator

A synthetic temperature CSV is passed through the same rowwise compatibility engine. Tests check row count, ordering and finite residual outputs.

## 9. Uncertainty engine

With all input standard deviations set to zero, Monte Carlo draws must reproduce the deterministic result exactly within the selected numerical thresholds. This is a basic but powerful regression test for uncertainty plumbing.

For real publication work, nonzero uncertainty should additionally be checked against independent calculations and, if available, a full covariance matrix from diffraction refinement.

## 10. Multi-step metric engine

The generic unit-cell/correspondence stage engine is checked against the specialized B2→B19′ stretch calculation. For the same parent/product metrics and correspondence, the generic and specialized stretch eigenvalues must agree.

## 11. Frontier applicability guard

The tests verify that the built-in monoclinic unique b axis maps to parent `<110>`, so cubic→monoclinic-II extreme-compatible reference tensors are not applied as a certification test. They also verify the current count of commuting stretch-variant pairs for the benchmark.

## 12. ML pipeline

A synthetic exactly linear dataset is used to test:

- target-column validation;
- cross-validation plumbing;
- expected high \(R^2\) for a linear model on linear data;
- candidate lattice prediction;
- post-prediction physics screening.

This validates software behavior, not the scientific accuracy of an arbitrary user dataset.

## 13. What validation does NOT prove

Passing the suite does not prove that:

- a candidate composition can be synthesized;
- exact crystallographic compatibility guarantees low hysteresis or infinite fatigue life;
- an ML model trained on a small/narrow dataset extrapolates correctly;
- independent Gaussian XRD errors are always appropriate;
- frontier compatibility theorems apply outside their stated transformation classes.

Those require separate experimental/theoretical justification in any paper using this software.

## 14. Parent/daughter reconstruction validation

The reconstruction layer is tested independently of the compatibility engine.

The suite verifies:

- Bunge Euler angle → orientation matrix → Bunge Euler round-trip to numerical precision;
- the proper cubic rotation group contains exactly 24 rotations of determinant +1;
- ideal OR variant counts: KS 24, NW 12, Bain 3, Pitsch 12 and Burgers 12;
- an ideal KS daughter orientation generates 24 candidate parent orientations, including the true parent modulo cubic symmetry;
- on a synthetic two-parent KS grain graph with angular noise, **all five reconstruction algorithms** recover two parent clusters with permutation-invariant clustering accuracy ≥99% and mean parent fit below 1°;
- bounded OR refinement reduces the robust neighbor-consistency objective and recovers a deliberately imposed 2° OR perturbation to below 0.4° on the sufficiently constrained synthetic benchmark.

The synthetic benchmark is not evidence that real EBSD reconstruction is always this accurate. Real maps can contain deformation gradients, segmentation errors, retained parent phase, twins, missing variants and intrinsically ambiguous parent solutions.

## 15. Independent compatibility-method validation

Additional tests verify:

- an exact rank-one Hadamard jump gives zero normalized residual, while a full-rank perturbation does not;
- the single-variant spectral routine reproduces the eigenvalues of a known stretch and identifies \(|\lambda_2-1|\);
- the pairwise variant routine reproduces the Ball–James/Bhattacharya \(\chi_1\le1,\chi_2=1,\chi_3\ge1\) criterion;
- for cofactor-passing benchmark domain systems, a direct 101-point volume-fraction scan keeps the laminate middle singular value at one to approximately machine precision (maximum residual below \(10^{-12}\));
- non-trivial algebraically constructed cubic→orthorhombic TC-I and TC-II examples evaluate to residuals below \(10^{-12}\).

These tests are intentionally theorem-level or synthetic checks. They do not replace experimental validation of a proposed alloy.

## 16. Current automated suite

The final academic-reconstruction release collects **96 tests**. They span the original metric/cofactor/twin engine, temperature, uncertainty, inverse/ML/multi-step utilities, five parent-reconstruction families, robust ANG/CTF parsing, raw-point segmentation, reconstruction reporting, cross-method ARI/NMI/boundary/orientation agreement, academic ZIP export, the metric-aware B19′→B2→B19′ cycle/retransformation workflow, the CT/Otsuka–Ren model-derived starting-OR bridge, cleaner four-workspace reconstruction navigation, and the independent compatibility methods. The complete suite was validated in deterministic batches because the execution harness times out long single test commands; every collected test passed.


## 17. Analytical equation-parity validation

The final release includes an independent deterministic audit in `scripts/verify_equation_parity.py` plus pytest coverage. The release requires:

- analytical B2/B19′ CMC matrix = general metric/correspondence CMC to <1e-13 Frobenius residual;
- analytical CMC eigenvalues = numerical eigenspectrum to <1e-13;
- explicit quadratic polynomial = `4 u^T CMC u`;
- C1/C2a/C2b/C3 constructed family checkpoints;
- D1/D2/E higher-degeneracy checkpoints;
- Appendix-C beta=98° closed-form values and invariant-length relation;
- explicit PTMC IPS reconstruction with rank-one and plane-invariance residuals <1e-12 on both branches;
- left-coset partition `48 = 12 x 4`;
- equation-provenance registry coverage.

## 18. Source-discrepancy regression tests

Three source inconsistencies are deliberately testable. The suite fails if the software silently erases them:

1. analytical CMC `q3` uses the algebraically required plus branch and must match the printed CMC matrix eigenspectrum;
2. C2b must expose both the analytical `b<=sqrt(2)` branch and the conflicting later table `b>=sqrt(2)` branch;
3. the later `0.269706` distance value and the visually/parsed squared expression must remain separately reproducible.

See `docs/SOURCE_DISCREPANCIES.md`.

## 19. Paper-ready record validation

Automated tests require the reproducibility record to remain JSON-readable and to include source DOI, build ID, equation provenance, claim boundaries and a SHA-256 fingerprint. The Markdown audit export is also tested for the presence of equation-provenance and claim-boundary sections.

## 20. Release-gate commands

The release was locally validated with the deterministic scientific audits plus the complete 89-test suite in bounded batches:

```bash
python -m pytest --collect-only -q
python -m pytest -q tests/test_compatibility_methods.py tests/test_core.py tests/test_design.py tests/test_equation_parity.py tests/test_ml.py tests/test_ptmc.py
python -m pytest -q tests/test_release_structure.py tests/test_reproducibility.py tests/test_research_modules.py tests/test_symbol_completeness.py tests/test_symmetry.py tests/test_twins.py
python -m pytest -q tests/test_reconstruction.py
python -m pytest -q tests/test_reconstruction_reporting.py
python -m pytest -q tests/test_reconstruction_workbench.py
python -m pytest -q tests/test_reconstruction_academic.py
python -m pytest -q tests/test_retransformation.py
python scripts/deployment_preflight.py
python scripts/verify_equation_parity.py
python scripts/verify_release.py
python scripts/verify_embedded_app.py
python research_selftest.py
```

`pytest --collect-only` reports **96 tests**. The final validated batches cover all 96 tests; the reconstruction file is split into two invocations only to stay below the execution harness limit. Every collected test passed. A single long `pytest -q` invocation exceeds the execution harness time limit, so the release audit does not falsely claim that one monolithic command completed here. The deterministic release audit, embedded-bundle byte-integrity audit and embedded-engine scientific self-test also passed. Streamlit itself is health-checked in GitHub Actions after dependencies are installed; the local packaging environment used for this release did not have Streamlit installed, so no local runtime-health result is claimed.


## 21. Notation/symbol completeness validation

The final release treats notation as part of scientific reproducibility. Automated tests require every registered equation to have explicit symbol definitions, units/scales and operator metadata; every primary lattice input to have a visible meaning and unit; every core output to have a notation entry; and every equation-provenance export to carry its symbol/operator dictionary. The Streamlit equation helper renders those definitions immediately beside the displayed relation.

The UI additionally exposes a plain-language acronym/crystallographic-notation legend and labels plot axes with physical names and units rather than raw internal column names wherever a technical quantity is plotted.


## Transformation-matched synthetic parent reconstruction validation

The built-in two-parent validation dataset is now generated from the **exact OR matrix and parent/daughter symmetry groups currently selected in the UI**. This closes a validation loophole in which a KS FCC→BCC synthetic dataset could previously be reconstructed under a different transformation such as NiTi B2→B19′. The numerical reconstruction algorithms themselves were not changed.

Known-truth validation no longer uses majority-remapped cluster accuracy as a primary metric because that quantity can report 100% for pathological singleton over-segmentation. The academic workbench instead reports: **Adjusted Rand Index (ARI), homogeneity, completeness, V-measure, parent-count error, singleton-parent fraction, fragmentation per true parent, and adjacency-boundary precision/recall/F1/Jaccard/false-boundary rate**.

A dedicated regression test encodes the former failure mode: two true parents split into 12 singleton reconstructed parents must produce truth ARI = 0, poor completeness, boundary precision = 1/15 and boundary F1 = 0.125 rather than being presented as a correct reconstruction. A transformation-matched CT/Otsuka–Ren NiTi synthetic benchmark with 6 daughters per parent and 0.35° noise recovers exactly two parents with truth ARI = 1, completeness = 1, homogeneity = 1 and boundary F1 = 1 for the tested graph, variant-graph and operator/groupoid routes.

Cross-method ARI/NMI remain in the software, but are explicitly labelled **agreement, not accuracy**. Near-zero OR residual and high heuristic support are also explicitly guarded against when most reconstructed parents are singletons, because a singleton can trivially fit its own best candidate.
