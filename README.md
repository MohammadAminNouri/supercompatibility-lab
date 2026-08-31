# Supercompatibility Lab

A responsive **Streamlit research platform for martensitic compatibility, cofactor conditions, transformation twins, analytical equation parity, uncertainty, temperature dependence, inverse design and paper-ready reproducibility**.

The built-in reference transformation is **B2 austenite → B19′ martensite**. The software keeps several theoretical questions separate instead of collapsing them into one vague “compatibility score”:

1. **Single-variant austenite/martensite compatibility** — stretch condition \(\lambda_2=1\) and CMC degeneracy.
2. **Classical cofactor supercompatibility** — full CC1–CC3 evaluation after classifying non-trivial systems as Type-I, Type-II or Compound domains.
3. **Metric/correspondence intercompatibility** — CMC habit planes, SMC shear and the shear/shear mismatch \(\varepsilon\).
4. **Transformation-twin structure** — symmetry reduction, double cosets and full two-fold-generated twin exploration for the built-in transformation.
5. **Research extensions** — temperature sweeps, XRD uncertainty propagation, literature benchmarks, microstructure metadata, elastocaloric/fatigue metrics, inverse design, user-trained composition→lattice ML, multi-step transformation chains and guarded frontier diagnostics.

The project is designed so that **physics calculations live in `src/` and the Streamlit UI only orchestrates them**. The final release additionally binds every major result to an equation/method provenance record: displayed relation, source or extension class, source location, implementation path, scope/caveat, tolerance and reproducible export. Numerical claims are backed by automated tests, analytical/general parity checks and explicit tolerances.

---

## Why this is more than a calculator

The app is structured as a reproducible research workflow:

```text
measured lattice parameters
        ↓
metrics + correspondence
        ↓
┌─────────────────────┬────────────────────────┐
│ metric route        │ classical stretch route│
│ CMC → habit plane   │ U → λ1, λ2, λ3         │
│ SMC → A/M shear     │ CC1, CC2, CC3          │
└─────────────────────┴────────────────────────┘
        ↓                       ↓
full symmetry/twin exploration + cross-checking
        ↓
explicit residuals / uncertainty / T-dependence
        ↓
inverse lattice design and data-driven screening
```

The program **does not** assign an arbitrary “95% compatible” score. It reports physical residuals and inequality margins with units/definitions.

---

## UX design

The Streamlit interface has two layers:

- **Guided overview** — default learning/screening path that answers four questions in plain language before exposing residuals: single-variant A/M fit, CMC habit-plane fit, cofactor domain fit, and metric shear/shear fit.
- **Research workspaces** — full matrices, twin/domain tables, temperature/uncertainty, reconstruction, inverse design, ML and theorem-level cross-checks.

The B2/B19′ lattice input form is hidden automatically in the parent-reconstruction workspace because it is irrelevant there. Physical input labels use full phase/parameter names and display units; dimensionless graph or ratio controls are explicitly identified. Reconstruction also asks for the orientation-matrix direction convention instead of silently assuming that all EBSD exports use the same convention.

**No-symbol-left-behind rule:** every mathematical symbol shown in the app is defined immediately beside the formula, input, output, or technical table where it appears. The same notation dictionary is embedded in paper-ready JSON/Markdown exports, so a result remains interpretable outside the UI.

**Final validation:** 93/93 automated tests pass across the complete suite when run in validated batches, including dedicated EBSD import/segmentation, academic reconstruction comparison and export tests. Equation/source parity, numerical release, embedded-bundle integrity, scientific self-test and deployment-preflight audits also pass. GitHub Actions performs the actual Streamlit server health check after installing dependencies.

**B19′ → B2 → B19′ cycle reconstruction:** the parent/daughter workspace can now reuse a completed B19′→B2 reconstruction and regenerate every symmetry-distinct B19′ orientation branch from each reconstructed B2 parent. It reports round-trip closure residuals, observed branch occupancy, branch-to-branch orientation-change matrices, and can match an independently measured later-cycle B19′ EBSD dataset. A metric-aware NiTi natural/AQ quick setup constructs the OR from `(010)B19′ ∥ (110)B2` and `[101]B19′ ∥ [−1 1 1]B2` using direct and reciprocal lattice metrics rather than treating monoclinic Miller indices as Cartesian vectors.

**NiTi CT + Otsuka–Ren reconstruction bridge:** the reconstruction workspace now also exposes the source-derived B2/B19′ correspondence matrix and measured lattice metrics used by Correspondence Theory. Because correspondence is not itself a unique experimental orientation relationship, the app transparently forms a metric-aware correspondence deformation and uses its polar rotational factor only as a reproducible **model-derived starting OR**. The natural/AQ route, KS/NW/Bain/Pitsch/Burgers presets, custom parallelisms and custom OR matrix remain available unchanged. The parent/daughter home page is grouped into four academic workspaces: reconstruct parent, NiTi cycle, forward/batch, and academic guide.

## Main workspaces

### 1. Compatibility dashboard

A compact view of:

- \(\lambda_1,\lambda_2,\lambda_3\)
- \(|\lambda_2-1|\)
- CMC degeneracy distance
- best normalized CC2 residual
- CC3 margin
- metric shear/shear mismatch \(\varepsilon\), when an exact CMC habit plane exists
- explicit warning when the classical cofactor test and the metric shear/shear construction do not give the same certification
- interactive 3D CMC zero surface

### 2. PTMC / cofactor cross-check

Implements the full nonlinear conditions

\[
\mathrm{CC1}:\;\lambda_2=1
\]

\[
\mathrm{CC2}:\;a\cdot U\,\operatorname{cof}(U^2-I)n=0
\]

\[
\mathrm{CC3}:\;\operatorname{tr}(U^2)-\det(U^2)-\frac{|a|^2|n|^2}{4}-2\ge0.
\]

Non-trivial cubic two-fold generators are grouped by the symmetry-related stretch they produce. Singleton generators are evaluated as Type-I/Type-II domains; paired perpendicular generators are evaluated with the dedicated Compound-domain rank-one solutions. The exact CC2 value, a normalized CC2 residual, applicable simplified Type-I/II checks, and the CC3 margin are shown.

### 3. Complete twin explorer

For the built-in B2→B19′ correspondence:

- full cubic parent point group: **48 operations**
- correspondence intersection subgroup: **4 operations**
- double-coset partition: **7 classes**
- symmetry-equivalent two-fold generators are retained and labeled
- Type-I and Type-II transformation-twin elements are calculated from correspondence, metrics and parent symmetry
- compound twins are detected when the Type-I and Type-II descriptions coincide for the same intercorrespondence class
- when CMC gives habit planes, every candidate is ranked by the metric shear/shear mismatch \(\varepsilon\)

### 4. Temperature & uncertainty

**Temperature sweep** accepts measured/fitted lattice parameters as a function of temperature and recomputes CMC/cofactor residuals row by row.

**Uncertainty propagation** performs Monte Carlo propagation of user-supplied 1σ lattice-parameter uncertainties. Reported fractions mean **fraction of samples inside explicitly chosen tolerances**. They are not presented as the probability of satisfying an exact equality.

### 5. Literature, microstructure & functional performance

A curated database separates:

- experimental benchmarks
- mechanism studies
- theory
- data-driven design methods
- device demonstrations

Microstructure inputs—grain size, distribution, precipitates, retained martensite, dislocation density and processing—are stored alongside the exact geometry calculations. The app intentionally **does not fabricate a universal fatigue-life equation** from these descriptors.

Optional functional measurements include transformation temperatures, thermal/stress hysteresis, \(\Delta T_{ad}\), COP, recoverable strain and cycle count.

### 6. Inverse lattice design

Differential-evolution and Latin-hypercube/Pareto search operate directly in normalized lattice space \((a,b,c,\beta)\). The objective is composed of explicit residuals:

- \(|\lambda_2-1|\)
- normalized CC2 residual
- CC3 violation only
- CMC degeneracy distance
- optional proximity penalty to the starting alloy

The objective weights are visible and user-controlled.

### 7. Composition → lattice ML

No universal pretrained model is shipped. Instead:

1. upload a real composition/processing dataset;
2. select numerical descriptors;
3. cross-validate linear regression, random forest or extra-trees regression;
4. inspect MAE/RMSE/\(R^2\) for every lattice target;
5. predict candidate lattice parameters;
6. pass the predicted lattices through the exact compatibility engine.

This prevents a low predicted compatibility residual from being treated as meaningful when the lattice model itself is poorly validated.

### 8. Multi-step transformations

A generic metric/stretch engine accepts phase chains such as B2→R→B19′. Each stage has its own general triclinic unit cell and 3×3 correspondence matrix. The generic engine reports stagewise stretch eigenvalues and CMC distance.

Full built-in twin classification is **not silently generalized** to arbitrary phase pairs; that requires a phase-specific symmetry/correspondence module.

### 9. Parent ↔ daughter orientation reconstruction

The reconstruction workbench accepts either pre-segmented daughter grains or raw `.ang` / `.ctf` / delimited EBSD point maps. Raw imports are audited before use: phase 0 is treated as non-indexed, the daughter phase is selected explicitly, vendor quality filters are optional and off by default, and the user must confirm the Euler/map reference-frame convention before segmentation. Small raw maps can be segmented in-app with explicit disorientation, minimum-point and spatial-neighbor controls; production-scale maps can instead supply externally segmented grains and measured adjacency.

Five independently testable reconstruction routes can be run on the same data:

- neighbor voting;
- grain graph + Markov clustering;
- candidate-level variant graph propagation;
- nucleation + growth;
- operator / groupoid consistency with per-edge operator IDs and residuals.

A bounded OR-refinement routine can reduce neighbor inconsistency near a supplied physical OR. Built-in KS, NW, Bain, Pitsch and Burgers OR families provide symmetry-generated variants and forward parent→daughter prediction. The output now exposes best/second-best candidate residuals, candidate separation, heuristic support decomposition, parent-level variant diversity, operator frequencies and full parent orientations. Cross-method evidence includes ARI, NMI, prior-parent-boundary Jaccard agreement, overlap-matched parent-orientation disagreement and a boundary-consensus table. A one-click academic evidence ZIP exports the analyzed input, exact OR, all controls, per-method parent/daughter tables, variant/operator statistics and cross-method matrices. Approximate centroid k-NN adjacency remains available only as a visibly labelled fallback when true shared-boundary adjacency is unavailable. See `docs/RECONSTRUCTION.md`.

### 10. Independent compatibility methods

The app now cross-checks the main compatibility engine against additional exact or theorem-based diagnostics:

- generic Hadamard rank-one jump condition;
- single-variant Ball–James middle-stretch condition;
- pairwise martensite-variant rank-one spectral criterion;
- direct all-volume-fraction laminate singular-value scan;
- cubic→orthorhombic triplet-condition residuals with an explicit applicability guard.

These methods answer different questions and are never collapsed into a fabricated universal percentage score.

### 11. Research frontier & methods

The app reports stretch-variant and commutator diagnostics and includes a guarded implementation of current extreme-compatibility reference tensors. Those tensors are applied only when the transformation belongs to their stated crystallographic applicability class. For the built-in B2→B19′ correspondence, the monoclinic unique axis maps to a cubic `<110>` direction, so cubic→monoclinic-II extreme targets are **not** used as a pass/fail criterion.

### 12. Equation explorer & analytical solutions

A dedicated audit workspace exposes the actual mathematics behind the computed numbers rather than treating equations as hidden implementation details. It includes:

- explicit PTMC IPS construction `v → R → F=RU=I+d m^T` with rank-one and plane-invariance residuals;
- analytical B2/B19′ CMC and closed-form eigenspectrum cross-checked against the general metric engine;
- C1/C2a/C2b/C3 and D1/D2/E degeneracy families;
- source-defined distance calculations with discrepancy auditing;
- left-coset correspondence variants, CMC zero-set symmetries and double-coset intercorrespondence classes;
- closed-form O2 and Appendix-C O4 supercompatibility families;
- searchable equation/method catalog with source, implementation and scope.

Known source inconsistencies are surfaced instead of silently corrected. See `docs/SOURCE_DISCREPANCIES.md`.

### 13. Paper-ready audit & export

Creates a research record containing physical inputs and units, normalized ratios, correspondence matrices, numerical tolerances, software build, core matrices, principal stretches, CMC eigenvalues/habit planes, analytical parity residuals, equation provenance, source discrepancies and claim boundaries. Exports are available as JSON, Markdown and CSV and are fingerprinted with SHA-256.

This is designed to make a Methods/Results calculation traceable; it does **not** turn crystallographic calculations into unsupported claims about fatigue, hysteresis, chemistry or experimental reversibility.

---

## Friendly physical input

The primary UI requests real unit-cell quantities with full names and units:

| Input | Unit | Meaning |
|---|---:|---|
| B2 austenite lattice parameter — `a_B2` | Å | cubic parent cell edge |
| B19′ martensite lattice parameter — `a_B19′` | Å | martensite a-axis |
| B19′ martensite lattice parameter — `b_B19′` | Å | martensite b-axis |
| B19′ martensite lattice parameter — `c_B19′` | Å | martensite c-axis |
| B19′ monoclinic angle — `β` | ° | angle between martensite a- and c-axes |

Normalized \(a,b,c\) ratios are calculated internally and shown as outputs, not ambiguous input symbols.

---

## Reproducibility and validation

The current test suite checks:

- published B2/B19′ normalized lattice ratios;
- CMC matrices and degeneracy;
- first-order habit-plane solutions;
- SMC shear vectors;
- transformation-twin shear amplitudes;
- metric shear/shear \(\varepsilon\) and angles;
- construction of \(U\) and \(\lambda_2\);
- classical Type-I/Type-II/Compound cofactor conditions;
- the deliberate cross-framework discrepancy in a C1 teaching benchmark;
- cubic symmetry-group size, correspondence subgroup and double-coset partition;
- temperature evaluation;
- zero-uncertainty Monte Carlo behavior;
- generic multi-step stretch equivalence;
- frontier applicability guard;
- user-trained ML and physics screening;
- integrity of the curated literature database;
- ideal KS/NW/Bain/Pitsch/Burgers variant counts;
- five-method synthetic parent-reconstruction accuracy and OR refinement;
- Hadamard rank-one diagnostics and pairwise stretch-variant compatibility;
- direct all-volume-fraction cofactor verification;
- exact constructed triplet-condition benchmarks.

Run:

```bash
pip install -r requirements-dev.txt
pytest -q
```

GitHub Actions runs the same checks on pushes and pull requests.

See:

- [`docs/EQUATION_PROVENANCE.md`](docs/EQUATION_PROVENANCE.md) — equation keys, source/extension classes and implementation traceability
- [`docs/NOTATION_GUIDE.md`](docs/NOTATION_GUIDE.md) — explicit meaning and units for every formula/input/output symbol
- [`docs/SOURCE_DISCREPANCIES.md`](docs/SOURCE_DISCREPANCIES.md) — explicit register of source inconsistencies and implementation policy
- [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md) — paper-ready record contract and minimum reporting fields
- [`docs/PAPER_METHODS_TEMPLATE.md`](docs/PAPER_METHODS_TEMPLATE.md) — manuscript-methods checklist
- [`docs/METHODS.md`](docs/METHODS.md) — equations, conventions and algorithmic details
- [`docs/VALIDATION.md`](docs/VALIDATION.md) — benchmark philosophy and expected numerical checkpoints
- [`docs/REFERENCES.md`](docs/REFERENCES.md) — primary literature and what each source supports
- [`docs/PAPER_WORKFLOW.md`](docs/PAPER_WORKFLOW.md) — a defensible workflow for a research study
- [`docs/DATA_SCHEMA.md`](docs/DATA_SCHEMA.md) — CSV schemas and units
- [`docs/RECONSTRUCTION.md`](docs/RECONSTRUCTION.md) — orientation conventions, OR families and parent-reconstruction algorithms

---

## Install and run

Recommended: Python 3.12.

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
pip install -r requirements.txt
streamlit run supercompatibility_final.py
```

### Streamlit Community Cloud

1. Upload this repository to GitHub.
2. Create a Streamlit Community Cloud app from the repository.
3. Set the entry point to `supercompatibility_final.py`.
4. Deploy. No secrets are required.

---

## Repository structure

```text
supercompatibility-lab/
├── supercompatibility_final.py # recommended fresh self-contained Streamlit entrypoint
├── supercompatibility_r7.py  # byte-identical compatibility fallback
├── app.py                    # byte-identical self-contained fallback
├── streamlit_app.py          # byte-identical self-contained fallback
├── src/
│   ├── core.py              # metric correspondence, CMC, SMC, metric twins
│   ├── ptmc.py              # stretch tensor and CC1–CC3
│   ├── symmetry.py          # point groups, double cosets, full twin explorer
│   ├── distances.py         # compatibility residual dashboard
│   ├── temperature.py       # temperature-series evaluation
│   ├── uncertainty.py       # Monte Carlo error propagation
│   ├── literature.py        # curated literature loader
│   ├── microstructure.py    # evidence-linked context, no fake fatigue model
│   ├── performance.py       # optional measured functional metrics
│   ├── design.py            # inverse and Pareto lattice design
│   ├── ml.py                # user-trained composition→lattice ML
│   ├── multistep.py         # general stage-wise metric/stretch engine
│   ├── frontier.py          # guarded extreme-compatibility diagnostics
│   ├── reconstruction.py    # parent/daughter ORs and five reconstruction methods
│   ├── compatibility_methods.py # independent rank-one, all-f and triplet checks
│   ├── presets.py
│   └── visualization.py
├── data/
├── docs/
├── scripts/
├── tests/
├── .github/workflows/tests.yml
├── requirements.txt
└── requirements-dev.txt
```

---

## Scientific scope and non-claims

This software is intended for **research screening, reproducible analysis and hypothesis generation**. It does not claim that crystallographic supercompatibility alone determines fatigue life, hysteresis, elastocaloric performance or synthesizability. Those outcomes also depend on microstructure, defects, processing, kinetics, loading and other material-specific factors.

Likewise, equality conditions such as \(\lambda_2=1\), CC2=0 and \(\varepsilon=0\) are exact mathematical statements. Experimental and floating-point calculations require tolerances; the software exposes those tolerances rather than hiding them.

## License

MIT — see [`LICENSE`](LICENSE).
