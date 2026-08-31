# Supercompatibility Lab

**Supercompatibility Lab** is a research-oriented computational workbench for crystallographic compatibility, martensitic transformations, correspondence variants, transformation twinning, interaction work, EBSD parent–daughter reconstruction, inverse lattice design, uncertainty analysis, independent theoretical cross-checks, and reproducible manuscript reporting.

The principal implementation focuses on the **B2 ↔ B19′ martensitic transformation in NiTi**, while several modules support more general transformation-matrix and compatibility problems.

The aim is to provide a single traceable research workflow in which lattice measurements, crystallographic theory, mechanical selection, reconstruction, numerical validation, and publication evidence remain connected without treating distinct theories as interchangeable.

---

## Scientific workflow

For a typical NiTi study:

```text
Experimental lattice parameters
        │
        ▼
Compatibility diagnosis
        │
        ├──────────► PTMC / cofactor verification
        │
        ▼
Variants, operators and twins
        │
        ├──────────► Interaction Work
        │
        ├──────────► EBSD parent reconstruction
        │
        ▼
Inverse lattice design
        │
        ▼
Uncertainty / independent checks
        │
        ▼
Manuscript audit and reproducible export
```

At the lattice level:


$$
\text{lattice parameters}
\rightarrow
\text{metrics}
\rightarrow
\text{transformation matrices}
\rightarrow
\text{stretch/distortion tensors}
$$


At the compatibility level:


$$
\text{A/M compatibility}
\rightarrow
\text{M/M twinning}
\rightarrow
\text{A/M--M/M shear intercompatibility}
$$


At the variant level:


$$
\text{symmetry}
\rightarrow
\text{correspondence variants}
\rightarrow
\text{intervariant operators}
\rightarrow
\text{twin relations}
$$


At the mechanical level:


$$
\text{orientation}
+
\text{stress}
+
\text{deformation}
\rightarrow
\text{Interaction Work}
$$


At the experimental level:


$$
B19'
\text{ EBSD}
\rightarrow
B2
\text{ parent reconstruction}
\rightarrow
\text{variant identification}
$$


At the design level:


$$
\text{compatibility objective}
\rightarrow
\text{inverse lattice search}
\rightarrow
\text{candidate target metrics}
$$


These branches are deliberately kept distinct. Agreement between them can strengthen an interpretation; one method is not silently substituted for another.

---

# Research workspaces

## 1 · Compatibility verdict & diagnostics

The first workspace provides the quickest quantitative assessment of the active lattice.

Typical outputs include:

- normalized lattice ratios;
- transformation metrics;
- principal stretch eigenvalues;
- $\lambda_2$;
- $|\lambda_2-1|$;
- CMC eigenvalues;
- normalized distance from exact CMC degeneracy;
- existence or absence of an exact A/M habit-plane solution;
- availability of downstream SMC and shear/shear calculations;
- explicit numerical PASS / FAIL / BLOCKED criteria.

A result is not reported only as a qualitative label. Paper-facing output is designed to retain:


$$
\boxed{\text{value}}
\qquad
\boxed{\text{criterion}}
\qquad
\boxed{\text{residual}}
\qquad
\boxed{\text{interpretation}}
$$


This workspace is a diagnostic summary. Detailed theory remains available in the specialized workspaces.

---

## 2 · Classical PTMC / cofactor verification

This workspace independently evaluates the classical stretch-tensor and cofactor-condition framework.

For the transformation stretch tensor


$$
U,
$$


the ordered principal stretches are


$$
\lambda_1 \leq \lambda_2 \leq \lambda_3.
$$


A central compatibility quantity is


$$
|\lambda_2-1|.
$$


Applicable Type-I, Type-II and compound twin domains are then evaluated through the individual cofactor conditions.

The software keeps:

- CC1;
- CC2;
- CC3;

separate.

A domain is not declared to satisfy the tested cofactor framework unless every required condition satisfies its declared numerical criterion.

Raw matrices and complete domain tables are available for audit but are not required in the main workflow.

---

## 3 · Variants, twins & Interaction Work

This workspace links crystallographic variant theory, transformation twins, A/M–M/M compatibility, and mechanical variant selection.

### Correspondence variants

The B2 parent symmetry is reduced by the correspondence subgroup.

For the implemented NiTi correspondence:


$$
N_C^M
=
\frac{|G^A|}{|H_C^A|}
=
12.
$$


The software therefore produces twelve B19′ correspondence variants.

The user can inspect an explicit pair:


$$
V_i \leftrightarrow V_j.
$$


### Intercorrespondence operators

Variant relationships are grouped into symmetry-equivalent operator classes through double-coset decomposition.

For a selected pair, the application reports:

- operator class;
- mirror content;
- 180° rotational content;
- Type-I twin capability;
- Type-II twin capability;
- complete $V_i\rightarrow V_j$ operator map when requested;
- reduced multivalued operator-composition table.

The full group-theoretical matrices are retained for audit but remain secondary to the physical interpretation.

### Transformation twins

For an applicable variant relationship the application can calculate M/M twin solutions including:

- Type I;
- Type II;
- compound-equivalent cases.

For each selected solution the program distinguishes explicitly between:

- parent symmetry generator;
- calculated twin shear direction;
- calculated twin-plane normal;
- shear magnitude;
- B2-coordinate representation;
- B19′ crystallographic representation where justified;
- class-representative geometry;
- pair-specific information.

A parent mirror or twofold generator is not automatically presented as if it were the conventional martensite twin-plane indexing.

### Supercompatibility bridge

The logical dependency is:


$$
\boxed{\text{A/M compatibility}}
$$


$$
\downarrow
$$


$$
\boxed{\text{A/M habit plane and shear}}
$$


$$
\downarrow
$$


$$
\boxed{\text{M/M twin compatibility}}
$$


$$
\downarrow
$$


$$
\boxed{\text{A/M--M/M shear matching}}
$$


$$
\downarrow
$$


$$
\boxed{\text{supercompatibility verdict}}
$$


If the exact A/M prerequisite does not exist, the final shear/shear quantity is reported as **NOT COMPUTED**, together with the reason. No artificial value is generated.

### Interaction Work

Interaction Work evaluates the mechanical work provided by an external stress to a declared martensitic deformation:


$$
IW
=
\boldsymbol{\sigma}:\boldsymbol{\varepsilon}.
$$


The deformation strain is defined as


$$
\boldsymbol{\varepsilon}
=
F-I.
$$


For martensite reorientation from an initial state $i$ to a destination state $j$,


$$
F_{i\rightarrow j}
=
F_jF_i^{-1},
$$


hence


$$
\varepsilon_{i\rightarrow j}
=
F_jF_i^{-1}-I.
$$


The software evaluates


$$
IW_{i\rightarrow j}
$$


for all candidate destination states and reports the maximum-work candidate.

Inputs can include:

- tension or compression;
- stress magnitude;
- parent crystallographic loading direction;
- specimen loading direction;
- B2 orientation;
- Bunge Euler angles;
- initial martensite state.

When an EBSD-derived parent orientation is supplied, the stress/loading representation is transformed consistently between specimen and crystal coordinates.

Interaction Work is a **mechanical selection criterion**.

It is not treated as:

- proof of transformation occurrence;
- a complete kinetic law;
- proof of supercompatibility;
- a replacement for an activation-energy model.

---

## 4 · Temperature & measurement uncertainty

This workspace addresses temperature dependence and uncertainty in experimental lattice measurements.

### Temperature dependence

Temperature-dependent lattice data can be supplied to evaluate quantities such as:


$$
|\lambda_2(T)-1|
$$


and


$$
\delta_{\mathrm{CMC}}(T).
$$


Outputs can include:

- supplied temperature range;
- temperature giving minimum $|\lambda_2-1|$;
- temperature giving minimum CMC distance;
- rows satisfying declared compatibility criteria;
- compatibility trends across the experimental range.

The application does not invent temperature dependence if it is absent from the supplied dataset.

### Measurement uncertainty

Monte-Carlo propagation can be performed using declared uncertainties in:


$$
a_{B2},
\quad
a_{B19'},
\quad
b_{B19'},
\quad
c_{B19'},
\quad
\beta.
$$


Typical outputs include:

- valid sample count;
- $|\lambda_2-1|$ distribution;
- CMC-distance distribution;
- quantiles;
- fraction satisfying a declared tolerance.

The reported fraction describes the supplied measurement/uncertainty model. It is not a universal probability that the material is compatible.

---

## 5 · Experimental context & literature

Crystallographic compatibility is only one part of martensitic materials behavior.

This workspace stores experimental information such as:

- composition;
- processing route;
- heat treatment;
- grain size;
- precipitates;
- texture;
- dislocation state;
- transformation temperatures;
- thermal hysteresis;
- mechanical hysteresis;
- fatigue cycles;
- functional degradation;
- experimental notes.

This information is retained as experimental evidence.

The application does not invent deterministic fatigue, hysteresis, or performance laws from crystallographic compatibility alone.

Literature records can be compared with the active study without mixing reported experimental observations with software predictions.

---

## 6 · Inverse lattice design

The inverse-design workspace asks:

> Which lattice parameters would move the current system toward a selected compatibility condition?

The optimizer may vary:


$$
a_{B19'},
\quad
b_{B19'},
\quad
c_{B19'},
\quad
\beta
$$


within explicit user-defined bounds.

Objectives can emphasize:

- middle-eigenvalue compatibility;
- CMC;
- cofactor conditions;
- combinations of compatibility residuals.

Outputs include:

- original lattice;
- optimized target lattice;
- absolute parameter changes;
- percentage changes;
- normalized lattice ratios;
- $|\lambda_2-1|$;
- cofactor diagnostics;
- CMC diagnostics;
- Pareto trade-offs when requested.

An optimized result is a **mathematical lattice target**.

It is not automatically a realizable alloy composition.

Experimental/compositional feasibility must be assessed independently.

---

## 7 · Composition → lattice modelling

This workspace connects empirical composition/processing datasets to the crystallographic physics engine.

A training dataset may contain:

- elemental concentrations;
- processing variables;
- heat-treatment descriptors;
- measured B2 lattice parameters;
- measured B19′ lattice parameters;
- measured monoclinic angle.

The workflow is:

```text
composition / processing data
        ↓
cross-validated predictive model
        ↓
predicted lattice parameters
        ↓
compatibility calculations
        ↓
candidate screening
```

Model quality is evaluated **before** compatibility screening.

The application reports metrics such as cross-validation error and $R^2$. Weak or negative predictive performance is flagged explicitly.

Compatibility calculations cannot convert an unreliable regression model into a reliable composition prediction.

No universal NiTi composition→lattice law is assumed.

---

## 8 · Multi-step transformations

This workspace supports transformation sequences containing intermediate phases, for example:


$$
A\rightarrow B\rightarrow C.
$$


For each transformation step the user supplies the required phase/cell and correspondence information.

Per-stage outputs can include:

- transformation matrices;
- stretch tensors;
- principal stretches;
- $|\lambda_2-1|$;
- CMC diagnostics where mathematically applicable.

Generic transformation stages are not automatically assigned NiTi-specific twin systems or monoclinic B19′ formulas.

Applicability is kept explicit.

---

## 9 · EBSD parent ↔ daughter reconstruction

This workspace reconstructs parent grains from measured daughter-phase EBSD data.

For the principal NiTi use case:


$$
B19'
\rightarrow
B2.
$$


The reconstruction workflow is conceptually:

```text
measured B19′ orientations
        ↓
daughter-grain segmentation
        ↓
daughter adjacency
        ↓
theoretical variant/operator relations
        ↓
parent clustering
        ↓
reconstructed B2 grains
        ↓
reconstructed B2 orientations
        ↓
daughter variant assignments
```

Depending on the selected reconstruction strategy, outputs may include:

- reconstructed parent ID;
- reconstructed parent orientation;
- daughter→parent assignment;
- candidate variant;
- OR residual;
- second-best residual;
- candidate separation;
- support score;
- review flag;
- variant statistics;
- operator statistics;
- parent summaries.

### Known-truth validation

When truth labels are available, reconstruction is assessed using metrics such as:

- Adjusted Rand Index;
- normalized mutual information;
- homogeneity;
- completeness;
- V-measure;
- reconstructed-parent count;
- singleton-parent fraction;
- fragments per true parent;
- boundary precision;
- boundary recall;
- boundary F1;
- boundary Jaccard;
- false-parent-boundary rate.

A small OR residual is not treated as proof of successful parent reconstruction.

Likewise:


$$
\boxed{\text{method agreement}}
\neq
\boxed{\text{accuracy}}.
$$


Cross-method agreement and truth-referenced validation remain separate.

---

## 10 · Independent theorem cross-checks

This workspace evaluates selected conclusions through additional mathematical frameworks.

Every method is labeled by applicability.

Possible classifications include:

- applicable to the current NiTi transformation;
- generic matrix-level condition;
- specialized for another transformation class;
- data-dependent.

Available checks can include:

- variant-pair spectral compatibility;
- laminate/all-volume-fraction tests;
- specialized triplet conditions where applicable;
- generic Hadamard rank-one jump diagnostics.

A theorem derived for another crystal-transformation class is not silently presented as a B2→B19′ NiTi result.

---

## 11 · Applicability, references & methods

This workspace acts as the methodological map of the application.

It answers:

- Which method applies?
- What inputs are required?
- What physical/mathematical question does it answer?
- What can be concluded?
- What cannot be concluded?
- Is the method lattice-level, variant-level, mechanical, reconstruction-level, or empirical?
- Which literature source/equation supports it?

This workspace is intended to assist the preparation of manuscript **Methods** sections and to reduce accidental overclaiming.

---

## 12 · Equation & analytical solution library

The equation library exposes the mathematics used throughout the application.

It includes, where available:

- metric tensors;
- normalized lattice relations;
- transformation/stretch relations;
- CMC equations;
- IPS relations;
- habit-plane formulas;
- cofactor conditions;
- twin equations;
- correspondence groups;
- left cosets;
- double cosets;
- operator relations;
- analytical compatibility families;
- supercompatibility equations;
- inverse-design equations.

Symbols are defined beside the equations whenever practical.

Where an independent analytical result exists, the software can compare it with the general numerical implementation.

These parity checks form part of the software validation strategy.

---

## 13 · Manuscript audit & reproducible export

The final workspace is the publication/reproducibility checkpoint.

It records what was actually calculated.

The audit may contain:

- lattice inputs;
- units;
- normalized ratios;
- temperature;
- uncertainty assumptions;
- numerical tolerances;
- selected correspondence;
- selected variant/operator/twin;
- loading/orientation inputs;
- equations used;
- numerical residuals;
- PASS / FAIL / BLOCKED criteria;
- applicability restrictions;
- source discrepancies;
- reconstruction settings;
- optimization settings;
- software build;
- record hash.

A complete evidence bundle can be exported for manuscript preparation, supplementary information, or internal archiving.

Depending on the analyses performed during the session, the archive may contain:

```text
study metadata
compatibility diagnostics
analytical checks
variant/operator results
twin results
Interaction Work results
temperature results
Monte-Carlo uncertainty results
inverse-design target
Pareto scan
ML model/screening results
multi-step calculations
reconstruction results
microstructure record
functional-performance record
equation provenance
references
software/build metadata
SHA-256 fingerprints
```

The export is intended to make later reproduction of the analysis possible.

---

# NiTi input convention

The principal NiTi lattice input consists of:


$$
a_{B2}
$$


for cubic B2 austenite, and:


$$
a_{B19'},
\quad
b_{B19'},
\quad
c_{B19'},
\quad
\beta
$$


for monoclinic B19′ martensite.

Lengths are entered in ångström:


$$
1\text{ Å}=10^{-10}\text{ m}.
$$


The monoclinic angle $\beta$ is entered in degrees and converted internally where necessary.

Normalized quantities used by analytical theory are derived from the physical lattice values.

They are not treated as independent experimental inputs.

---

# Compatibility concepts

Several related but distinct compatibility concepts occur in martensitic-transformation theory.

The application does not treat them as synonyms.

## Middle-eigenvalue condition

A commonly used compatibility condition is:


$$
\lambda_2=1.
$$


The software therefore reports:


$$
|\lambda_2-1|.
$$


A small middle-eigenvalue residual is important, but it is not automatically sufficient for every stronger cofactor or supercompatibility theorem.

## CMC

The compatibility matrix condition is evaluated through its matrix/eigenvalue formulation.

Exact compatibility requires the relevant zero-eigenvalue degeneracy.

The application reports:

- all CMC eigenvalues;
- closest eigenvalue to zero;
- normalized zero-set distance;
- numerical tolerance;
- resulting status.

A `FAIL` therefore retains a quantitative measure of how far the lattice lies from the exact condition.

## SMC

Once the required A/M compatibility solution exists, the secondary construction supplies the corresponding shear/interface geometry used in later intercompatibility analysis.

SMC is not evaluated as if its prerequisite existed when CMC fails.

## Cofactor conditions

CC1, CC2 and CC3 are evaluated independently.

Passing CC1 alone does not imply that the complete tested cofactor conditions have been satisfied.

## Shear/shear intercompatibility

The final transformation-supercompatibility check compares the A/M shear geometry with an M/M twin shear.

The calculation is only performed when its mathematical prerequisites have been satisfied.

---

# Numerical tolerances

Floating-point calculations cannot establish literal symbolic equality.

Every numerical classification therefore depends on a declared tolerance.

Paper-facing results should retain:

- raw quantity;
- residual;
- normalization;
- tolerance;
- status.

Changing a tolerance changes the numerical classification criterion.

Tolerance values should therefore be archived with the calculation.

Default software tolerances are numerical defaults, not universal experimental limits.

---

# Interpretation and claim boundaries

The application deliberately separates mathematical calculations from broader physical conclusions.

## Compatibility calculations do not automatically predict

- low hysteresis;
- fatigue resistance;
- reversibility;
- transformation temperature;
- microstructural stability;
- functional life;
- experimentally realizable composition.

## A twin solution does not automatically prove

- that the twin will form experimentally;
- that it is kinetically preferred;
- that A/M compatibility exists;
- that full supercompatibility exists.

## A high Interaction Work does not automatically prove

- transformation occurrence;
- fastest kinetics;
- lowest total activation barrier;
- full crystallographic compatibility.

## A low reconstruction OR residual does not automatically prove

- correct parent clustering;
- correct parent count;
- correct boundaries;
- absence of fragmentation;
- agreement with known truth.

## A machine-learning prediction does not automatically prove

- physical realizability;
- reliable extrapolation;
- compatibility of a real alloy.

These claim boundaries should also be respected in publications using software-generated results.

---

# Installation

Clone the repository:

```bash
git clone <repository-url>
cd supercompatibility-lab
```

Create an environment:

```bash
python -m venv .venv
```

Activate on Linux/macOS:

```bash
source .venv/bin/activate
```

Activate on Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
```

Install dependencies:

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

For development/testing:

```bash
python -m pip install -r requirements-dev.txt
```

---

# Running the application

Primary entrypoint:

```bash
python -m streamlit run supercompatibility_final.py
```

The repository also maintains deterministic self-contained Streamlit entrypoints used by deployment/release verification.

These entrypoints are expected to remain synchronized according to the repository tests.

---

# Testing and release verification

Run the full regression suite:

```bash
python -m pytest -q
```

Additional verification scripts may include:

```bash
python scripts/deployment_preflight.py
python scripts/verify_embedded_app.py
python scripts/verify_equation_parity.py
python scripts/verify_release.py
python research_selftest.py
```

These checks address different failure modes.

## Deployment preflight

Checks deployment/package consistency.

## Embedded-app verification

Checks the embedded self-contained research engine and deterministic entrypoints.

## Equation parity

Checks selected analytical and numerical implementations against one another.

## Release verification

Checks repository/release structure and required assets.

## Scientific self-test

Runs numerical sanity checks on the scientific engine.

## Pytest

Runs the complete regression and policy suite.

A release is not considered validated merely because the Streamlit UI opens.

---

# Repository structure

The repository is organized around self-contained deployment entrypoints, scientific validation scripts, data assets, and tests.

A typical structure is:

```text
supercompatibility-lab/
│
├── app.py
├── streamlit_app.py
├── supercompatibility_r7.py
├── supercompatibility_final.py
│
├── requirements.txt
├── requirements-dev.txt
├── research_selftest.py
│
├── scripts/
│   ├── deployment_preflight.py
│   ├── verify_embedded_app.py
│   ├── verify_equation_parity.py
│   └── verify_release.py
│
├── data/
│   ├── templates
│   ├── demonstration datasets
│   └── research assets
│
├── tests/
│   └── regression and release tests
│
└── README.md
```

---

# Reproducibility

For publication-quality work, record at minimum:

- Git commit hash;
- software/build identifier;
- complete lattice input;
- units;
- temperature;
- uncertainty model where relevant;
- correspondence convention;
- symmetry convention;
- selected variant;
- selected operator;
- selected twin;
- loading direction;
- stress magnitude and sign;
- orientation convention;
- numerical tolerances;
- optimization bounds;
- ML dataset/version where relevant;
- reconstruction parameters;
- equation/source identifiers;
- exported evidence archive.

For long-term reproducibility, archive the evidence ZIP together with the exact Git commit used to create it.

---

# Literature basis

The implementation is based on published crystallographic and martensitic-transformation theory.

Important references include:

## Correspondence theory for NiTi

**The Correspondence Theory and Its Application to NiTi Shape Memory Alloys**  
*Crystals* 12 (2022) 130  
DOI: `10.3390/cryst12020130`

## NiTi martensite reorientation and EBSD

**An investigation on reorientation and textural evolution in a martensitic NiTi rolled sheet using EBSD**  
*International Journal of Plasticity* 159 (2022) 103468  
DOI: `10.1016/j.ijplas.2022.103468`

## Interaction Work

**The role of interaction work in martensite deformation**  
*Scripta Materialia* 256 (2025) 116433  
DOI: `10.1016/j.scriptamat.2024.116433`

## Correspondence/metric/symmetry supercompatibility

**Compatibilities and supercompatibility conditions in shape memory alloys determined from correspondence, metrics and symmetries**  
*Acta Materialia* 316 (2026) 122399  
DOI: `10.1016/j.actamat.2026.122399`

## Orientational variants and operator/groupoid methodology

**GenOVa: a computer program to generate orientational variants**  
*Journal of Applied Crystallography* 40 (2007) 1179–1182  
DOI: `10.1107/S0021889807048741`

## Parent-grain reconstruction methodology

**ARPGE: a computer program to automatically reconstruct the parent grains from electron backscatter diffraction data**  
*Journal of Applied Crystallography* 40 (2007) 1183–1188  
DOI: `10.1107/S0021889807048777`

The software is an independent implementation of published mathematical concepts. It does not reproduce source code or GUI assets from software described in the literature.

Where sources use different signs, bases, correspondence conventions, analytical forms, or interpretations, the adopted convention should be stated explicitly rather than silently combining incompatible expressions.

---

# Recommended academic workflow

For a paper using measured NiTi lattice parameters:

1. record the experimental lattice parameters and temperature;
2. run the compatibility verdict;
3. inspect the raw CMC and $\lambda_2$ residuals;
4. perform the PTMC/cofactor cross-check;
5. inspect relevant variant/operator/twin relationships;
6. use EBSD reconstruction when experimental parent orientations are required;
7. evaluate Interaction Work when mechanical variant selection is part of the research question;
8. propagate measurement uncertainty where relevant;
9. use inverse design only after documenting the original lattice state;
10. apply independent theorems only when their assumptions match the transformation;
11. export the complete manuscript evidence record;
12. archive the Git commit and evidence package.

Results should be reported using numerical quantities and explicit criteria rather than only software-generated PASS/FAIL labels.

---

# Suggested Methods statement

A concise manuscript description may be adapted from:

> Crystallographic calculations were performed using Supercompatibility Lab. Measured parent and martensite lattice parameters were supplied directly to the transformation-metric engine. Austenite/martensite compatibility, stretch-eigenvalue conditions, transformation-twin relations and applicable supercompatibility conditions were evaluated using explicitly defined correspondence, metric and symmetry conventions. Numerical classifications were based on recorded residual tolerances. Where independent analytical expressions were available, they were compared against the general numerical formulation. Input parameters, residuals, calculation settings, equation provenance and software version were retained in a reproducibility evidence record.

Only methods actually used in the study should be included in the final manuscript.

---

# Citation

Research using Supercompatibility Lab should cite:

1. the repository release, archived software version or commit used for the calculation; and
2. the original theoretical publications corresponding to the methods used.

The software should not be cited as a replacement for the underlying crystallographic theory.

---

# Software status

Supercompatibility Lab is research software intended for:

- crystallographic analysis;
- numerical verification;
- hypothesis testing;
- transformation-variant analysis;
- EBSD reconstruction;
- method comparison;
- uncertainty analysis;
- inverse design;
- experimental interpretation;
- reproducible manuscript preparation.

All results remain conditional on the assumptions of the selected mathematical framework and the quality of the supplied experimental data.

