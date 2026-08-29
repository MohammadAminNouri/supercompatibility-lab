# Research workflow for a publication-grade study

This is a suggested workflow for using Supercompatibility Lab in a paper without overclaiming what the software establishes.

## 1. Define the scientific question before optimizing

Examples:

- Does a composition series approach single-variant A/M compatibility as temperature or composition changes?
- Which conventional cofactor domain system is closest to CC1–CC3?
- Does metric shear/shear intercompatibility track measured hysteresis better than \(|\lambda_2-1|\) alone?
- Can inverse lattice design identify experimentally reachable target geometries?
- Does grain size or precipitate state explain residual performance variation after crystallographic metrics are accounted for?

Do not begin by optimizing a composite score and only later assign physical meaning to it.

## 2. Acquire lattice data with uncertainty

For each alloy/condition record:

- composition and processing history
- temperature of diffraction measurement
- \(a_{B2}\)
- \(a_{B19'}\), \(b_{B19'}\), \(c_{B19'}\), \(\beta\)
- estimated standard uncertainties, ideally covariance from refinement
- phase fraction / confidence that the assigned phases and lattice model are appropriate

The UI uses Å and degrees; preserve original significant figures in the raw dataset.

## 3. Freeze crystallographic conventions

Before batch analysis, report:

- parent/product phase labels
- unit-cell conventions
- correspondence matrix
- whether a different correspondence model was tested
- reciprocal/direct-space convention for planes and directions

For the built-in specialized model, the correspondence matrix is printed in `docs/METHODS.md` and `src/core.py`.

## 4. Report three compatibility layers separately

### Layer A — A/M compatibility

Report at least

\[
\lambda_1,\lambda_2,\lambda_3,\quad |\lambda_2-1|,
\]

and the CMC eigenvalues / CMC degeneracy distance.

When an exact/within-tolerance CMC degeneracy occurs, report the habit-plane solution(s).

### Layer B — conventional cofactor supercompatibility

For every non-trivial tested two-fold domain axis and Type I/II system, preserve:

- axis \(\hat e\)
- Type I, Type II or Compound classification
- CC1 residual
- raw and normalized CC2 residual
- simplified CC2 cross-check only when applicable to non-compound Type-I/II systems
- CC3 margin
- explicit tolerances

Do not report only the best row unless the selection rule is predeclared.

### Layer C — metric shear/shear intercompatibility

When CMC supplies a compatible habit plane, report:

- habit plane \(m_A\)
- SMC shear \(d_A\)
- twin plane/normal and twin shear vector
- twin shear amplitude \(s\)
- \(\varepsilon\)
- shear-direction angular mismatch

If classical CC1–CC3 and metric \(\varepsilon\) disagree, report the discrepancy rather than choosing one silently.

If a **Compound** domain satisfies CC1–CC3, state explicitly that this is not automatically a proof of transition-layer elimination for austenite/compound laminates. Treat extreme compatibility as a separate theorem with its own transformation-class assumptions.

## 5. Temperature studies

Use `data/temperature_sweep_template.csv`.

Recommended figures:

- \(|\lambda_2-1|\) vs T
- normalized CC2 vs T
- CC3 margin vs T
- CMC distance vs T
- measured hysteresis / transformation temperature on a separate axis or aligned panel

Avoid taking a minimum at one temperature as representative of the full operating window.

## 6. Propagate measurement uncertainty

At minimum perform a sensitivity study around the refined lattice parameters.

If only independent standard errors are available, the built-in Gaussian Monte Carlo is a transparent first model. Report:

- assumed σ values
- number of samples
- random seed
- CC1/CC2/CMC tolerances
- fraction within each tolerance

If the diffraction refinement covariance matrix is available, extend `src/uncertainty.py` to sample the multivariate covariance rather than treating parameters as independent.

## 7. Correlate with performance without assuming causation

Store measured:

- thermal hysteresis
- stress hysteresis
- \(\Delta T_{ad}\)
- transformation temperatures
- cycle life
- recoverable strain
- COP, when applicable

Then test associations with compatibility residuals statistically. Possible analyses include Spearman correlation, robust regression or mixed effects if repeated measurements exist. These statistical layers are intentionally not hard-coded as “laws” in the current app.

## 8. Include microstructure as an independent descriptor block

Recommended descriptors:

- mean grain size and distribution shape
- precipitate type, size and fraction
- retained martensite
- dislocation density
- texture / orientation distribution if available
- processing history

A useful paper question is whether geometry explains performance after controlling for microstructure, rather than assuming one replaces the other.

## 9. Inverse design

When using the lattice optimizer:

1. define bounds from experimentally plausible lattice variations;
2. predeclare the target domain axis/type;
3. report all objective weights;
4. repeat with multiple seeds;
5. inspect the Pareto set rather than a single scalar optimum;
6. verify the final candidate directly with all exact residual calculations;
7. compare target lattice shifts with known alloying/temperature trends before proposing synthesis.

An optimizer minimum is a **crystallographic target**, not evidence that a chemistry producing that target exists.

## 10. ML study

For a composition→lattice paper:

- keep the exact compatibility engine fixed;
- split data by chemistry family or processing group when leakage is possible;
- report per-target MAE/RMSE/R²;
- compare against a simple linear baseline;
- avoid extrapolating far outside the training descriptor range;
- propagate lattice-prediction uncertainty into compatibility if possible;
- rank candidates only after physics screening.

A good research contribution would be **uncertainty-aware inverse composition design whose objective is a vector of exact crystallographic residuals**, not a black-box binary label.

## 11. Multi-step transformations

For B2→R→B19′ or other chains, analyze each stage separately first. The generic module evaluates metric/stretch compatibility but does not infer stage-specific transformation twins without a supplied symmetry/correspondence model.

A paper should state clearly which stages have full cofactor/twin certification and which have metric-only screening.

## 12. Frontier theory

Before using a frontier compatibility theorem:

- verify transformation-class assumptions;
- verify variant construction and symmetry class;
- verify the theorem's exact tensor parameterization;
- only then compare to target tensors/conditions.

The built-in guard demonstrates this philosophy by refusing to certify the built-in B2→B19′ correspondence with a cubic→monoclinic-II exact target that does not match its axis mapping.

## 13. Reproducibility package for a paper

Archive with the manuscript:

- exact Git commit/tag
- raw lattice/performance CSVs
- processed result CSVs/JSONs downloaded from the app
- optimizer bounds/weights/seeds
- ML feature list/model/hyperparameters/folds/seed
- uncertainty model and seed
- software environment (`requirements.txt` or lock file)
- `pytest -q` output

The repository's `scripts/reproduce_benchmarks.py` can be included in supplementary material as an executable numerical audit.

## Parent-reconstruction study workflow

For a reconstruction-focused paper or methods section:

1. State the parent and daughter point groups and the initial physical OR.
2. Report the EBSD grain-reconstruction/segmentation settings used before this app; those preprocessing choices are outside this repository.
3. Supply the measured grain adjacency graph rather than centroid k-NN whenever possible.
4. Run at least two reconstruction families (for example, variant graph and operator/groupoid consistency) and report disagreements at boundaries instead of selecting the visually preferred answer.
5. If OR refinement is used, report the starting OR, correction bound, objective before/after and correction angle.
6. Report angular fit distributions, confidence, reconstructed parent count and variant frequencies, not only a colored map.
7. Validate against retained-parent regions, independently reconstructed regions, synthetic benchmarks or another established implementation when available.
8. Export the method-comparison CSV and retain it with the manuscript analysis archive.

## Cross-compatibility study workflow

A defensible compatibility paper can report a hierarchy rather than one scalar score:

1. CMC degeneracy and habit-plane result;
2. \(\lambda_1,\lambda_2,\lambda_3\) and \(|\lambda_2-1|\);
3. complete pairwise stretch-variant rank-one table;
4. classified CC1–CC3 domain results;
5. direct all-volume-fraction singular-value scan for the selected domain system;
6. SMC/shear-shear residual when a CMC habit plane exists;
7. triplet-condition residuals only for a transformation class for which the implemented TC equations are valid;
8. uncertainty and temperature dependence for any experimentally claimed near-zero residual.

If two frameworks disagree, report the disagreement. Do not tune numerical tolerances post hoc to force agreement.
