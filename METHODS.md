# Methods and mathematical conventions

This document defines the equations implemented by **Supercompatibility Lab**. It is written to make the numerical workflow auditable and suitable for a Methods/Supplementary Methods section.

## 1. Built-in crystallographic system

The fully specialized compatibility/twin module uses the B2 parent → B19′ product transformation with monoclinic B19′ unique b axis.

The user supplies physical lattice parameters

\[
a_{B2},\;a_{B19'},\;b_{B19'},\;c_{B19'},\;\beta.
\]

The UI uses Å for lengths and degrees for \(\beta\). Internally, the B19′ lengths are normalized by \(a_{B2}\):

\[
a=\frac{a_{B19'}}{a_{B2}},\qquad
b=\frac{b_{B19'}}{a_{B2}},\qquad
c=\frac{c_{B19'}}{a_{B2}}.
\]

The normalized parent and product metrics are

\[
M_A=I,
\]

\[
M_M=
\begin{pmatrix}
a^2 & 0 & ac\cos\beta\\
0 & b^2 & 0\\
ac\cos\beta & 0 & c^2
\end{pmatrix}.
\]

The built-in correspondence matrices are

\[
C_{A\to M}=\begin{pmatrix}
0&1&-1\\
0&1&1\\
1&0&0
\end{pmatrix},\qquad
C_{M\to A}=C_{A\to M}^{-1}.
\]

The code keeps the active/passive coordinate convention fixed throughout the specialized engine. Tests detect accidental transposes/inversions because the benchmark habit planes, SMC vectors and twin shears change immediately if the convention is altered.

## 2. Compatibility by Metric Correspondence (CMC)

For an austenite direction \(u_A\), correspondence gives a product-phase crystallographic direction. The change in squared norm can be written

\[
u_A^T\,CMC\,u_A,
\]

with

\[
CMC=C_{M\to A}^{T}M_M C_{M\to A}-M_A.
\]

Therefore

\[
u_A^T CMC\,u_A=0
\]

is the zero-length-change quadratic form. Its solution set is generally a homogeneous quadratic cone.

Let \(q_1,q_2,q_3\) be the eigenvalues of the symmetric CMC matrix. The first-order A/M compatibility condition is implemented as

\[
q_i=0,\qquad q_jq_k\le 0,
\]

with the two nonzero eigenvalues of opposite sign. The quadratic form then factorizes into two linear plane equations, yielding two habit-plane candidates. Two zero eigenvalues give a second-order degeneracy (one plane); three zero eigenvalues give complete metric coincidence.

### Numerical CMC distance

For screening away from exact degeneracy the app reports

\[
d_{CMC}=\frac{\min_i |q_i|}{\max(1,\max_i|q_i|)}.
\]

This is a numerical distance-to-zero-eigenvalue diagnostic, **not a theorem that a finite threshold is physically universal**.

## 3. Shear by Metric Correspondence (SMC)

After a CMC habit plane \(m_A\) has been identified, the specialized metric formulation computes

\[
SMC=M_A^{-1}-C_{A\to M}M_M^{-1}C_{A\to M}^{T}
\]

and

\[
d_A=SMC\,m_A.
\]

The plane is a reciprocal-space vector; normalization uses the reciprocal metric when a normalized shear/shear comparison is required.

## 4. Metric transformation twins

For Type-I twins, a parent mirror operation is mapped through the correspondence. For Type-II twins, a parent 180° rotation is mapped through the correspondence. Twin shear amplitudes and the complementary twin elements are calculated directly from the product metric and intercorrespondence.

The metric shear/shear mismatch is

\[
\varepsilon=
\frac{\left\|2(m_A^T n)d_A-a\right\|_{M_A}}{s},
\]

where \(n\) is the direct-space normal of the twin shear plane, \(a\) is the twin shear vector and \(s=\|a\|\) is the twin shear amplitude under the adopted normalization.

Exact metric intercompatibility corresponds to

\[
\varepsilon=0.
\]

The app never converts a nonzero \(\varepsilon\) into a categorical “exact” result without an explicit numerical tolerance.

## 5. Stretch tensor construction

The classical cofactor route is calculated independently of the CMC/SMC route.

The metrics and correspondence satisfy

\[
U^T M_A U=C_{M\to A}^{T}M_M C_{M\to A}.
\]

For the normalized cubic parent, \(M_A=I\), so

\[
U^2=C_{M\to A}^{T}M_M C_{M\to A}.
\]

The code takes the principal symmetric positive-definite square root and diagonalizes \(U\):

\[
0<\lambda_1\le\lambda_2\le\lambda_3.
\]

For a general positive-definite parent metric the implementation first maps the target Gram matrix into the parent orthonormalized frame using \(M_A^{-1/2}\). Classical cofactor equations are evaluated in that orthonormalized frame.

## 6. Classical Type-I / Type-II domains and cofactor conditions

A parent two-fold axis \(\hat e\) defines

\[
Q=-I+2\hat e\otimes\hat e,
\]

and a symmetry-related stretch

\[
\hat U=Q U Q^T.
\]

For non-trivial domains, the conventional rank-one vectors are implemented as follows.

### Type I

\[
n_I=\hat e,
\]

\[
a_I=2\left(\frac{U^{-1}\hat e}{|U^{-1}\hat e|^2}-U\hat e\right).
\]

### Type II

\[
a_{II}=U\hat e,
\]

\[
n_{II}=2\left(\hat e-\frac{U^2\hat e}{|U\hat e|^2}\right).
\]

### Compound domains

If two perpendicular non-parallel generators \(\hat e_1,\hat e_2\) produce the same non-trivial related stretch, the pair is classified as a **Compound domain** rather than being double-counted as four separate Type-I/Type-II labels. The two rank-one solutions are

\[
n_C^1=\hat e_1,\qquad a_C^1=\xi U\hat e_2,
\]

\[
\xi=2\frac{\hat e_2\cdot U^{-2}\hat e_1}{\hat e_1\cdot U^{-2}\hat e_1},
\]

\[
n_C^2=\hat e_2,\qquad a_C^2=\eta U\hat e_1,
\]

\[
\eta=-2\frac{\hat e_2\cdot U^2\hat e_1}{\hat e_1\cdot U^2\hat e_1}.
\]

For Compound domains the software evaluates the **full CC2 equation directly**. It does not display the simplified Type-I/Type-II \(|U^{-1}\hat e|=1\) or \(|U\hat e|=1\) tests, because those simplified forms assume a non-compound Type-I/Type-II domain classification.

The app evaluates the full nonlinear cofactor conditions:

\[
CC1:\quad \lambda_2=1,
\]

\[
CC2:\quad a\cdot U\,\operatorname{cof}(U^2-I)n=0,
\]

\[
CC3:\quad \operatorname{tr}(U^2)-\det(U^2)-\frac{|a|^2|n|^2}{4}-2\ge0.
\]

A raw CC2 scalar and a scale-normalized residual are both retained. The normalized quantity is

\[
r_{CC2}=\frac{|CC2|}{\|a\|\,\|U\|_2\,\|\operatorname{cof}(U^2-I)\|_2\,\|n\|},
\]

with a small numerical denominator floor only to prevent division by machine zero.

When CC1 is satisfied, the familiar simplified cross-checks are also evaluated:

\[
\text{Type I: }\;|U^{-1}\hat e|=1,
\qquad
\text{Type II: }\;|U\hat e|=1.
\]

The full CC2 expression remains the certification quantity in the software.

For Compound domains, satisfaction of CC1–CC3 should not be conflated with the stronger statement that the transition layer between austenite and a compound laminate is eliminated. Current extreme-compatibility theory treats that stronger question separately. The app therefore reports Compound cofactor satisfaction in the PTMC workspace and reserves transition-layer-free compound diagnostics for the guarded frontier workspace.

## 7. Why the two compatibility routes are not silently merged

The software distinguishes:

- **CMC/SMC + metric shear/shear construction**, and
- **classical cofactor conditions**.

The 2026 metric/correspondence formulation demonstrates that its A/M, M/M and shear/shear conditions imply the classical supercompatibility conditions. It does not provide a proof of the full converse. Consequently, the app treats a case where CC1–CC3 pass but the specialized metric \(\varepsilon\) is nonzero as a **cross-framework discrepancy to inspect**, not as a software error to conceal.

This distinction is especially important in the C1-compatible benchmark encoded in the tests.

## 8. Symmetry reduction and full twin explorer

The parent cubic \(m\bar 3m\) point group is represented as the 48 signed permutation matrices.

The product monoclinic 2/m group is mapped by correspondence. The intersection subgroup

\[
H_C^A=G_A\cap C_{A\to M}G_M C_{M\to A}
\]

contains 4 operations for the built-in correspondence.

The parent point group is partitioned into double cosets

\[
H_C^A g H_C^A.
\]

The current implementation obtains 7 double-coset classes with sizes

\[
4,4,8,8,8,8,8.
\]

Two-fold-generating mirrors and 180° rotations are identified inside the classes and used to calculate all specialized Type-I/Type-II twin candidates while preserving their equivalence labels.

This avoids a hand-selected twin list in the research explorer while retaining a compact benchmark subset in the original validation layer.

## 9. Compatibility-distance dashboard

The app reports individual residuals rather than a composite percentage:

- \(|\lambda_2-1|\)
- CMC relative zero distance
- normalized CC2 residual
- CC3 margin
- metric shear/shear \(\varepsilon\), when CMC habit planes exist
- analytic C1 equality/inequality distances for the built-in transformation

No universal scalar combination is claimed to be physically privileged.

## 10. Temperature sweep

A temperature series supplies

\[
T,\;a_{B2}(T),\;a_{B19'}(T),\;b_{B19'}(T),\;c_{B19'}(T),\;\beta(T).
\]

Every row is independently reevaluated. The software does not interpolate unless the user has already supplied interpolated/fitted lattice parameters in the CSV.

This is intentionally a **compatibility-vs-temperature evaluator**, not a thermodynamic phase-stability model.

## 11. Uncertainty propagation

Input measurement uncertainties are modeled as independent Gaussian standard deviations unless the code is extended by the user. Monte Carlo samples are rejected if they produce non-positive lengths or angles outside \((0,180^\circ)\).

For each valid draw the chosen cofactor domain system and CMC metric are evaluated. The software reports fractions inside explicit tolerances for CC1, CC2, CC3 and CMC.

Important limitation: correlated refinement errors are common in diffraction. If covariance information is available, an independent-Gaussian model should be replaced by multivariate sampling before publication-level uncertainty claims are made.

## 12. Inverse design

The single-objective optimizer uses differential evolution. Its objective is

\[
J=w_1|\lambda_2-1|+w_2r_{CC2}+w_3\max(0,-CC3)+w_4d_{CMC}+w_5d_{prox},
\]

where \(d_{prox}\) is an optional scaled Euclidean distance from the starting normalized lattice geometry.

Every term and weight is exposed in the UI. The returned optimum is therefore reproducible and interpretable; it is not a hidden ML score.

The Pareto workbench uses Latin-hypercube sampling and flags nondominated points in the objective space \((|\lambda_2-1|,r_{CC2},CC3\ penalty,d_{CMC})\).

## 13. Composition → lattice machine learning

The ML module requires user-supplied training data. Target variables are

- `a_B2_A`
- `a_B19p_A`
- `b_B19p_A`
- `c_B19p_A`
- `beta_deg`

and the user selects numeric composition/processing descriptors as features.

Available regressors:

- linear regression with feature standardization
- random forest
- extra trees

K-fold cross-validation reports MAE, RMSE and \(R^2\) separately for each target. Only after lattice prediction are candidates passed through the exact physics screen.

No claim of causal composition→lattice mapping is made from model fit alone.

## 14. Multi-step transformations

The generic stage engine constructs the triclinic metric

\[
M=\begin{pmatrix}
a^2 & ab\cos\gamma & ac\cos\beta\\
ab\cos\gamma & b^2 & bc\cos\alpha\\
ac\cos\beta & bc\cos\alpha & c^2
\end{pmatrix}
\]

for parent and product phases and accepts a user-defined 3×3 correspondence matrix for every stage.

It reports the CMC distance and stretch eigenvalues stage by stage. It does **not** call specialized B2→B19′ twin formulas for arbitrary phase pairs.

## 15. Frontier compatibility guard

Recent extreme-compatibility theory provides exact target structures for particular transformation classes, including cubic→monoclinic-II. Those results are not universal tensors to compare against every monoclinic transformation.

The app first checks how the product monoclinic symmetry axis maps into the parent basis. Under the built-in correspondence the B19′ unique b axis maps to a cubic `<110>` direction, not `<100>`. Therefore the stored cubic→monoclinic-II exact target tensors are shown only as literature/reference objects and are **not used to certify the built-in transformation**.

The generic commutator diagnostic

\[
\|U_iU_j-U_jU_i\|_F
\]

is still evaluated among unique parent-symmetry stretch variants.

## 16. Numerical tolerances

Exact theory uses exact equalities and inequalities. Numerical implementation requires thresholds because of floating-point arithmetic and rounded experimental data.

The UI exposes CC1 and CC2 thresholds. Other internal tolerances are documented in function signatures and tests. A publication should report the chosen thresholds and, where experimental data are used, perform a sensitivity/uncertainty analysis rather than claiming an exact equality from rounded values.

## 16. Parent ↔ daughter orientation reconstruction

### 16.1 Orientation convention

Every grain orientation is a proper rotation matrix \(G\in SO(3)\) mapping crystal-frame coordinates into the specimen frame. The orientation relationship (OR) is stored as a child→parent crystal-frame rotation \(R_{cp}\). For a parent orientation \(G_p\), ideal daughter variants are generated as

\[
G_c^{(i)} = G_p S_p^{(i)} R_{cp},
\]

with \(S_p^{(i)}\) drawn from the proper parent point group and duplicate daughter orientations removed modulo the daughter point group. For a measured daughter orientation \(G_c\), candidate parent orientations are generated by the inverse relation, including daughter symmetry equivalents, and are deduplicated modulo parent symmetry.

The built-in proper rotation groups contain 24 cubic, 12 hexagonal, 4 orthorhombic, 2 monoclinic and 1 triclinic rotations. Improper rotations/reflections are not used as orientation matrices.

### 16.2 Built-in orientation relationships

The reconstruction workbench includes representative members of five standard OR families:

- Kurdjumov–Sachs, FCC→BCC: \(\{111\}_p\parallel\{110\}_c\), \(\langle110\rangle_p\parallel\langle111\rangle_c\);
- Nishiyama–Wassermann, FCC→BCC: \(\{111\}_p\parallel\{110\}_c\), \(\langle011\rangle_p\parallel\langle001\rangle_c\);
- Bain, FCC→BCC: \(\{100\}_p\parallel\{100\}_c\), \(\langle100\rangle_p\parallel\langle110\rangle_c\);
- Pitsch, FCC→BCC: \(\{100\}_p\parallel\{110\}_c\), \(\langle110\rangle_p\parallel\langle111\rangle_c\);
- Burgers, BCC→HCP: \(\{110\}_p\parallel(0001)_c\), \(\langle111\rangle_p\parallel\langle11\bar20\rangle_c\).

The exact representative matrix is less important than the OR family because crystallographically equivalent representatives are generated by the point groups. The automated tests require 24 KS, 12 NW, 3 Bain, 12 Pitsch and 12 Burgers daughter variants under the conventions above.

### 16.3 Reconstruction methods

The software intentionally provides several independent reconstruction routes.

**Neighbor voting.** Every daughter grain has a discrete set of candidate parent orientations. A candidate receives support from adjacent grains according to the closest parent-symmetry-reduced angular mismatch. A Gaussian-like angular weight is summed over neighbors, the best-supported candidate is retained, and selected candidates are merged spatially when their parent misorientation is below the chosen merge tolerance.

**Grain graph + Markov clustering.** A graph node represents a daughter grain. An edge weight is determined by the smallest angular mismatch between any parent candidates of the two adjacent daughter grains. The weighted graph is partitioned with Markov clustering (MCL). The inflation parameter is exposed because it changes cluster granularity and must therefore be reported in a publication.

**Variant graph.** A graph state is retained for every candidate-parent orientation of every daughter grain rather than collapsing a grain to one candidate at the outset. Neighbor candidate-pair compatibility is converted to weights, probability-like support is propagated iteratively, and only then is the best candidate selected. This preserves ambiguity information that a one-node-per-grain graph loses.

**Nucleation + growth.** Parent grains are nucleated from daughter grains with strong local OR support using a strict angular criterion. Unassigned neighbors are then grown into the existing parents with a larger tolerance. The two angular tolerances are independent inputs and are exported with the result.

**Operator / groupoid consistency.** The ideal daughter–daughter operators induced by the parent symmetry and the OR are calculated first. A measured neighboring daughter pair is eligible only if its misorientation is close to one of these theoretical operators modulo daughter symmetry. Eligible edges are then required to admit a mutually consistent parent candidate before they are clustered. This is deliberately stricter than chaining all locally operator-like edges.

### 16.4 Orientation-relationship refinement

A supplied physical OR can be refined by a bounded rotation-vector correction. For every trial correction, the candidate-parent sets are recomputed and neighboring grains are scored by the minimum common-parent angular disagreement. A pseudo-Huber objective reduces sensitivity to boundary outliers. The correction is explicitly bounded in degrees. The routine is intended to refine a known OR; it is not an unconstrained OR-discovery algorithm.

A publication should report the starting OR, correction bound, objective before/after refinement and final correction angle.

The per-grain `confidence` field in reconstruction exports is a **heuristic support score**, constructed from the absolute parent-candidate fit and separation from the next candidate. It is bounded between 0 and 1 for convenient ranking but is not a calibrated statistical probability; papers should call it a support score unless an external calibration is performed.

### 16.5 Important reconstruction limitations

- Grain-boundary adjacency is preferred. A centroid k-nearest-neighbor graph is provided only as an explicitly labelled approximation.
- Different reconstruction methods can disagree at prior-parent boundaries, twins, sparse-variant regions and highly deformed regions; the app therefore supports method comparison instead of concealing disagreement.
- A reconstruction from orientations alone can be intrinsically ambiguous when distinct parent states share the same observed daughter variants.
- The implementations are transparent research implementations of established algorithmic families; they are not binary reproductions of external packages.

## 17. Compatibility methods beyond the cofactor checklist

### 17.1 Hadamard jump condition

Two deformation gradients \(A\) and \(B\) can meet across a coherent planar interface if

\[
A-B=a\otimes m.
\]

The matrix jump is therefore rank one (rank zero is the trivial identical-state case). The generic diagnostic computes the singular values \(\sigma_1\ge\sigma_2\ge\sigma_3\) of \(A-B\) and reports \(\sigma_2/\sigma_1\) as a scale-free numerical rank-one residual.

### 17.2 Single-variant IPS spectral condition

For a positive-definite transformation stretch \(U\), the classical single-variant austenite/martensite invariant-plane condition is

\[
\lambda_1\le1,\qquad \lambda_2=1,\qquad \lambda_3\ge1.
\]

The app reports the ordered principal stretches and \(|\lambda_2-1|\) separately from all twin and cofactor diagnostics.

### 17.3 Pairwise martensite-variant rank-one compatibility

For two stretch variants \(U_i,U_j\), define

\[
B_{ij}=U_j^{-T}U_i^TU_iU_j^{-1}.
\]

If its ordered eigenvalues are \(\chi_1\le\chi_2\le\chi_3\), the rank-one compatibility problem admits a twinning solution when

\[
\chi_1\le1,\qquad\chi_2=1,\qquad\chi_3\ge1.
\]

The pairwise explorer evaluates this criterion for every symmetry-generated stretch-variant pair and exports the residual \(|\chi_2-1|\) and bracket margins.

### 17.4 Direct all-volume-fraction laminate verification

For a twin relation represented by rank-one vectors \(a,n\), the average laminate gradient is sampled as

\[
F(f)=U+f\,a\otimes n,\qquad 0\le f\le1.
\]

At every sampled volume fraction the singular values of \(F(f)\) are computed. Under the full cofactor conditions, the middle singular value remains one for all \(f\). The app therefore reports both the maximum and RMS values of \(|\sigma_2(F(f))-1|\). This is a direct numerical cross-check of the all-volume-fraction implication rather than a re-evaluation of the algebraic CC equations.

### 17.5 Cubic→orthorhombic triplet-condition diagnostic

The triplet condition is a distinct martensitic supercompatibility principle involving three variants. The current implementation is intentionally restricted to the cubic→orthorhombic principal-stretch parameterization

\[
\alpha=\frac{a_o}{a_c},\qquad
\beta=\frac{b_o}{\sqrt2a_c},\qquad
\gamma=\frac{c_o}{\sqrt2a_c}.
\]

It reports residuals for the specialized TC-I and TC-II algebraic conditions implemented in `src/compatibility_methods.py`. The general theory contains additional triplet branches; this module therefore does **not** claim to be a complete general triplet solver. It is guarded in the UI: it is **not** used as a certification test for the built-in cubic→monoclinic transformation. The 2023 corrigendum to the original triplet-condition paper is listed in the references and should be cited with the primary paper when these results are used academically.

## 18. Cross-method interpretation

The app deliberately does not create a universal scalar “compatibility percentage”. Each method asks a different mathematical question:

| Method | Question |
|---|---|
| CMC degeneracy | Does correspondence preserve an entire plane of lengths? |
| SMC + shear/shear | Does the A/M IPS shear interlock with a specific M/M twin shear? |
| \(\lambda_2=1\) | Can one stretch variant form an invariant-plane interface after rotation? |
| pairwise \(\chi_2=1\) | Can two stretch variants be rank-one connected? |
| CC1–CC3 | Does a selected twin system admit transition-layer-free A/M compatibility for every twin fraction? |
| direct all-\(f\) scan | Does the middle singular value numerically remain one for every sampled twin fraction? |
| triplet condition | Can a compatible three-variant martensitic junction exist in the supported transformation class? |
| Hadamard rank-one test | Is a supplied deformation-gradient jump planar-compatible? |

A paper should report the individual residuals and applicability assumptions rather than replacing them by an invented score.
