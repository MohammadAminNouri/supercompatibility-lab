# Notation and symbol guide

This project follows a **no-symbol-left-behind** rule: a reader should never have to guess what a mathematical letter, subscript, superscript, operator, crystallographic bracket, unit, acronym, input field, or technical output means. The Streamlit UI shows local definitions beside formulas/results; this file is the consolidated manuscript/code-review reference.

## Acronyms

- **PTMC** = phenomenological theory of martensite crystallography.
- **CC1 / CC2 / CC3** = first / second / third cofactor conditions.
- **CMC** = compatibility metric-change matrix.
- **SMC** = shear metric-change matrix.
- **IPS** = invariant-plane strain.
- **OR** = orientation relationship.
- **A/M** = austenite/martensite interface; **M/M** = martensite/martensite interface.
- **XRD** = X-ray diffraction; **EBSD** = electron backscatter diffraction.
- **ML** = machine learning; **RMS** = root-mean-square.

## Primary input symbols

| Symbol | Meaning | Unit/scale | Note |
|---|---|---|---|
| `a_{B2}` | B2 austenite cubic unit-cell edge length | Å | Å (ångström) = 10⁻¹⁰ m |
| `a_{B19'}` | B19′ martensite a-axis unit-cell length | Å | physical measured lattice parameter |
| `b_{B19'}` | B19′ martensite b-axis unit-cell length | Å | physical measured lattice parameter |
| `c_{B19'}` | B19′ martensite c-axis unit-cell length | Å | physical measured lattice parameter |
| `\beta` | monoclinic angle between the B19′ a and c axes | ° | degrees in the UI; converted internally to radians for trigonometric functions |
| `\mathrm{tol}_{CC1}` | numerical decision tolerance applied to \|λ₂−1\| | dimensionless | not a physical law; must be reported in a paper |
| `\mathrm{tol}_{CC2}` | numerical decision tolerance applied to the normalized CC2 residual | dimensionless | not a physical law; must be reported in a paper |

## Canonical formula symbols

| Symbol | Meaning | Unit/scale | Note |
|---|---|---|---|
| `\lambda_1` | smallest principal stretch factor of the transformation stretch tensor U | dimensionless | <1 means contraction along its principal direction |
| `\lambda_2` | middle principal stretch factor of U | dimensionless | \lambda_2=1 is the classical single-variant invariant-plane compatibility condition |
| `\lambda_3` | largest principal stretch factor of U | dimensionless | >1 means extension along its principal direction |
| `U` | symmetric positive-definite transformation stretch tensor | dimensionless | pure stretch part of the deformation gradient in PTMC |
| `I` | identity matrix | dimensionless |  |
| `\mathrm{cof}(\cdot)` | matrix cofactor operator | same matrix power as its argument |  |
| `\mathrm{tr}(\cdot)` | matrix trace: sum of diagonal entries | dimensionless here |  |
| `\det(\cdot)` | matrix determinant | dimensionless here | for U it measures transformation volume ratio |
| `a` | rank-one twin shear/displacement vector in the cofactor or twin relation | dimensionless in deformation-gradient form | Do not confuse with normalized lattice ratio a in the NiTi analytical section |
| `n` | normal/covector of the martensite twin or domain plane | reciprocal-direction scale | only direction matters; normalization convention is recorded |
| `e_1` | principal direction (eigenvector) associated with \lambda_1 | unit vector |  |
| `e_2` | principal direction (eigenvector) associated with \lambda_2 | unit vector |  |
| `e_3` | principal direction (eigenvector) associated with \lambda_3 | unit vector |  |
| `\hat e` | unit symmetry axis used to classify a Type-I or Type-II cofactor domain | unit vector |  |
| `v` | direction whose length is unchanged by U in the PTMC construction | direction vector |  |
| `v'` | rotated image of v used to construct the IPS rotation | direction vector |  |
| `R` | proper rotation matrix used with U to form the deformation gradient F | dimensionless |  |
| `F` | deformation gradient mapping the parent configuration to the transformed configuration | dimensionless |  |
| `d` | IPS displacement/shear vector | dimensionless strain-like vector |  |
| `m` | normal/covector of the invariant habit plane in F=I+d m^T | reciprocal-direction scale |  |
| `m_A` | habit-plane normal expressed in the austenite/parent crystallographic basis | reciprocal-direction scale |  |
| `m_B` | habit-plane normal of the second compatible martensite variant | reciprocal-direction scale |  |
| `\tau` | magnitude of the tangential shear part of d | dimensionless strain |  |
| `\tau_\perp` | component of IPS shear perpendicular to the line of intersection used in the geometric construction | dimensionless strain |  |
| `\delta` | normal dilatation across the habit plane; 1+\delta is the normal stretch contribution | dimensionless strain |  |
| `\phi` | angle between the individual habit plane and the relevant twin mirror/junction plane | radian or degree when displayed |  |
| `P` | volume-fraction-averaged laminate deformation gradient | dimensionless |  |
| `f` | volume fraction of martensite variant 1 in a two-variant laminate | dimensionless, 0 to 1 |  |
| `F_1` | deformation gradient of martensite variant 1 | dimensionless |  |
| `F_2` | deformation gradient of martensite variant 2 | dimensionless |  |
| `G^A` | point-symmetry group of the austenite/parent phase | finite group |  |
| `G^M` | point-symmetry group of the martensite/daughter phase | finite group |  |
| `H_C^A` | correspondence subgroup: parent symmetries preserved by the chosen lattice correspondence | finite subgroup of G^A |  |
| `C^{A\to M}` | correspondence matrix mapping parent crystallographic coordinates to daughter/martensite coordinates | dimensionless matrix |  |
| `C^{M\to A}` | inverse correspondence matrix mapping daughter/martensite crystallographic coordinates to parent coordinates | dimensionless matrix |  |
| `C^{A\to M_i}` | correspondence matrix for martensite correspondence variant i | dimensionless matrix |  |
| `g_i^A` | representative parent symmetry operation for coset or double-coset class i | orthogonal/symmetry matrix in crystal metric |  |
| `N_C^M` | number of distinct martensite correspondence variants generated from the parent | count |  |
| `\|G\|` | order/cardinality of finite group G: number of symmetry operations | count |  |
| `O_k` | k-th double-coset/intercorrespondence class | set of symmetry operations |  |
| `M_A` | metric tensor of the austenite/parent lattice | length squared; normalized form is dimensionless |  |
| `M_M` | metric tensor of the martensite/daughter lattice | length squared; normalized form is dimensionless |  |
| `CMC` | compatibility metric-change matrix comparing squared lengths before and after correspondence | length squared; normalized form dimensionless |  |
| `u_A` | parent/austenite crystallographic direction vector being tested for preserved length | direction coordinates |  |
| `u` | generic direction vector | direction coordinates |  |
| `q(\cdot)` | quadratic form associated with CMC, usually u^T CMC u | length squared; normalized form dimensionless |  |
| `S_{CMC}` | zero set of the CMC quadratic form | set of directions |  |
| `G_{CMC}` | parent symmetries that preserve the CMC zero set | finite group/subgroup |  |
| `q_i,q_j,q_k` | eigenvalues of the symmetric CMC matrix | same scale as CMC | indices i,j,k denote three distinct eigenvalues |
| `SMC` | shear metric-change matrix used to obtain the IPS displacement vector from a compatible habit plane | inverse metric scale; normalized implementation dimensionless |  |
| `d_A` | austenite-basis IPS displacement/shear vector | dimensionless strain-like vector |  |
| `\varepsilon` | normalized mismatch between the A/M IPS shear requirement and an M/M twin shear | dimensionless | 0 is exact matching |
| `s` | magnitude (amplitude) of the martensite twin shear vector a | dimensionless shear |  |
| `a_{B2}` | physical cubic B2 austenite lattice parameter (unit-cell edge length) | Å |  |
| `a` | normalized B19′ a-axis ratio a_B19′/a_B2 | dimensionless | not chemical composition and not the twin-shear vector a |
| `b` | normalized B19′ b-axis ratio b_B19′/a_B2 | dimensionless |  |
| `c` | normalized B19′ c-axis ratio c_B19′/a_B2 | dimensionless |  |
| `\beta` | monoclinic B19′ angle between the a and c lattice axes | degree in inputs; converted to radian internally |  |
| `K` | analytical abbreviation K=2a²+c²−4 used in the closed-form CMC eigenvalues | dimensionless |  |
| `\Delta` | analytical discriminant 4a⁴+c⁴+4a²c²cos(2β) | dimensionless |  |
| `q_1` | CMC eigenvalue associated with the normalized b-ratio branch | dimensionless in normalized CMC |  |
| `q_2,q_3` | the other two analytical CMC eigenvalues | dimensionless in normalized CMC |  |
| `\sqrt2` | square root of 2, approximately 1.41421356 | dimensionless |  |
| `T^{M\to A}` | passive orientation/coordinate-transformation matrix from martensite coordinates to austenite coordinates | dimensionless matrix |  |
| `F_A` | active distortion matrix of the austenite/parent lattice | dimensionless |  |
| `p_M` | martensite twin-plane covector/normal | reciprocal direction |  |
| `p_A` | parent symmetry-plane covector/normal | reciprocal direction |  |
| `C_{int}` | intercorrespondence operator relating a pair of martensite correspondence variants | dimensionless matrix |  |
| `a_M` | martensite Type-II shear/rotation direction inherited by correspondence | direct-lattice direction |  |
| `a_A` | parent 180° symmetry axis/direction used to generate a Type-II twin | direct-lattice direction |  |
| `T` | temperature | K or °C exactly as supplied in the dataset |  |
| `\mathcal{R}(T)` | vector of calculated compatibility results evaluated at temperature T | mixed outputs |  |
| `x^{(k)}` | k-th Monte-Carlo sampled input vector | same units as sampled inputs |  |
| `\mu` | mean vector of uncertain measured inputs | same units as inputs |  |
| `\Sigma` | input covariance matrix | products of input units |  |
| `y^{(k)}` | calculated output for Monte-Carlo sample k | depends on output |  |
| `x` | vector of design variables optimized by the inverse-design solver | mixed; explicitly listed with bounds |  |
| `J(x)` | scalar optimization objective assembled from chosen compatibility residuals/weights | dimensionless after normalization |  |
| `\hat y` | machine-learning-predicted lattice/output vector | target-specific units |  |
| `f_\theta` | trained machine-learning mapping parameterized by learned parameters θ | data-driven function |  |
| `\theta` | learned model parameters | model dependent |  |
| `A` | first deformation gradient/matrix in the Hadamard jump test | dimensionless |  |
| `B` | second deformation gradient/matrix in the Hadamard jump test | dimensionless |  |
| `m` | interface normal/covector in the Hadamard rank-one jump A−B=a⊗m | reciprocal-direction scale |  |
| `\chi_1,\chi_2,\chi_3` | ordered eigenvalues of the pairwise twin-compatibility spectral matrix χ | dimensionless |  |
| `U_i` | stretch tensor of martensite variant i | dimensionless |  |
| `U_j` | stretch tensor of martensite variant j | dimensionless |  |
| `\sigma_2` | middle singular value of a deformation gradient | dimensionless |  |
| `\alpha` | normalized principal stretch/lattice parameter used by the specialized cubic→orthorhombic triplet condition | dimensionless |  |
| `\beta` | second normalized principal stretch/lattice parameter in the specialized triplet condition | dimensionless | not the monoclinic angle β in the B19′ model |
| `\gamma` | third normalized principal stretch/lattice parameter in the specialized triplet condition | dimensionless |  |
| `g_P` | parent crystal orientation matrix | rotation matrix |  |
| `g_D` | daughter crystal orientation matrix | rotation matrix |  |
| `OR` | orientation relationship mapping parent and daughter crystallographic frames | rotation/orientation relation |  |
| `variant` | symmetry-equivalent crystallographic transformation choice | discrete index/class |  |

## Core output symbols

| Symbol | Meaning | Unit/scale | Note |
|---|---|---|---|
| `\lambda_2` | middle principal stretch factor of U | dimensionless | \lambda_2=1 is the classical single-variant invariant-plane compatibility condition |
| `d_{CMC}` | relative numerical distance of the CMC eigenspectrum from exact degeneracy | dimensionless | software residual; not identical to the paper-defined family distances |
| `r_{CC2}` | scale-normalized absolute residual of the full CC2 equality for the best tested domain system | dimensionless | 0 is exact |
| `\varepsilon` | normalized mismatch between the A/M IPS shear requirement and an M/M twin shear | dimensionless | 0 is exact matching |
| `\det U` | determinant of U; transformation volume ratio in the stretch description | dimensionless | 1 means zero net volume change |
| `q_i,q_j,q_k` | eigenvalues of the symmetric CMC matrix | same scale as CMC | indices i,j,k denote three distinct eigenvalues |
| `m_A=(h,k,l)` | compatible habit-plane normal/covector in parent crystallographic coordinates | reciprocal direction | components are directional coordinates; overall scaling does not change the plane |
| `d_A` | austenite-basis IPS displacement/shear vector | dimensionless strain-like vector |  |
| `s` | magnitude (amplitude) of the martensite twin shear vector a | dimensionless shear |  |

## Formula operators and notation

| Symbol | Meaning | Unit/scale | Note |
|---|---|---|---|
| `(\cdot)^T` | transpose: swap matrix/vector rows and columns | operator |  |
| `(\cdot)^{-1}` | matrix inverse: the matrix that undoes the original matrix | operator |  |
| `(\cdot)^{-T}` | inverse transpose: transpose of the matrix inverse | operator |  |
| `\cdot` | dot/inner product: combines two vectors into a scalar | operator |  |
| `\otimes` | outer/dyadic product: combines two vectors into a rank-one matrix | operator |  |
| `\mathrm{cof}(\cdot)` | cofactor-matrix operator | operator |  |
| `\mathrm{tr}(\cdot)` | trace: sum of matrix diagonal entries | operator |  |
| `\det(\cdot)` | determinant: scalar matrix volume/scale factor | operator |  |
| `\|x\|\;\text{or}\;\\|x\\|` | absolute value for a scalar or norm/length for a vector, as stated by the equation | operator |  |
| `\sqrt{x}` | positive square root of x | operator |  |
| `\pm` | two mathematical branches: one with plus and one with minus | operator |  |
| `\cap` | set intersection: elements common to both sets | set operator |  |
| `\forall` | for every value in the stated set | logical symbol |  |
| `\in` | is an element/member of | set relation |  |
| `\ge` | greater than or equal to | comparison |  |
| `\le` | less than or equal to | comparison |  |
| `\to` | maps/transforms from the object on the left to the object on the right | mapping symbol |  |
| `\parallel` | parallel to | geometric relation |  |
| `=` | equals: the left and right sides have the same value | relation |  |
| `+` | addition: add the quantities on each side of the sign | arithmetic operator |  |
| `-` | subtraction or a negative sign, according to its position | arithmetic operator |  |
| `/\;\text{or}\;\frac{a}{b}` | division: numerator divided by denominator | arithmetic operator |  |
| `x^p` | power/exponent: multiply x by itself according to exponent p; negative exponents mean reciprocal powers | operator |  |
| `x_i` | subscript: label/index identifying a component, phase, branch, or ordered value; it is not multiplication | notation |  |
| `AB` | adjacent compatible scalars/matrices/vectors mean multiplication or linear action in the order written | multiplication notation |  |
| `\sin` | sine trigonometric function of the stated angle | operator |  |
| `\cos` | cosine trigonometric function of the stated angle | operator |  |
| `\tan` | tangent trigonometric function of the stated angle | operator |  |
| `\mathrm{rank}(A)` | matrix rank: number of linearly independent directions carried by matrix A | operator |  |
| `\mathrm{eig}(A)` | eigenvalue operation: returns principal scalar factors associated with matrix A | operator |  |
| `\Longleftrightarrow` | if and only if: both statements imply each other | logical relation |  |
| `\sim` | is sampled/distributed according to the probability distribution on the right | statistical relation |  |
| `\mathcal{N}(\mu,\Sigma)` | multivariate normal (Gaussian) distribution with mean vector μ and covariance matrix Σ | probability distribution |  |
| `\min` | minimize: search for variable values that make the following objective as small as possible | optimization operator |  |
| `^\circ` | degree symbol: the angle is measured in degrees | angle unit |  |
| `\begin{pmatrix}\cdots\end{pmatrix}` | matrix: rectangular array; rows and columns refer to the axes stated beside the output | notation |  |
| `\{x: condition\}` | set-builder notation: the set of x values satisfying the condition after the colon | set notation |  |
| `[a,b]` | closed interval: every value from a through b, including both endpoints | set notation |  |

## Direction and plane notation

- `[u v w]` denotes a **direct-lattice direction**. Its components are relative crystallographic direction coordinates. Multiplying all components by the same nonzero scalar does not change the direction.
- `(h k l)` denotes a **crystallographic plane** through its reciprocal-lattice normal/covector. Multiplying all indices by the same nonzero scalar describes the same plane.
- `A` in a phase superscript/subscript means austenite/parent phase; `M` means martensite/daughter phase.
- A superscript `T` on a vector/matrix means transpose; plain scalar `T` in the temperature workspace means temperature.

## Units

- `Å` (ångström) = `10^-10 m`, used for lattice parameters.
- `°` = degree, used for user-facing angular inputs/outputs unless radians are explicitly stated.
- `K` = kelvin when attached to temperature.
- `MPa` = megapascal = `10^6 Pa`, used for stress.
- `m^-2` is used for dislocation density.
- Principal stretches, normalized lattice ratios, most compatibility residuals and deformation-gradient shear amplitudes are dimensionless.

## Ambiguous letters

- `a` can be the normalized lattice ratio `a_B19′/a_B2` **or** the rank-one/twin shear vector `a`. The UI defines the active meaning beside each relation.
- `β` is the B19′ monoclinic angle in the B2→B19′ model, while the specialized cubic→orthorhombic triplet workspace also uses a conventional dimensionless ratio called `β`; the UI explicitly distinguishes them.
- `T` can be temperature or part of the coordinate-transformation matrix `T^{M→A}`; the local formula/input definition states which one applies.

## Paper-writing rule

When copying a result into a manuscript, copy the associated equation/provenance block or export the paper-ready record. Do not quote a bare number without its symbol definition, unit/scale, source/theorem, residual/tolerance and applicability caveat.
