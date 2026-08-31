from __future__ import annotations

"""Human-readable notation dictionary used by both the Streamlit UI and exports.

The project deliberately treats notation as part of the scientific result.  Every
registered equation is mapped to explicit symbol definitions so a reader never
has to infer what a letter, superscript, subscript, operator, or unit means.
"""

from dataclasses import asdict, dataclass
from typing import Iterable


@dataclass(frozen=True)
class SymbolDefinition:
    symbol: str
    meaning: str
    unit: str = "dimensionless / context dependent"
    note: str = ""


# Canonical symbols.  Keys are stable machine identifiers; ``symbol`` is display text.
SYMBOLS: dict[str, SymbolDefinition] = {
    "lambda1": SymbolDefinition(r"\lambda_1", "smallest principal stretch factor of the transformation stretch tensor U", "dimensionless", "<1 means contraction along its principal direction"),
    "lambda2": SymbolDefinition(r"\lambda_2", "middle principal stretch factor of U", "dimensionless", r"\lambda_2=1 is the classical single-variant invariant-plane compatibility condition"),
    "lambda3": SymbolDefinition(r"\lambda_3", "largest principal stretch factor of U", "dimensionless", ">1 means extension along its principal direction"),
    "U": SymbolDefinition(r"U", "symmetric positive-definite transformation stretch tensor", "dimensionless", "pure stretch part of the deformation gradient in PTMC"),
    "I": SymbolDefinition(r"I", "identity matrix", "dimensionless"),
    "cof": SymbolDefinition(r"\mathrm{cof}(\cdot)", "matrix cofactor operator", "same matrix power as its argument"),
    "tr": SymbolDefinition(r"\mathrm{tr}(\cdot)", "matrix trace: sum of diagonal entries", "dimensionless here"),
    "det": SymbolDefinition(r"\det(\cdot)", "matrix determinant", "dimensionless here", "for U it measures transformation volume ratio"),
    "a_twin": SymbolDefinition(r"a", "rank-one twin shear/displacement vector in the cofactor or twin relation", "dimensionless in deformation-gradient form", "Do not confuse with normalized lattice ratio a in the NiTi analytical section"),
    "n_twin": SymbolDefinition(r"n", "normal/covector of the martensite twin or domain plane", "reciprocal-direction scale", "only direction matters; normalization convention is recorded"),
    "e1": SymbolDefinition(r"e_1", r"principal direction (eigenvector) associated with \lambda_1", "unit vector"),
    "e2": SymbolDefinition(r"e_2", r"principal direction (eigenvector) associated with \lambda_2", "unit vector"),
    "e3": SymbolDefinition(r"e_3", r"principal direction (eigenvector) associated with \lambda_3", "unit vector"),
    "ehat": SymbolDefinition(r"\hat e", "unit symmetry axis used to classify a Type-I or Type-II cofactor domain", "unit vector"),
    "v": SymbolDefinition(r"v", "direction whose length is unchanged by U in the PTMC construction", "direction vector"),
    "vprime": SymbolDefinition(r"v'", "rotated image of v used to construct the IPS rotation", "direction vector"),
    "R": SymbolDefinition(r"R", "proper rotation matrix used with U to form the deformation gradient F", "dimensionless"),
    "F": SymbolDefinition(r"F", "deformation gradient mapping the parent configuration to the transformed configuration", "dimensionless"),
    "d": SymbolDefinition(r"d", "IPS displacement/shear vector", "dimensionless strain-like vector"),
    "m": SymbolDefinition(r"m", "normal/covector of the invariant habit plane in F=I+d m^T", "reciprocal-direction scale"),
    "mA": SymbolDefinition(r"m_A", "habit-plane normal expressed in the austenite/parent crystallographic basis", "reciprocal-direction scale"),
    "mB": SymbolDefinition(r"m_B", "habit-plane normal of the second compatible martensite variant", "reciprocal-direction scale"),
    "tau": SymbolDefinition(r"\tau", "magnitude of the tangential shear part of d", "dimensionless strain"),
    "tauperp": SymbolDefinition(r"\tau_\perp", "component of IPS shear perpendicular to the line of intersection used in the geometric construction", "dimensionless strain"),
    "delta": SymbolDefinition(r"\delta", r"normal dilatation across the habit plane; 1+\delta is the normal stretch contribution", "dimensionless strain"),
    "phi": SymbolDefinition(r"\phi", "angle between the individual habit plane and the relevant twin mirror/junction plane", "radian or degree when displayed"),
    "P": SymbolDefinition(r"P", "volume-fraction-averaged laminate deformation gradient", "dimensionless"),
    "f": SymbolDefinition(r"f", "volume fraction of martensite variant 1 in a two-variant laminate", "dimensionless, 0 to 1"),
    "F1": SymbolDefinition(r"F_1", "deformation gradient of martensite variant 1", "dimensionless"),
    "F2": SymbolDefinition(r"F_2", "deformation gradient of martensite variant 2", "dimensionless"),
    "GA": SymbolDefinition(r"G^A", "point-symmetry group of the austenite/parent phase", "finite group"),
    "GM": SymbolDefinition(r"G^M", "point-symmetry group of the martensite/daughter phase", "finite group"),
    "HC": SymbolDefinition(r"H_C^A", "correspondence subgroup: parent symmetries preserved by the chosen lattice correspondence", "finite subgroup of G^A"),
    "C_A_M": SymbolDefinition(r"C^{A\to M}", "correspondence matrix mapping parent crystallographic coordinates to daughter/martensite coordinates", "dimensionless matrix"),
    "C_M_A": SymbolDefinition(r"C^{M\to A}", "inverse correspondence matrix mapping daughter/martensite crystallographic coordinates to parent coordinates", "dimensionless matrix"),
    "Ci": SymbolDefinition(r"C^{A\to M_i}", "correspondence matrix for martensite correspondence variant i", "dimensionless matrix"),
    "gi": SymbolDefinition(r"g_i^A", "representative parent symmetry operation for coset or double-coset class i", "orthogonal/symmetry matrix in crystal metric"),
    "NC": SymbolDefinition(r"N_C^M", "number of distinct martensite correspondence variants generated from the parent", "count"),
    "card": SymbolDefinition(r"|G|", "order/cardinality of finite group G: number of symmetry operations", "count"),
    "Ok": SymbolDefinition(r"O_k", "k-th double-coset/intercorrespondence class", "set of symmetry operations"),
    "MA": SymbolDefinition(r"M_A", "metric tensor of the austenite/parent lattice", "length squared; normalized form is dimensionless"),
    "MM": SymbolDefinition(r"M_M", "metric tensor of the martensite/daughter lattice", "length squared; normalized form is dimensionless"),
    "CMC": SymbolDefinition(r"CMC", "compatibility metric-change matrix comparing squared lengths before and after correspondence", "length squared; normalized form dimensionless"),
    "uA": SymbolDefinition(r"u_A", "parent/austenite crystallographic direction vector being tested for preserved length", "direction coordinates"),
    "u": SymbolDefinition(r"u", "generic direction vector", "direction coordinates"),
    "q": SymbolDefinition(r"q(\cdot)", "quadratic form associated with CMC, usually u^T CMC u", "length squared; normalized form dimensionless"),
    "SCMC": SymbolDefinition(r"S_{CMC}", "zero set of the CMC quadratic form", "set of directions"),
    "GCMC": SymbolDefinition(r"G_{CMC}", "parent symmetries that preserve the CMC zero set", "finite group/subgroup"),
    "qi": SymbolDefinition(r"q_i,q_j,q_k", "eigenvalues of the symmetric CMC matrix", "same scale as CMC", "indices i,j,k denote three distinct eigenvalues"),
    "SMC": SymbolDefinition(r"SMC", "shear metric-change matrix used to obtain the IPS displacement vector from a compatible habit plane", "inverse metric scale; normalized implementation dimensionless"),
    "dA": SymbolDefinition(r"d_A", "austenite-basis IPS displacement/shear vector", "dimensionless strain-like vector"),
    "epsilon": SymbolDefinition(r"\varepsilon", "normalized mismatch between the A/M IPS shear requirement and an M/M twin shear", "dimensionless", "0 is exact matching"),
    "s": SymbolDefinition(r"s", "magnitude (amplitude) of the martensite twin shear vector a", "dimensionless shear"),
    "aB2": SymbolDefinition(r"a_{B2}", "physical cubic B2 austenite lattice parameter (unit-cell edge length)", "Å"),
    "aRatio": SymbolDefinition(r"a", "normalized B19′ a-axis ratio a_B19′/a_B2", "dimensionless", "not chemical composition and not the twin-shear vector a"),
    "bRatio": SymbolDefinition(r"b", "normalized B19′ b-axis ratio b_B19′/a_B2", "dimensionless"),
    "cRatio": SymbolDefinition(r"c", "normalized B19′ c-axis ratio c_B19′/a_B2", "dimensionless"),
    "beta": SymbolDefinition(r"\beta", "monoclinic B19′ angle between the a and c lattice axes", "degree in inputs; converted to radian internally"),
    "K": SymbolDefinition(r"K", "analytical abbreviation K=2a²+c²−4 used in the closed-form CMC eigenvalues", "dimensionless"),
    "Delta": SymbolDefinition(r"\Delta", "analytical discriminant 4a⁴+c⁴+4a²c²cos(2β)", "dimensionless"),
    "q1": SymbolDefinition(r"q_1", "CMC eigenvalue associated with the normalized b-ratio branch", "dimensionless in normalized CMC"),
    "q23": SymbolDefinition(r"q_2,q_3", "the other two analytical CMC eigenvalues", "dimensionless in normalized CMC"),
    "sqrt2": SymbolDefinition(r"\sqrt2", "square root of 2, approximately 1.41421356", "dimensionless"),
    "T_M_A": SymbolDefinition(r"T^{M\to A}", "passive orientation/coordinate-transformation matrix from martensite coordinates to austenite coordinates", "dimensionless matrix"),
    "FA": SymbolDefinition(r"F_A", "active distortion matrix of the austenite/parent lattice", "dimensionless"),
    "pM": SymbolDefinition(r"p_M", "martensite twin-plane covector/normal", "reciprocal direction"),
    "pA": SymbolDefinition(r"p_A", "parent symmetry-plane covector/normal", "reciprocal direction"),
    "Cint": SymbolDefinition(r"C_{int}", "intercorrespondence operator relating a pair of martensite correspondence variants", "dimensionless matrix"),
    "aM": SymbolDefinition(r"a_M", "martensite Type-II shear/rotation direction inherited by correspondence", "direct-lattice direction"),
    "aA": SymbolDefinition(r"a_A", "parent 180° symmetry axis/direction used to generate a Type-II twin", "direct-lattice direction"),
    "T": SymbolDefinition(r"T", "temperature", "K or °C exactly as supplied in the dataset"),
    "RofT": SymbolDefinition(r"\mathcal{R}(T)", "vector of calculated compatibility results evaluated at temperature T", "mixed outputs"),
    "xk": SymbolDefinition(r"x^{(k)}", "k-th Monte-Carlo sampled input vector", "same units as sampled inputs"),
    "mu": SymbolDefinition(r"\mu", "mean vector of uncertain measured inputs", "same units as inputs"),
    "Sigma": SymbolDefinition(r"\Sigma", "input covariance matrix", "products of input units"),
    "yk": SymbolDefinition(r"y^{(k)}", "calculated output for Monte-Carlo sample k", "depends on output"),
    "xopt": SymbolDefinition(r"x", "vector of design variables optimized by the inverse-design solver", "mixed; explicitly listed with bounds"),
    "J": SymbolDefinition(r"J(x)", "scalar optimization objective assembled from chosen compatibility residuals/weights", "dimensionless after normalization"),
    "yhat": SymbolDefinition(r"\hat y", "machine-learning-predicted lattice/output vector", "target-specific units"),
    "ftheta": SymbolDefinition(r"f_\theta", "trained machine-learning mapping parameterized by learned parameters θ", "data-driven function"),
    "theta": SymbolDefinition(r"\theta", "learned model parameters", "model dependent"),
    "A": SymbolDefinition(r"A", "first deformation gradient/matrix in the Hadamard jump test", "dimensionless"),
    "B": SymbolDefinition(r"B", "second deformation gradient/matrix in the Hadamard jump test", "dimensionless"),
    "m_jump": SymbolDefinition(r"m", "interface normal/covector in the Hadamard rank-one jump A−B=a⊗m", "reciprocal-direction scale"),
    "chi": SymbolDefinition(r"\chi_1,\chi_2,\chi_3", "ordered eigenvalues of the pairwise twin-compatibility spectral matrix χ", "dimensionless"),
    "Ui": SymbolDefinition(r"U_i", "stretch tensor of martensite variant i", "dimensionless"),
    "Uj": SymbolDefinition(r"U_j", "stretch tensor of martensite variant j", "dimensionless"),
    "sigma2": SymbolDefinition(r"\sigma_2", "middle singular value of a deformation gradient", "dimensionless"),
    "alpha": SymbolDefinition(r"\alpha", "normalized principal stretch/lattice parameter used by the specialized cubic→orthorhombic triplet condition", "dimensionless"),
    "beta_trip": SymbolDefinition(r"\beta", "second normalized principal stretch/lattice parameter in the specialized triplet condition", "dimensionless", "not the monoclinic angle β in the B19′ model"),
    "gamma": SymbolDefinition(r"\gamma", "third normalized principal stretch/lattice parameter in the specialized triplet condition", "dimensionless"),
    "gP": SymbolDefinition(r"g_P", "parent crystal orientation matrix", "rotation matrix"),
    "gD": SymbolDefinition(r"g_D", "daughter crystal orientation matrix", "rotation matrix"),
    "OR": SymbolDefinition(r"OR", "orientation relationship mapping parent and daughter crystallographic frames", "rotation/orientation relation"),
    "variant": SymbolDefinition(r"variant", "symmetry-equivalent crystallographic transformation choice", "discrete index/class"),
    "SPk": SymbolDefinition(r"S_P^{(k)}", "k-th proper parent crystal-symmetry rotation used to generate a daughter orientation branch", "rotation matrix"),
    "Rcp": SymbolDefinition(r"R_{cp}", "proper daughter-to-parent crystal-frame orientation-relationship rotation", "rotation matrix"),
    "BA": SymbolDefinition(r"B_A", "Cartesian direct-lattice basis matrix of the parent/austenite B2 crystal", "length matrix", "columns are the parent direct-lattice basis vectors in an orthonormal Cartesian crystal frame"),
    "BM": SymbolDefinition(r"B_M", "Cartesian direct-lattice basis matrix of the daughter/martensite B19′ crystal", "length matrix", "contains the measured monoclinic lattice lengths and angle"),
    "FMA": SymbolDefinition(r"F_{M\leftarrow A}", "metric-aware correspondence deformation mapping parent Cartesian lattice vectors into daughter Cartesian lattice vectors", "dimensionless deformation gradient", "software bridge used to obtain a model-derived starting OR from the CT/Otsuka–Ren correspondence"),
    "RAM": SymbolDefinition(r"R_{A\to M}", "proper rotational factor of the polar decomposition of the correspondence deformation", "rotation matrix", "its transpose is the daughter→parent rotation R_cp consumed by reconstruction"),
    "gDk": SymbolDefinition(r"g_D^{(k)}", "k-th symmetry-distinct regenerated daughter orientation generated from one parent", "rotation matrix"),
    "delta_cycle": SymbolDefinition(r"\delta_i", "round-trip daughter→parent→daughter closure misorientation for measured daughter grain i", "degrees"),
    "dGD": SymbolDefinition(r"d_{G_D}", "minimum disorientation angle after reducing by daughter crystal symmetry group G_D", "degrees"),
    "gDmeas": SymbolDefinition(r"g_{D,i}^{meas}", "measured daughter orientation of grain i", "rotation matrix"),
    "Pi": SymbolDefinition(r"P_i", "reconstructed parent assigned to daughter grain i", "discrete parent ID"),
    "DeltaSwitch": SymbolDefinition(r"\Delta\theta_{j\to k}", "daughter-symmetry-reduced orientation change from regenerated branch j to branch k", "degrees"),
}


# Each registered equation must have an explicit notation list.  The test suite
# checks coverage against the equation registry.
EQUATION_SYMBOLS: dict[str, tuple[str, ...]] = {
    "PTMC-SC1": ("lambda2",),
    "PTMC-SC2": ("a_twin", "U", "cof", "I", "n_twin"),
    "PTMC-SC3": ("tr", "U", "det", "a_twin", "n_twin"),
    "PTMC-V": ("v", "lambda1", "lambda3", "e1", "e3"),
    "PTMC-IPS": ("F", "I", "d", "m"),
    "PTMC-TRACE": ("lambda1", "lambda3", "delta", "tau"),
    "PTMC-DET": ("lambda1", "lambda3", "delta"),
    "SHEAR-PHI": ("phi", "delta", "tauperp"),
    "SHEAR-SHEAR-GEO": ("a_twin", "phi", "d"),
    "LAMINATE-F": ("P", "f", "F1", "F2", "I", "d", "mA", "mB"),
    "HC": ("HC", "GA", "C_A_M", "GM", "C_M_A"),
    "LEFT-COSET": ("Ci", "gi", "HC", "C_A_M"),
    "N-CORR": ("NC", "card", "GA", "HC"),
    "DOUBLE-COSET": ("Ok", "HC", "gi"),
    "CMC": ("CMC", "C_M_A", "MM", "MA"),
    "CMC-Q": ("uA", "CMC"),
    "G-CMC": ("GCMC", "GA", "q", "u", "SCMC"),
    "CMC-DEGEN": ("qi",),
    "SMC": ("SMC", "MA", "C_A_M", "MM"),
    "SMC-D": ("dA", "SMC", "mA"),
    "SHEAR-SHEAR": ("mA", "n_twin", "dA", "a_twin"),
    "EPSILON": ("epsilon", "mA", "n_twin", "dA", "a_twin", "s"),
    "NITI-CMC": ("CMC", "aB2", "aRatio", "bRatio", "cRatio", "beta"),
    "NITI-Q": ("q1", "q23", "bRatio", "K", "Delta", "aRatio", "cRatio", "beta"),
    "CT-C1": ("bRatio", "sqrt2", "aRatio", "cRatio", "beta"),
    "CT-C2a": ("cRatio", "aRatio", "beta", "bRatio", "sqrt2"),
    "CT-C2b-EQ": ("aRatio", "beta", "bRatio", "sqrt2"),
    "CT-C2b-TABLE": ("aRatio", "beta", "bRatio", "sqrt2"),
    "CT-C3": ("cRatio", "aRatio", "beta", "bRatio", "sqrt2"),
    "CT-D1": ("bRatio", "sqrt2", "cRatio", "aRatio", "beta"),
    "CT-D2": ("aRatio", "cRatio", "sqrt2", "beta"),
    "CT-E": ("aRatio", "bRatio", "cRatio", "sqrt2", "beta", "CMC"),
    "APP-C-O4": ("aRatio", "cRatio", "beta"),
    "FTC": ("C_M_A", "T_M_A", "FA"),
    "TWIN-I": ("pM", "C_M_A", "pA", "s", "tr", "Cint", "MM"),
    "TWIN-II": ("aM", "C_M_A", "aA", "s", "tr", "Cint", "MM"),
    "APP-TEMP": ("RofT", "T", "MA", "MM"),
    "APP-MC": ("xk", "mu", "Sigma", "yk"),
    "APP-OPT": ("xopt", "J"),
    "APP-ML": ("yhat", "ftheta", "theta"),
    "IND-HADAMARD": ("A", "B", "a_twin", "m_jump"),
    "IND-TWIN-SPECTRAL": ("chi", "Ui", "Uj"),
    "IND-LAMINATE": ("F", "f", "U", "a_twin", "n_twin", "sigma2"),
    "IND-TRIPLET-I": ("alpha", "beta_trip", "gamma"),
    "IND-TRIPLET-II": ("alpha", "beta_trip", "gamma"),
    "APP-RECON": ("gP", "OR", "variant", "gD"),
    "APP-CT-OR": ("FMA", "BM", "C_M_A", "BA", "RAM", "U", "Rcp"),
    "APP-REGEN": ("gP", "SPk", "Rcp", "gDk"),
    "APP-CYCLE": ("delta_cycle", "dGD", "gDmeas", "Pi", "gDk"),
    "APP-SWITCH": ("DeltaSwitch", "dGD", "gDk"),
}


INPUT_SYMBOLS: dict[str, SymbolDefinition] = {
    "a_b2": SymbolDefinition(r"a_{B2}", "B2 austenite cubic unit-cell edge length", "Å", "Å (ångström) = 10⁻¹⁰ m"),
    "a_b19p": SymbolDefinition(r"a_{B19'}", "B19′ martensite a-axis unit-cell length", "Å", "physical measured lattice parameter"),
    "b_b19p": SymbolDefinition(r"b_{B19'}", "B19′ martensite b-axis unit-cell length", "Å", "physical measured lattice parameter"),
    "c_b19p": SymbolDefinition(r"c_{B19'}", "B19′ martensite c-axis unit-cell length", "Å", "physical measured lattice parameter"),
    "beta_deg": SymbolDefinition(r"\beta", "monoclinic angle between the B19′ a and c axes", "°", "degrees in the UI; converted internally to radians for trigonometric functions"),
    "cc1_tol": SymbolDefinition(r"\mathrm{tol}_{CC1}", "numerical decision tolerance applied to |λ₂−1|", "dimensionless", "not a physical law; must be reported in a paper"),
    "cc2_tol": SymbolDefinition(r"\mathrm{tol}_{CC2}", "numerical decision tolerance applied to the normalized CC2 residual", "dimensionless", "not a physical law; must be reported in a paper"),
}


OUTPUT_SYMBOLS: dict[str, SymbolDefinition] = {
    "lambda2": SYMBOLS["lambda2"],
    "cmc_distance": SymbolDefinition(r"d_{CMC}", "relative numerical distance of the CMC eigenspectrum from exact degeneracy", "dimensionless", "software residual; not identical to the paper-defined family distances"),
    "cc2_residual": SymbolDefinition(r"r_{CC2}", "scale-normalized absolute residual of the full CC2 equality for the best tested domain system", "dimensionless", "0 is exact"),
    "epsilon": SYMBOLS["epsilon"],
    "detU": SymbolDefinition(r"\det U", "determinant of U; transformation volume ratio in the stretch description", "dimensionless", "1 means zero net volume change"),
    "q": SYMBOLS["qi"],
    "habit_plane": SymbolDefinition(r"m_A=(h,k,l)", "compatible habit-plane normal/covector in parent crystallographic coordinates", "reciprocal direction", "components are directional coordinates; overall scaling does not change the plane"),
    "dA": SYMBOLS["dA"],
    "s": SYMBOLS["s"],
}


def definitions_for_equation(key: str) -> tuple[SymbolDefinition, ...]:
    ids = EQUATION_SYMBOLS.get(key, ())
    return tuple(SYMBOLS[i] for i in ids)


def symbol_rows(definitions: Iterable[SymbolDefinition]) -> list[dict[str, str]]:
    return [asdict(d) for d in definitions]


def symbol_text(definition: SymbolDefinition) -> str:
    unit = f" Unit: {definition.unit}." if definition.unit else ""
    note = f" {definition.note}" if definition.note else ""
    return f"{definition.symbol} = {definition.meaning}.{unit}{note}".strip()

OPERATOR_SYMBOLS: dict[str, SymbolDefinition] = {
    "transpose": SymbolDefinition(r"(\cdot)^T", "transpose: swap matrix/vector rows and columns", "operator"),
    "inverse": SymbolDefinition(r"(\cdot)^{-1}", "matrix inverse: the matrix that undoes the original matrix", "operator"),
    "inverse_transpose": SymbolDefinition(r"(\cdot)^{-T}", "inverse transpose: transpose of the matrix inverse", "operator"),
    "dot": SymbolDefinition(r"\cdot", "dot/inner product: combines two vectors into a scalar", "operator"),
    "outer": SymbolDefinition(r"\otimes", "outer/dyadic product: combines two vectors into a rank-one matrix", "operator"),
    "cof": SymbolDefinition(r"\mathrm{cof}(\cdot)", "cofactor-matrix operator", "operator"),
    "trace": SymbolDefinition(r"\mathrm{tr}(\cdot)", "trace: sum of matrix diagonal entries", "operator"),
    "determinant": SymbolDefinition(r"\det(\cdot)", "determinant: scalar matrix volume/scale factor", "operator"),
    "norm": SymbolDefinition(r"|x|\;\text{or}\;\|x\|", "absolute value for a scalar or norm/length for a vector, as stated by the equation", "operator"),
    "sqrt": SymbolDefinition(r"\sqrt{x}", "positive square root of x", "operator"),
    "plusminus": SymbolDefinition(r"\pm", "two mathematical branches: one with plus and one with minus", "operator"),
    "intersection": SymbolDefinition(r"\cap", "set intersection: elements common to both sets", "set operator"),
    "forall": SymbolDefinition(r"\forall", "for every value in the stated set", "logical symbol"),
    "in": SymbolDefinition(r"\in", "is an element/member of", "set relation"),
    "ge": SymbolDefinition(r"\ge", "greater than or equal to", "comparison"),
    "le": SymbolDefinition(r"\le", "less than or equal to", "comparison"),
    "arrow": SymbolDefinition(r"\to", "maps/transforms from the object on the left to the object on the right", "mapping symbol"),
    "proportional": SymbolDefinition(r"\parallel", "parallel to", "geometric relation"),
    "equals": SymbolDefinition(r"=", "equals: the left and right sides have the same value", "relation"),
    "plus": SymbolDefinition(r"+", "addition: add the quantities on each side of the sign", "arithmetic operator"),
    "minus": SymbolDefinition(r"-", "subtraction or a negative sign, according to its position", "arithmetic operator"),
    "division": SymbolDefinition(r"/\;\text{or}\;\frac{a}{b}", "division: numerator divided by denominator", "arithmetic operator"),
    "power": SymbolDefinition(r"x^p", "power/exponent: multiply x by itself according to exponent p; negative exponents mean reciprocal powers", "operator"),
    "subscript": SymbolDefinition(r"x_i", "subscript: label/index identifying a component, phase, branch, or ordered value; it is not multiplication", "notation"),
    "juxtaposition": SymbolDefinition(r"AB", "adjacent compatible scalars/matrices/vectors mean multiplication or linear action in the order written", "multiplication notation"),
    "sin": SymbolDefinition(r"\sin", "sine trigonometric function of the stated angle", "operator"),
    "cos": SymbolDefinition(r"\cos", "cosine trigonometric function of the stated angle", "operator"),
    "tan": SymbolDefinition(r"\tan", "tangent trigonometric function of the stated angle", "operator"),
    "rank": SymbolDefinition(r"\mathrm{rank}(A)", "matrix rank: number of linearly independent directions carried by matrix A", "operator"),
    "eig": SymbolDefinition(r"\mathrm{eig}(A)", "eigenvalue operation: returns principal scalar factors associated with matrix A", "operator"),
    "iff": SymbolDefinition(r"\Longleftrightarrow", "if and only if: both statements imply each other", "logical relation"),
    "distributed": SymbolDefinition(r"\sim", "is sampled/distributed according to the probability distribution on the right", "statistical relation"),
    "normal_dist": SymbolDefinition(r"\mathcal{N}(\mu,\Sigma)", "multivariate normal (Gaussian) distribution with mean vector μ and covariance matrix Σ", "probability distribution"),
    "minimize": SymbolDefinition(r"\min", "minimize: search for variable values that make the following objective as small as possible", "optimization operator"),
    "degrees": SymbolDefinition(r"^\circ", "degree symbol: the angle is measured in degrees", "angle unit"),
    "matrix": SymbolDefinition(r"\begin{pmatrix}\cdots\end{pmatrix}", "matrix: rectangular array; rows and columns refer to the axes stated beside the output", "notation"),
    "set_builder": SymbolDefinition(r"\{x: condition\}", "set-builder notation: the set of x values satisfying the condition after the colon", "set notation"),
    "closed_interval": SymbolDefinition(r"[a,b]", "closed interval: every value from a through b, including both endpoints", "set notation"),
}


def operators_for_latex(latex: str) -> tuple[SymbolDefinition, ...]:
    """Return plain-language definitions for non-trivial operators present in a relation."""
    checks = [
        ("inverse_transpose", ("^{-T}", "^-T")),
        ("transpose", ("^T", "^{T}")),
        ("inverse", ("^{-1}",)),
        ("dot", ("\\cdot",)),
        ("outer", ("\\otimes",)),
        ("cof", ("\\mathrm{cof}", "\\operatorname{cof}")),
        ("trace", ("\\mathrm{tr}", "\\operatorname{tr}")),
        ("determinant", ("\\det",)),
        ("sqrt", ("\\sqrt",)),
        ("plusminus", ("\\pm", "\\mp")),
        ("intersection", ("\\cap",)),
        ("forall", ("\\forall",)),
        ("in", ("\\in",)),
        ("ge", ("\\ge", "\\geq")),
        ("le", ("\\le", "\\leq")),
        ("arrow", ("\\to", "\\rightarrow")),
        ("proportional", ("\\parallel",)),
        ("equals", ("=",)),
        ("division", ("\\frac", "/")),
        ("sin", ("\\sin",)),
        ("cos", ("\\cos",)),
        ("tan", ("\\tan",)),
        ("rank", ("\\mathrm{rank}", "\\operatorname{rank}")),
        ("eig", ("\\mathrm{eig}", "\\operatorname{eig}")),
        ("iff", ("\\Longleftrightarrow", "\\Leftrightarrow")),
        ("normal_dist", ("\\mathcal{N}",)),
        ("distributed", ("\\sim",)),
        ("minimize", ("\\min",)),
        ("degrees", ("^\\circ", "^{\\circ}")),
        ("matrix", ("\\begin{pmatrix}", "\\begin{bmatrix}")),
        ("set_builder", ("\\{", ":")),
    ]
    out: list[SymbolDefinition] = []
    seen: set[str] = set()
    for key, needles in checks:
        if any(n in latex for n in needles) and key not in seen:
            out.append(OPERATOR_SYMBOLS[key])
            seen.add(key)
    # Generic arithmetic/notation. These are included because the UI is designed
    # for readers who should not need prior mathematical notation knowledge.
    if "+" in latex and "plus" not in seen:
        out.append(OPERATOR_SYMBOLS["plus"]); seen.add("plus")
    # A hyphen can occur in \text, so require common mathematical minus contexts.
    if any(token in latex for token in ("-", "^{-", "-1", "-2", "-3")) and "minus" not in seen:
        out.append(OPERATOR_SYMBOLS["minus"]); seen.add("minus")
    if "^" in latex and "power" not in seen:
        out.append(OPERATOR_SYMBOLS["power"]); seen.add("power")
    if "_" in latex and "subscript" not in seen:
        out.append(OPERATOR_SYMBOLS["subscript"]); seen.add("subscript")
    if "[0,1]" in latex and "closed_interval" not in seen:
        out.append(OPERATOR_SYMBOLS["closed_interval"]); seen.add("closed_interval")
    # Norm/absolute-value notation is intentionally broad; include when vertical bars occur.
    if "|" in latex and "norm" not in seen:
        out.append(OPERATOR_SYMBOLS["norm"]); seen.add("norm")
    return tuple(out)
