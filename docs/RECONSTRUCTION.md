# Parent ↔ daughter reconstruction workbench

This document describes the final academic reconstruction workflow implemented in `src/reconstruction.py`, `src/reconstruction_workbench.py`, `src/reconstruction_reporting.py` and `src/reconstruction_academic.py`.

## What the module solves

The workspace separates two different crystallographic problems.

1. **Daughter → parent reconstruction:** measured daughter orientations, crystallographic symmetry, an orientation relationship (OR), and daughter-grain adjacency are used to reconstruct prior-parent clusters and parent orientations.
2. **Parent → daughter variant generation:** a known parent orientation and OR are used to enumerate all symmetry-distinct daughter orientations. This predicts allowed orientations only; it does not predict which variant nucleates, volume fraction, morphology or spatial position.

## Minimum inverse-reconstruction input

The reconstruction engines require:

- one orientation per daughter grain;
- unique daughter `grain_id` values;
- parent point-group symmetry;
- daughter point-group symmetry;
- a parent/daughter orientation relationship;
- a daughter-grain neighbor graph.

A daughter orientation can be supplied as Bunge Euler angles `phi1_deg, Phi_deg, phi2_deg` or a scalar-first unit quaternion `qw,qx,qy,qz`. The internal matrix convention is **crystal → specimen**. A specimen → crystal input is accepted only when explicitly selected and is transposed before analysis.

## Raw EBSD import

The final workbench accepts:

- EDAX/TSL-style `.ang` point maps;
- Oxford/HKL-style `.ctf` point maps;
- delimited `.csv`, `.tsv` and `.txt` point tables;
- pre-segmented grain-level delimited tables.

The lightweight importer is intentionally **fail-closed**: it refuses to continue when required columns cannot be identified instead of silently guessing missing crystallography.

### ANG handling

The standard ANG v5 layout is interpreted as:

`phi1 PHI phi2 x y IQ CI phase SEM fit`

Standard ANG Euler angles are treated as radians and converted to degrees. If angle magnitudes are clearly inconsistent with a radian file, the importer marks that degrees were inferred and displays a warning. Quality columns such as IQ, CI, SEM signal and fit are preserved when present.

### CTF handling

The point-data header must contain `Phase`, `X`, `Y`, `Euler1`, `Euler2` and `Euler3`. Euler angles are interpreted as degrees. Common quality fields including Bands, Error, MAD, BC and BS are preserved when present.

### Reference-frame guard

Vendor text formats may encode specimen axes, scan axes and Euler reference frames differently. The app therefore requires the user to confirm that the exported Euler/map reference frame has been checked before raw-point segmentation is enabled. This is a deliberate research-quality guard: a wrong reference frame can produce a plausible-looking but physically wrong reconstruction.

## Raw-point phase and quality filtering

Before segmentation the workbench shows:

- parsed row count;
- indexed row count and indexed fraction;
- available phase IDs and point counts;
- inferred point spacing;
- available vendor quality fields;
- the detected Euler-unit convention.

Phase 0 is treated as non-indexed. The daughter phase must be selected explicitly. Optional CI, MAD or ANG-fit filters are **off by default** and are never applied silently. If enabled, the threshold is carried into the methods record.

## Daughter-grain segmentation

Raw EBSD pixels are converted to daughter grains only after explicit user choices.

Pixels are connected when:

1. they are spatial nearest neighbors;
2. they belong to the same indexed phase;
3. their symmetry-reduced disorientation is below the chosen grain-boundary threshold.

The exposed segmentation controls are:

- grain-boundary disorientation threshold: UI range 0.5–15°;
- minimum indexed points per retained daughter grain: 2–1000;
- spatial-neighbor radius divided by the median point spacing: 1.05–1.70.

These ranges are interface safety ranges, not universal material constants. The default 3° grain threshold is only a starting value and must be justified for the dataset.

The grain mean orientation is computed after aligning crystal-symmetry-equivalent point orientations. The output also records the RMS and maximum intragranular orientation spread. Approximate area and shared-boundary contact measures are retained from the point grid.

For very large production maps, a dedicated EBSD package is still preferred for segmentation. The in-app raw-point segmenter is intentionally capped and designed for transparent, reproducible analyses and cross-checks.

## Adjacency

Measured/shared-boundary adjacency is preferred. A table requires:

`grain_id_1, grain_id_2`

Optional boundary metadata may be retained. Manual adjacency is supported for small examples. If only grain centroids exist, an explicitly labelled centroid k-nearest-neighbor graph can be generated, but it is marked approximate and must not be described as measured boundary topology.

## Orientation relationship

The OR is stored internally as a proper daughter→parent crystal-frame rotation `R_cp`.

The user can choose:

- a built-in KS/NW/Bain/Pitsch/Burgers family;
- custom parent/daughter plane and direction parallelisms;
- a custom 3×3 OR rotation matrix.

The exact matrix used by the calculation is displayed and exported.

### OR refinement

A supplied physical OR can be refined within a bounded rotation neighborhood. The objective uses neighboring daughter-grain parent-candidate consistency with a robust loss. The search is deterministic and bounded.

The refinement reports:

- objective before;
- objective after;
- OR correction angle;
- refined OR matrix.

An improved numerical objective is not proof that the refined OR is physically correct. Independent retained-parent measurements, pole parallelisms or another reconstruction package should be used when available.

## Five inverse-reconstruction routes

### 1. Neighbor voting

Each daughter grain considers all parent candidates generated from its orientation and the OR. Neighboring grains vote for locally compatible candidates. Selected candidates are merged spatially when their parent disorientation is below the merge tolerance.

Use: fast, transparent local-consensus baseline.

### 2. Grain graph + Markov clustering

Daughter grains are graph nodes. Edge weights decrease with the best common-parent candidate mismatch. Markov clustering partitions the graph into prior-parent regions.

Use: global partitioning of a connected daughter network.

Important control: graph inflation, which affects cluster granularity.

### 3. Variant graph candidate propagation

Each daughter grain keeps several possible parent candidates instead of selecting one immediately. Candidate compatibility is propagated through a graph before a parent is selected.

Use: ambiguous transformations where preserving alternate candidates helps resolve prior-parent boundaries.

The implementation is a transparent research member of the variant-graph family, not a binary reproduction of external MTEX code.

### 4. Nucleation + growth

Strict parent seeds are created first; neighboring daughter grains are then added using a looser growth tolerance.

Use: microstructures with reliable local nuclei and a physically meaningful seed/growth hierarchy.

The nucleation and growth tolerances are reported separately.

### 5. Operator / groupoid consistency

The OR and parent symmetry generate theoretical daughter–daughter transformation operators. Every measured neighboring daughter pair is matched to its nearest theoretical operator modulo daughter symmetry. Operator-consistent edges are then subjected to common-parent candidate consistency before parent regions are grown.

This route is the closest workbench analogue to the classic ARPGE/GenOVa operator-groupoid reconstruction idea. It additionally exports the operator residual for every daughter-neighbor edge and operator-frequency statistics.

External method reference: J. Appl. Cryst. 40 (2007) 1183–1188, DOI `10.1107/S0021889807048777`.

## Daughter-level output

Every measured daughter grain can be exported with:

- `parent_id` — reconstructed prior-parent cluster;
- `candidate_variant_id` — selected symmetry-generated candidate index;
- `best_OR_fit_deg` — best candidate-to-parent angular residual;
- `second_best_candidate_fit_deg` — next-best candidate residual;
- `candidate_separation_deg` — second-best minus best residual;
- `candidate_count` — number of crystallographically distinct parent candidates considered;
- `absolute_fit_support_0_to_1`;
- `candidate_separation_support_0_to_1`;
- `support_score_0_to_1`;
- reconstructed parent Euler angles and quaternion;
- original position, area, phase and grain-spread metadata when available.

### Heuristic support score

The exact score is deliberately exposed:

\[
s_{abs}=\exp\left[-\frac12\left(\frac{\theta_{best}}{3^\circ}\right)^2\right]
\]

\[
s_{sep}=1-\exp\left[-\frac{\theta_{2nd}-\theta_{best}}{3^\circ}\right]
\]

\[
s=0.55s_{abs}+0.45s_{sep}.
\]

This is a deterministic ambiguity/support diagnostic, **not a calibrated probability**.

## Parent-level output

Each reconstructed parent row contains:

- parent ID;
- number of supporting daughter grains;
- distinct selected variants;
- dominant selected variant and fraction;
- supporting area and area fraction when available;
- mean, median, P95 and maximum OR residual;
- mean candidate separation;
- mean and minimum support;
- weak/ambiguous daughter count;
- parent centroid when coordinates exist;
- reconstructed parent Bunge Euler angles;
- reconstructed parent quaternion;
- a descriptive quality flag.

## Variant and operator statistics

Every reconstruction route exports selected candidate-variant frequency statistics. The operator/groupoid route additionally exports:

- neighboring daughter grain IDs;
- nearest theoretical operator ID;
- operator residual in degrees;
- whether the residual passes the operator threshold;
- whether the reconstructed result places the pair in the same parent;
- accepted operator frequencies and residual statistics.

Candidate/variant numbers are internal enumeration labels unless a transformation-specific canonical packet/Bain/variant mapping is explicitly supplied. They must not be silently renamed as canonical labels.

## Cross-method comparison

The final workspace deliberately avoids an opaque “best reconstruction score”. Instead it reports several independent questions.

### Method comparison table

For each method:

- reconstructed parent count;
- singleton-parent fraction;
- mean, median and P95 OR residual;
- fraction with OR residual >5°;
- mean support and low-support fraction;
- number/fraction of adjacency edges called prior-parent boundaries;
- runtime;
- known-truth clustering accuracy and ARI when truth exists.

### Adjusted Rand Index (ARI)

ARI compares the parent partition from two methods independent of arbitrary parent IDs. `ARI = 1` means identical partitions; values near zero indicate chance-level agreement.

### Normalized Mutual Information (NMI)

NMI compares clustering information. `NMI = 1` means identical clustering information. NMI is not chance-adjusted, so it is complementary to ARI.

### Boundary Jaccard agreement

Each method classifies the same daughter adjacency edges as “same parent” or “prior-parent boundary”. The Jaccard score compares the sets of boundary edges. `1` means identical boundary calls.

### Matched parent-orientation disagreement

Prior-parent clusters from two methods are first matched by maximum daughter-grain overlap using a Hungarian assignment. The parent orientations are then compared modulo parent symmetry. The resulting matrix is reported in degrees; smaller is stronger orientation agreement.

### Boundary consensus table

Every daughter adjacency edge receives the fraction of selected methods that call it a prior-parent boundary. Unanimous boundary or unanimous same-parent edges are the strongest cross-method evidence; disagreement edges should be inspected on the map.

## Independent retained-parent validation

Independently measured parent orientations can be uploaded after reconstruction. The nearest parent-symmetry-reduced misorientation is reported. Without a known spatial/ID correspondence, this is only a nearest-reference diagnostic, not proof that two IDs refer to the same physical grain.

## Academic evidence ZIP

The workspace can generate a single reconstruction evidence ZIP containing:

- analyzed daughter-grain table;
- adjacency table;
- exact OR matrix and metadata;
- orientation convention;
- selected methods and all controls;
- method-comparison CSV;
- ARI and NMI matrices;
- boundary Jaccard matrix;
- matched parent-orientation disagreement matrix;
- boundary consensus table;
- parent summary for every method;
- daughter evidence table for every method;
- variant frequency tables;
- operator edge/frequency tables when the operator route is used;
- per-method diagnostics JSON;
- reconstructed parent orientations;
- a human-readable Methods/interpretation README.

This bundle is intended to make a manuscript result traceable without relying on the Streamlit session state.

## Relation to established tools

The workbench follows the same broad research workflow used by established packages: verify EBSD reference frames, segment measured daughter grains, define/refine an OR, reconstruct prior parents, evaluate fit, and transfer parent/variant information back to the measured microstructure.

Relevant external method families include:

- ARPGE/GenOVa operator-groupoid reconstruction — DOI `10.1107/S0021889807048777`;
- MTEX parent-grain reconstruction — DOI `10.1107/S1600576721011560`;
- variant graph parent reconstruction — DOI `10.1016/j.mtla.2022.101399`.

The local implementations are transparent and testable but are **not** claimed to be binary-identical reproductions of MTEX, ARPGE, OIM Analysis, AZtecCrystal or vendor software.

## Minimum defensible Methods report

At minimum record:

1. daughter phase and symmetry;
2. parent phase and symmetry;
3. exact OR definition and OR matrix;
4. whether/how the OR was refined;
5. EBSD file/source software and verified reference frame;
6. daughter-grain segmentation threshold and minimum point count;
7. daughter adjacency source;
8. reconstruction method(s);
9. every method threshold/graph control;
10. daughter-grain count and reconstructed-parent count;
11. OR residual distribution, not only the mean;
12. support/ambiguity definition;
13. cross-method cluster/boundary/orientation agreement when multiple methods are used;
14. retained-parent or known-truth validation when available;
15. exact variant/operator ID convention before interpreting selection statistics.


## B19′ → B2 → B19′ round-trip and re-transformation

After a daughter→parent reconstruction, the workbench can reuse the exact reconstructed B2 parent orientations and the exact daughter→parent OR matrix to regenerate every symmetry-distinct B19′ orientation branch. No OR is silently re-entered or transposed between the inverse and forward steps.

For parent orientation \(g_P\), parent proper-symmetry operation \(S_P^{(k)}\), and daughter→parent crystal-frame OR rotation \(R_{cp}\), the regenerated daughter branch is

\[
g_D^{(k)}=g_P S_P^{(k)}R_{cp}.
\]

For measured daughter grain \(i\), already assigned to reconstructed parent \(P_i\), round-trip closure is

\[
\delta_i=\min_k d_{G_D}\left(g_{D,i}^{\mathrm{meas}},g_{D,P_i}^{(k)}\right),
\]

where \(d_{G_D}\) is the minimum disorientation modulo daughter crystal symmetry. The best regenerated branch is the corresponding \(\operatorname*{arg\,min}\). The software reports the best residual, second-best residual, branch separation and a user-thresholded interpretation band. These thresholds are reporting choices, not material constants.

The pairwise regenerated-branch orientation-change catalogue uses

\[
\Delta\theta_{j\rightarrow k}=d_{G_D}\left(g_D^{(j)},g_D^{(k)}\right).
\]

This is purely orientation geometry. It does **not** predict nucleation probability, transformed volume fraction, morphology, stress selection or energy barriers.

### Metric-aware NiTi natural/AQ quick setup

For NiTi the dedicated quick setup implements the reported natural/AQ orientation relation

\[
(010)_{B19'}\parallel(110)_{B2},\qquad
[101]_{B19'}\parallel[\bar 1 1 1]_{B2}.
\]

The implementation does not treat monoclinic indices as Euclidean vectors. Direct-lattice directions are converted using \(\mathbf d=B\,u\), plane normals using the reciprocal basis \(\mathbf n\propto B^{-T}h\), and the two orthonormal crystal triads are then related by a proper rotation \(R_{cp}\). The exact B2/B19′ lattice parameters used to form the OR are stored with the calculation. With cubic B2 and monoclinic B19′ proper symmetries this construction yields 12 symmetry-distinct regenerated B19′ orientation branches per B2 parent.

Literature basis for the natural/AQ relation: *What EBSD and TKD Tell Us about the Crystallography of the Martensitic B2-B19′ Transformation in NiTi Shape Memory Alloys*, Crystals 10 (2020) 562, DOI `10.3390/cryst10070562`.

### Independent later-cycle EBSD

A separately measured later-cycle B19′ map can be uploaded as an already segmented grain CSV or raw ANG/CTF text EBSD. Raw maps use the same explicit phase selection, reference-frame confirmation and grain-segmentation controls used elsewhere in the reconstruction workbench. If later-cycle grains carry a known `parent_id`, each is tested only against that reconstructed B2 parent. Otherwise all reconstructed parents are compared and the best-versus-second-parent angular separation is reported so orientation-only parent-assignment ambiguity is visible.

This later-cycle comparison is stronger evidence than the round-trip closure because the new B19′ orientations were not used to reconstruct the original B2 parent.

### Academic cycle export

The cycle evidence ZIP contains:

- exact regenerated daughter orientation library (Euler, quaternion and 3×3 matrices);
- measured round-trip closure per daughter grain;
- parent-level closure and branch-coverage summary;
- observed regenerated-branch occupancy;
- complete branch-to-branch orientation-change catalogue;
- exact OR matrix, symmetry labels and interpretation thresholds;
- independently measured later-cycle matches when provided;
- a scope README explicitly distinguishing internal consistency from independent validation and allowed orientation branches from variant-selection probability.
