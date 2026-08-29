# Data schemas and units

## 1. Temperature sweep

Template: `data/temperature_sweep_template.csv`

Required columns:

| Column | Unit | Meaning |
|---|---:|---|
| `temperature_K` | K | measurement/evaluation temperature |
| `a_B2_A` | Å | B2 parent cubic lattice parameter |
| `a_B19p_A` | Å | B19′ a-axis |
| `b_B19p_A` | Å | B19′ b-axis |
| `c_B19p_A` | Å | B19′ c-axis |
| `beta_deg` | ° | B19′ monoclinic angle |

The app recomputes compatibility for every row. It does not assume a thermal-expansion law.

## 2. Composition→lattice ML training data

Template: `data/ml_template.csv`

Mandatory targets:

| Column | Unit |
|---|---:|
| `a_B2_A` | Å |
| `a_B19p_A` | Å |
| `b_B19p_A` | Å |
| `c_B19p_A` | Å |
| `beta_deg` | ° |

All other numeric columns can be selected as features. Suggested examples:

- atomic % of alloying elements
- aging temperature [K or °C, but use one convention consistently]
- aging time [h]
- solution treatment temperature
- cold-work fraction [%]
- grain size [nm]

Do not mix units within a feature column.

Candidate composition grids need the same feature columns used to train the model; lattice target columns are not required for prediction.

## 3. Multi-step transformation chain

Template: `data/multistep_template.csv`

Each row describes one parent→product stage.

Cell parameters:

- `parent_a`, `parent_b`, `parent_c`
- `parent_alpha_deg`, `parent_beta_deg`, `parent_gamma_deg`
- `product_a`, `product_b`, `product_c`
- `product_alpha_deg`, `product_beta_deg`, `product_gamma_deg`

Lengths must use one consistent unit within a stage (Å recommended). Angles are degrees.

Correspondence matrix columns:

```text
C11,C12,C13,
C21,C22,C23,
C31,C32,C33
```

The matrix maps the parent crystallographic direction coordinates into the product crystallographic direction coordinates under the convention used by the generic stage engine.

## 4. Literature database

File: `data/literature.csv`

This is a curated evidence table, not a homogeneous meta-analysis dataset. Records intentionally contain missing values where the source used for curation did not report a quantity.

`record_type` distinguishes experimental, mechanism, method, device and theory records.

## 5. Missing values

For research data use empty cells/NaN for unknown quantities. Do **not** encode an unknown physical value as zero. Some UI fields use zero as an interactive “not supplied” convenience, but exported research CSVs should use proper missing values.

## Parent/daughter reconstruction orientation table

Required column:

- `grain_id` — unique integer grain identifier, dimensionless.

Supply **one** orientation representation. The app also asks explicitly whether the exported orientation maps crystal→specimen or specimen→crystal; the internal convention is crystal→specimen and specimen→crystal inputs are transposed on import. Do not guess this setting—use the convention of the originating EBSD software/export.


### Bunge Euler representation

- `phi1_deg` — Bunge \(\phi_1\), degrees;
- `Phi_deg` — Bunge \(\Phi\), degrees;
- `phi2_deg` — Bunge \(\phi_2\), degrees.

### Quaternion representation

- `qw`, `qx`, `qy`, `qz` — unit quaternion components, dimensionless, scalar-first convention.

Optional columns:

- `x`, `y` — grain-centroid coordinates in any consistent spatial unit, used only for plotting or the explicitly approximate k-nearest-neighbor fallback;
- `true_parent_id` — optional ground-truth label for synthetic/benchmark data only.

Template: `data/reconstruction_orientations_template.csv`.

## Parent/daughter reconstruction adjacency table

Required columns:

- `grain_id_1`
- `grain_id_2`

Each row means the two daughter grains share a measured/segmented neighborhood relation. Duplicate and self edges are removed internally. Template: `data/reconstruction_adjacency_template.csv`.

For academic reconstruction work, true grain-boundary adjacency from the segmented EBSD map is preferred over centroid k-NN adjacency.
