# Data directory

- `literature.csv` — curated evidence/reference table; missing values are intentionally blank.
- `temperature_sweep_template.csv` — required schema for temperature-dependent lattice evaluation.
- `ml_template.csv` — example schema for user-trained composition/processing→lattice regression.
- `multistep_template.csv` — general cell/correspondence schema for stage-wise transformation analysis.

See `docs/DATA_SCHEMA.md` for units and interpretation.

## Parent reconstruction files

- `reconstruction_orientations_template.csv` — grain-level daughter orientation input template.
- `reconstruction_adjacency_template.csv` — neighboring daughter-grain pairs.
- `reconstruction_demo_orientations.csv` / `reconstruction_demo_adjacency.csv` — deterministic synthetic validation example with known parent labels.

Prefer measured/segmented grain-boundary adjacency for publication work. Centroid k-nearest-neighbor adjacency in the app is an explicitly approximate fallback.
