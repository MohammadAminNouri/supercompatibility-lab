# Contributing

The main rule is: **do not let UI convenience silently change the crystallography.**

## Before a pull request

```bash
python -m pip install -r requirements-dev.txt
pytest -q
python -m compileall -q app.py src tests
ruff check app.py src tests scripts
```

## Calculation changes

If you change a formula, coordinate convention, normalization or correspondence matrix:

1. state the mathematical equation and source/derivation in the PR;
2. add a benchmark test with expected numerical values;
3. explain whether previous saved results change;
4. avoid replacing exact residuals with undocumented scalar scores.

## New transformations

Do not reuse the built-in B2→B19′ twin machinery for a new transformation merely because both phases have similar names or lattice systems. A full transformation module should define:

- parent/product point groups;
- correspondence matrix;
- subgroup/variant construction;
- applicable Type-I/II/compound domain equations;
- numerical benchmark data.

## Literature database

Every new row in `data/literature.csv` should have:

- a DOI or stable source;
- a concise evidence note stating what is actually supported;
- blanks instead of invented values;
- a record type distinguishing experiment, mechanism, theory, method or device.

## ML contributions

Do not commit a pretrained model unless its training data, license, split strategy, validation metrics and applicability domain are documented. A physics-screened prediction is still only as credible as the lattice model producing it.
