# Release audit

Release audit date: 2026-08-29

This file records the checks performed on the packaged research release.

## Automated scientific/software tests

Command:

```bash
python -m pytest -q
```

Result:

```text
43 passed
```

The suite covers the metric/CMC/SMC engine, PTMC/cofactor cross-checks, Type-I/Type-II/Compound domain classification, symmetry and double-coset invariants, independent compatibility methods, temperature evaluation, uncertainty propagation, inverse design, ML screening, multi-step transformations, parent/daughter orientation reconstruction, OR refinement and release structure.

## Independent numerical benchmark audit

Command:

```bash
python scripts/verify_release.py
```

Result:

```text
Binary benchmark
  ratios = (0.9627907, 1.3647841, 1.5435216)
  lambda2 = 0.9650480588
  stretch variants = 12
  rank-one-compatible variant pairs = 42/66
C1 all-volume-fraction audit
  selected domain = Compound 1
  max |sigma2-1| = 8.882e-16
Synthetic parent-reconstruction audit
  Neighbor voting: parents=2, accuracy=1.000, mean_fit=0.174 deg
  Grain graph + Markov clustering: parents=2, accuracy=1.000, mean_fit=0.174 deg
  Variant graph: parents=2, accuracy=1.000, mean_fit=0.174 deg
  Nucleation + growth: parents=2, accuracy=1.000, mean_fit=0.174 deg
  Operator / groupoid consistency: parents=2, accuracy=1.000, mean_fit=0.174 deg
Release audit: PASS
```

## Syntax/compilation audit

Command:

```bash
python -m py_compile app.py src/*.py scripts/*.py
```

Result: PASS.

## Naming/privacy audit

A case-insensitive recursive scan was performed for the excluded author name requested for this repository. No occurrence was found in source, documentation, data or configuration files included in the release.

## Streamlit runtime note

The build container used for this release does not contain the `streamlit` package and cannot be used to launch the full Streamlit server. The application entry point and all Python modules compile successfully, and the numerical/reconstruction engines are covered by the automated test suite. A real browser-level Streamlit rendering check should still be performed after installing `requirements.txt` in the target environment.

## Scientific scope note

Passing these tests establishes internal numerical consistency and reproduction of the encoded benchmark cases. It does not prove that a proposed composition is synthesizable, that a machine-learning model extrapolates reliably outside its training domain, or that crystallographic compatibility alone determines hysteresis, fatigue life or functional performance. Those claims require independent experimental validation.
