# Reproducibility contract

The **Paper-ready audit & export** workspace creates a machine-readable record containing the physical inputs and units, normalized ratios, correspondence matrices, numerical tolerances, software build, Python/NumPy versions, core matrices, principal stretches, CMC eigenspectrum and habit planes, analytical parity checks, compatibility-family residuals, source discrepancies, equation provenance, and claim boundaries.

Every JSON record is fingerprinted with SHA-256. Changing an input, result, provenance field, equation record, tolerance, or caveat changes the hash. The hash is a data-integrity fingerprint; it is not a cryptographic authorship claim.

## Minimum information to report in a paper

For any pass/fail or near-target claim, report:

- the original lattice parameters and units;
- temperature/state/measurement context of those parameters when applicable;
- phase/correspondence convention;
- equation or theorem used;
- numerical residual or inequality margin;
- tolerance used for a numerical label;
- software build ID and random seed for stochastic workflows;
- source discrepancy/caveat if the calculation touches one;
- whether the output is a source equation, independent cross-check, or software extension.

A green UI state by itself is not a publishable result.

## Claim boundary

The engine can justify the crystallographic calculations it performs. It cannot by itself establish fatigue life, hysteresis, chemical realizability of a target lattice, or experimental reversibility. Those claims require independent data and appropriate validation.


## Notation is part of the record

The exported record contains a notation dictionary for physical inputs, normalized ratios, core outputs and each registered equation. Symbol definitions include units/scales and disambiguation notes for overloaded letters such as `a`, `beta` and `T`. Equation provenance additionally records the meaning of non-trivial operators such as transpose, inverse, cofactor, trace, determinant, norm, outer product, set relations, trigonometric functions and inequalities. This prevents a numerical result from becoming uninterpretable when copied from the UI into a notebook or manuscript.
