# Manuscript methods template

Use this as a checklist, not as text to paste blindly.

## Crystallographic inputs

Report the parent and product phase, lattice parameters with units and measurement temperature/state, monoclinic angle, and correspondence matrix/convention.

## A/M compatibility

State whether the calculation used the PTMC middle-stretch condition, the CMC degeneracy condition, or both. Report the principal stretches, `|lambda2-1|`, CMC eigenvalues, degeneracy order, habit-plane solution(s), and numerical tolerance. If a closed-form B2/B19′ family is invoked, identify C1/C2a/C2b/C3 or D1/D2/E and report its equality residual and inequality margins.

## M/M compatibility

State how correspondence variants/intercorrespondence classes were generated and which parent symmetry element generated the twin. Report the twin type, shear plane/direction, amplitude, and the relevant group/coset counts if they are part of the result.

## A/M–M/M intercompatibility

Report the SMC-derived `d_A`, selected habit-plane normalization, twin `a,n`, the shear/shear relation, and epsilon. Do not replace epsilon with a qualitative adjective without the numerical value.

## Independent cross-checks

If used, identify the independent method separately: Hadamard rank-one jump, pairwise twin spectral criterion, all-volume-fraction numerical scan, or transformation-class-specific triplet condition. A pass in one framework must not be reported as a pass in another framework unless the implication has actually been established.

## Computational reproducibility

Report the build ID, software version, tolerances, stochastic seed, equation/provenance export, and the SHA-256 fingerprint of the exact result record. Attach the JSON/CSV audit files as supplementary data when practical.

## Source discrepancies

If the result touches the q3 branch, C2b inequality, or the later distance-table equality, state which interpretation was used and why. Never hide the alternative printed source statement.
