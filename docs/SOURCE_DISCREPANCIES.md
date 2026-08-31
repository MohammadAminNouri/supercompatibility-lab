# Source discrepancy register

Reproducibility requires preserving source inconsistencies rather than choosing one silently. This build records three discrepancies from the primary 2026 source.

## 1. Closed-form CMC q3 branch

The available PDF rendering/extraction is ambiguous and can appear to print the same minus branch for both `q2` and `q3`. The printed CMC matrix has two analytical branches. The implementation uses

`q2=(K-sqrt(Delta))/4`, `q3=(K+sqrt(Delta))/4`

because those values reproduce the numerical eigenspectrum of the printed CMC matrix to floating-point precision. The app reports the ambiguity every time the analytical eigenspectrum is used.

## 2. C2b inequality direction

The analytical first-order C2b equation/derivation and the plotted C2b example use `a=1`, `beta=90°`, `b<=sqrt(2)`. A later distance table prints `b>=sqrt(2)`. The implementation exposes both interpretations. The analytical equation is the default compatibility branch; the table branch is stored as a **source discrepancy**, not treated as a second hidden truth.

## 3. Table distance expression versus printed numerical value

For the shared C2a/C3 equality distance, the later table's printed numerical value `0.269706` is reproduced by the linear lattice-ratio distance `|c-c_target|`. The visually/parsed squared expression produces a different value for the same raw benchmark lattice constants. Both values are exported with explicit labels.

## Policy

No discrepancy is resolved solely by convenience. A resolution must be supported by algebraic parity, an unambiguous corrected source, or an explicit research decision reported in the manuscript. The app therefore keeps analytical parity residuals and source-discrepancy text in the reproducibility record.
