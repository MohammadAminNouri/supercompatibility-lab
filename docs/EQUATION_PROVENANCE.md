# Equation-level provenance

This release treats equation provenance as part of the numerical result, not as UI decoration. Each registered calculation has a stable key, displayed relation, source or extension class, source location, implementation path, scope, and caveat. The Streamlit **Math used** expanders and the paper-ready exports are generated from the same registry in `src/provenance.py`.

The primary source for the metric/correspondence and supercompatibility workflow is **Acta Materialia 316 (2026) 122399**, DOI `10.1016/j.actamat.2026.122399`.

## Core chain

The built-in B2 → B19′ workflow records the following chain explicitly:

1. physical lattice inputs and units → phase metrics `M_A`, `M_M`;
2. correspondence → `CMC = C^T M_M C - M_A`;
3. `u^T CMC u = 0` → preserved-length directions;
4. CMC eigenspectrum + degeneracy rule `q_i=0`, `q_j q_k<=0` → A/M habit-plane candidates;
5. `SMC = M_A^-1 - C M_M^-1 C^T`, then `d_A=SMC m_A`;
6. symmetry/correspondence double cosets → transformation-twin candidates;
7. `2(m_A^T n)d_A=a` and epsilon → metric A/M–M/M intercompatibility;
8. independent PTMC stretch construction → SC1/SC2/SC3 cross-check.

The analytical B2/B19′ formulas, C1/C2a/C2b/C3 families, higher-order D1/D2/E degeneracies, and Appendix-C closed-form family are separate registered equations and are cross-checked against the general matrix engine where possible.

## Provenance classes

- **source equation**: directly implements a relation in the primary source;
- **source-derived implementation**: algebraic/numerical implementation derived from a source construction;
- **source discrepancy**: retained conflicting source statements; never silently reconciled;
- **independent theorem implementation**: cross-check from a distinct compatibility theorem;
- **software extension**: temperature sweeps, uncertainty propagation, optimization, ML, reconstruction, or other workflow logic that is not presented as a source theorem.

A paper should cite the underlying scientific source for the equation and identify this software/build for the numerical implementation. The software output is not itself a replacement for a scientific citation.
