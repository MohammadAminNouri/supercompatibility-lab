# Final release manifest

**Build ID:** `2026-08-31-final-truth-valid-reconstruction`

**Recommended Streamlit entrypoint:** `supercompatibility_final.py`

## Scientific release checks completed before packaging

- 96/96 collected automated tests passed in validated batches after the CT/Otsuka–Ren and UI additions.
- Equation/source parity audit passed.
- Deployment preflight passed.
- Embedded-engine scientific self-test passed.
- Numerical release audit passed.
- Deterministic embedded-bundle audit passed: all 27 `src` Python modules and all 9 data assets match the self-contained Streamlit bundle byte-for-byte.
- `app.py`, `streamlit_app.py`, `supercompatibility_r7.py`, and `supercompatibility_final.py` are byte-identical deployment aliases.
- Obsolete upgrade/bootstrap backups and stale R6 deployment instructions were removed from the GitHub-ready release root.
- Parent/daughter reconstruction validation recovered the deterministic two-parent benchmark with all five implemented methods.
- Reconstruction-specific tests cover robust ANG/CTF parsing, explicit phase/reference-frame guards, raw-point segmentation, candidate ambiguity, ARI/NMI/boundary/orientation agreement, operator statistics and the academic evidence ZIP.
- Cycle tests cover the metric-aware NiTi B2↔B19′ natural/AQ OR, the CT/Otsuka–Ren model-derived polar starting OR, 12 regenerated B19′ branches for both NiTi routes, daughter→parent→daughter closure, branch occupancy, branch-switch geometry, independent later-cycle daughter matching and the cycle evidence ZIP.
- The dedicated CT/Otsuka–Ren option exposes the exact correspondence matrix, metric-aware deformation, polar rotation, principal stretches and the explicit caveat that correspondence is not itself a unique experimental OR.
- The cycle workflow explicitly distinguishes **internal round-trip consistency** from **independent later-cycle EBSD validation** and never turns an allowed daughter orientation into a claimed nucleation probability.
- The source-discrepancy register is preserved in `docs/SOURCE_DISCREPANCIES.md`; discrepancies are surfaced rather than silently reconciled.
- Formula notation, operators, primary inputs, and core outputs are covered by the explicit notation/provenance layer and tested in `tests/test_symbol_completeness.py`.

## GitHub Actions

The files `.github/workflows/ci.yml` and `.github/workflows/tests.yml` trigger automatically on pushes and pull requests. The main CI installs the pinned dependencies and performs a real Streamlit health check.

## Streamlit Community Cloud

Use:

- Branch: `main`
- Main file path: `supercompatibility_final.py`
- Python: `3.12`

The visible build identifier must be `2026-08-31-final-truth-valid-reconstruction`.

## Environment note

The packaging container used for this release does not provide the Streamlit browser runtime. GitHub Actions installs the pinned Streamlit dependency and requires the `/_stcore/health` endpoint to pass before CI is green.


## Reconstruction truth-validation hardening

- Built-in synthetic daughter data are generated from the **currently selected OR and parent/daughter symmetries**; no hidden KS/FCC→BCC dataset is reused for NiTi CT/Otsuka–Ren validation.
- Known-truth validation uses **Adjusted Rand Index (ARI), homogeneity, completeness, V-measure, and parent-boundary precision/recall/F1**. The legacy majority-remapped clustering accuracy is not used in the academic workbench because singleton over-segmentation can make it falsely equal 100%.
- Cross-method ARI/NMI are explicitly labelled as **agreement, not accuracy**.
- Detailed method evidence reports singleton-parent fraction and warns when low OR residuals are trivial singleton self-fits.
