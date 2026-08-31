# Final release manifest

**Build ID:** `2026-08-31-final-notation-explicit`

**Recommended Streamlit entrypoint:** `supercompatibility_final.py`

## Scientific release checks completed before packaging

- 63/63 automated tests passed.
- Equation/source parity audit passed.
- Deployment preflight passed.
- Embedded-engine scientific self-test passed.
- Numerical release audit passed.
- Deterministic embedded-bundle audit passed: all 23 `src` Python modules and all 9 data assets match the self-contained Streamlit bundle byte-for-byte.
- `app.py`, `streamlit_app.py`, `supercompatibility_r7.py`, and `supercompatibility_final.py` are byte-identical deployment aliases.
- Parent/daughter reconstruction validation recovered the deterministic two-parent benchmark with all five implemented methods.
- The source-discrepancy register is preserved in `docs/SOURCE_DISCREPANCIES.md`; discrepancies are surfaced rather than silently reconciled.
- Formula notation, operators, primary inputs, and core outputs are covered by the explicit notation/provenance layer and tested in `tests/test_symbol_completeness.py`.

## GitHub Actions

The files `.github/workflows/ci.yml` and `.github/workflows/tests.yml` trigger automatically on pushes and pull requests. The main CI installs the pinned dependencies and performs a real Streamlit health check.

## Streamlit Community Cloud

Use:

- Branch: `main`
- Main file path: `supercompatibility_final.py`
- Python: `3.12`

The visible build identifier must be `2026-08-31-final-notation-explicit`.

## Environment note

The packaging container used for this release did not have Streamlit installed, so the Streamlit server itself was not launched locally during packaging. The final GitHub Actions workflow installs the pinned Streamlit dependency and requires the `/_stcore/health` endpoint to pass before CI is green.
