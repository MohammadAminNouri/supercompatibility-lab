# Upload the final build to GitHub

Upload the **contents of this folder** to the repository root. Do not upload the ZIP as a single file and do not create an extra enclosing directory.

## Required final Streamlit entry point

`supercompatibility_final.py`

This filename is deliberately new so a stale Streamlit deployment cannot accidentally reuse an older `app.py` revision.

## Required root files

- `supercompatibility_final.py`
- `app.py`
- `streamlit_app.py`
- `requirements.txt`
- `BUILD_ID.txt`
- `STREAMLIT_ENTRYPOINT.txt`
- `src/`
- `data/`
- `docs/`
- `tests/`
- `scripts/`
- `.github/workflows/ci.yml`
- `.github/workflows/tests.yml`

## Automatic GitHub Actions

Uploading/committing these files automatically triggers both workflows. The main workflow installs dependencies, verifies equation/source parity, runs the complete scientific tests, checks the embedded bundle byte-for-byte, starts `supercompatibility_final.py`, and requires the Streamlit health endpoint to respond successfully.

## Streamlit Community Cloud

Create a **new deployment** with:

- Repository: your GitHub repository
- Branch: `main`
- Main file path: `supercompatibility_final.py`
- Python: 3.12

The visible build identifier must be `2026-08-31-final-ct-otsuka-reconstruction-ui`.
