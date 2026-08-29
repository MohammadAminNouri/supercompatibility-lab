# R6 deployment — use a fresh app

This package deliberately includes a brand-new Streamlit entrypoint:

`supercompatibility_r6.py`

Use that exact file when creating the Streamlit app. Do not point the new deployment at an old `app.py` app instance.

Recommended clean deployment:

1. Create a fresh GitHub repository, e.g. `supercompatibility-lab-r6`.
2. Upload the **contents** of this ZIP to the repository root.
3. Confirm `.github/workflows/ci.yml` appears in GitHub.
4. Wait for `R6 Scientific + Streamlit CI` to turn green.
5. In Streamlit Community Cloud create a **new app** with:
   - Repository: the fresh repository
   - Branch: `main`
   - Main file path: `supercompatibility_r6.py`
   - Python: 3.12
6. The page must display: `Build 2026-08-29-r6 · FRESH REPO`.

Why a new file name? A stale Streamlit deployment that keeps executing an old `app.py` cannot accidentally execute `supercompatibility_r6.py`, because that path did not exist in the old revision.
