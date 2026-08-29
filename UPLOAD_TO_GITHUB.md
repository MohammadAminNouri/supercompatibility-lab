# Upload / deployment checklist

This release is intentionally self-contained. `app.py` embeds the research engine as a deployment fallback, while the normal `src/` tree remains in the repository for tests, review and reproducibility.

## Required repository root

After extraction, upload the **contents of this folder**, not the ZIP file itself and not an extra enclosing folder.

The GitHub repository root must visibly contain:

- `app.py`
- `requirements.txt`
- `BUILD_ID.txt`
- `src/`
- `data/`
- `.streamlit/`
- `.github/workflows/ci.yml`
- `.github/workflows/tests.yml`

The expected build identifier is:

`2026-08-29-r4-deployment-rebuild`

The app displays this build identifier directly under the title. If it does not appear, Streamlit is deploying an older commit/branch/repository.

## GitHub Actions

Both workflows trigger automatically on every push. `ci.yml` runs the scientific suite and a real Streamlit health check. `tests.yml` is a fast repository-layout sanity check and intentionally replaces the older workflow filename too.

If your operating system hides `.github`, a visible backup copy is provided as `GITHUB_ACTIONS_ci.yml`. The workflow only works from `.github/workflows/ci.yml`, so ensure the hidden directory was actually uploaded.

## Streamlit Cloud

Configure:

- Repository: the repository containing this release
- Branch: the branch where these files were committed
- Main file path: `app.py`
- Python: 3.12 (the repository includes `.python-version`)

Then reboot the app. A successful deployment shows the `r4` build identifier below the app title.
