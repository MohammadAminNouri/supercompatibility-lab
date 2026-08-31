#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
echo "[1/5] Patching workspace 9 in app.py..."
python tools_upgrade_reconstruction.py
echo "[2/5] Rebuilding all self-contained Streamlit entrypoints..."
python scripts/build_standalone.py
echo "[3/5] Compiling Python sources..."
python -m compileall -q src app.py streamlit_app.py supercompatibility_r7.py supercompatibility_final.py
echo "[4/5] Running reconstruction tests..."
python -m pytest -q tests/test_reconstruction.py tests/test_reconstruction_workbench.py
echo "[5/5] Running deployment preflight..."
python scripts/deployment_preflight.py
printf '\nSUCCESS: reconstruction workbench upgraded.\n'
printf 'Next: git add . && git commit -m "Upgrade parent daughter reconstruction workbench" && git push\n'
