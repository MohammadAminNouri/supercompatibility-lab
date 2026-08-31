#!/usr/bin/env bash
set -euo pipefail
printf '\n[1/6] Checking repo root...\n'
test -f supercompatibility_final.py || { echo 'Run this from the repository root.'; exit 1; }
printf '[2/6] Installing dev dependencies...\n'
python -m pip install -r requirements-dev.txt
printf '[3/6] Checking required folders...\n'
test -f src/core.py
test -f data/literature.csv
test -f scripts/verify_release.py
test -f tests/test_core.py
printf '[4/6] Pytest collection check...\n'
python -m pytest --collect-only -q
printf '[5/6] Scientific smoke checks...\n'
python research_selftest.py || true
python scripts/verify_equation_parity.py
python scripts/verify_release.py
printf '[6/6] Running complete test suite...\n'
python -m pytest -q

echo
echo 'REPAIR COMPLETE. Now commit:'
echo '  git add src data scripts tests .github requirements*.txt pyproject.toml'
echo '  git commit -m "Restore research package and CI folders"'
echo '  git push'
