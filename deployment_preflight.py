from __future__ import annotations

import importlib
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

REQUIRED = [
    "app.py",
    "requirements.txt",
    "src/__init__.py",
    "src/core.py",
    "src/ptmc.py",
    "src/reconstruction.py",
    ".streamlit/config.toml",
    ".github/workflows/ci.yml",
]

missing = [p for p in REQUIRED if not (ROOT / p).is_file()]
if missing:
    raise SystemExit("DEPLOYMENT PRECHECK FAILED: missing: " + ", ".join(missing))

for package in ("numpy", "pandas", "scipy", "plotly", "sklearn"):
    importlib.import_module(package)

for module in (
    "src.core",
    "src.ptmc",
    "src.symmetry",
    "src.compatibility_methods",
    "src.reconstruction",
    "src.design",
):
    importlib.import_module(module)

build = (ROOT / "BUILD_ID.txt").read_text(encoding="utf-8").strip()
print(f"DEPLOYMENT PRECHECK PASS · build={build}")
