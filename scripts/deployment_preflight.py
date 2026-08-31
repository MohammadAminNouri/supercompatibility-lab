from __future__ import annotations

import importlib
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

REQUIRED = [
    "app.py",
    "streamlit_app.py",
    "supercompatibility_r7.py",
    "supercompatibility_final.py",
    "requirements.txt",
    "src/__init__.py",
    "src/core.py",
    "src/ptmc.py",
    "src/equation_engine.py",
    "src/provenance.py",
    "src/reproducibility.py",
    "src/symbols.py",
    "src/reconstruction.py",
    "scripts/verify_equation_parity.py",
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
    "src.equation_engine",
    "src.provenance",
    "src.reproducibility",
    "src.symbols",
):
    importlib.import_module(module)

build = (ROOT / "BUILD_ID.txt").read_text(encoding="utf-8").strip()
print(f"DEPLOYMENT PRECHECK PASS · build={build}")
