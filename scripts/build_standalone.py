from __future__ import annotations

"""Build deterministic single-file Streamlit entrypoints from src/ + data/.

The generated entrypoint embeds the complete local research engine and required
CSV assets.  This is a deployment safeguard: Streamlit can run the fresh build
even if a browser upload accidentally omits the src/ directory.
"""

import base64
from pathlib import Path
import pprint
import re
import zlib

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "app.py"
OUTPUTS = [ROOT / "app.py", ROOT / "streamlit_app.py", ROOT / "supercompatibility_r7.py", ROOT / "supercompatibility_final.py"]
BUILD_MARKER = "2026-08-31-final-truth-valid-reconstruction"


def enc(raw: bytes) -> str:
    return base64.b64encode(zlib.compress(raw, level=9)).decode("ascii")


def module_payloads() -> dict[str, str]:
    out: dict[str, str] = {}
    for path in sorted((ROOT / "src").glob("*.py")):
        name = "src" if path.name == "__init__.py" else f"src.{path.stem}"
        out[name] = enc(path.read_bytes())
    return out


def data_payloads() -> dict[str, str]:
    return {p.name: enc(p.read_bytes()) for p in sorted((ROOT / "data").iterdir()) if p.is_file()}


def main() -> None:
    source = TEMPLATE.read_text(encoding="utf-8")
    modules = pprint.pformat(module_payloads(), width=10_000, compact=True, sort_dicts=True)
    data = pprint.pformat(data_payloads(), width=10_000, compact=True, sort_dicts=True)

    # Replace the complete embedded-payload region by sentinel markers.  Do not
    # use a line-only regex here: pprint may emit a multi-line dictionary.
    mod_start = source.find("_EMBEDDED_MODULES =")
    root_start = source.find("_EMBED_ROOT =", mod_start)
    if mod_start < 0 or root_start < 0:
        raise SystemExit("Standalone replacement failure: embedded payload markers not found")
    source = (
        source[:mod_start]
        + f"_EMBEDDED_MODULES = {modules}\n\n"
        + f"_EMBEDDED_DATA = {data}\n\n"
        + source[root_start:]
    )
    source, n3 = re.subn(r"(?m)^# Build: .*?$", f"# Build: {BUILD_MARKER}", source, count=1)
    source, n4 = re.subn(
        r'(?m)^_EMBED_ROOT = Path\(_tempfile\.gettempdir\(\)\) / ".*?"$',
        '_EMBED_ROOT = Path(_tempfile.gettempdir()) / "supercompatibility_lab_embedded_final_truth_valid_reconstruction"',
        source,
        count=1,
    )
    if (n3, n4) != (1, 1):
        raise SystemExit(f"Standalone replacement failure: build={n3}, root={n4}")
    for out in OUTPUTS:
        out.write_text(source, encoding="utf-8")
        print(f"wrote {out.name} ({out.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
