from __future__ import annotations

import ast
import base64
import hashlib
from pathlib import Path
import sys
import zlib

ROOT = Path(__file__).resolve().parents[1]
APP = ROOT / "app.py"


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def assignments() -> dict[str, dict[str, str]]:
    tree = ast.parse(APP.read_text(encoding="utf-8"))
    out: dict[str, dict[str, str]] = {}
    for node in tree.body:
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        target = node.targets[0]
        if isinstance(target, ast.Name) and target.id in {"_EMBEDDED_MODULES", "_EMBEDDED_DATA"}:
            out[target.id] = ast.literal_eval(node.value)
    return out


def main() -> int:
    values = assignments()
    modules = values.get("_EMBEDDED_MODULES")
    data = values.get("_EMBEDDED_DATA")
    if not modules or data is None:
        print("FAIL: app.py is not the self-contained build.")
        return 2

    errors: list[str] = []
    for name, payload in modules.items():
        decoded = zlib.decompress(base64.b64decode(payload))
        path = ROOT / "src" / ("__init__.py" if name == "src" else name.split(".")[-1] + ".py")
        if not path.is_file():
            errors.append(f"missing source file: {path.relative_to(ROOT)}")
            continue
        actual = path.read_bytes()
        if decoded != actual:
            errors.append(
                f"embedded/source mismatch: {path.relative_to(ROOT)} "
                f"embedded={sha(decoded)[:12]} source={sha(actual)[:12]}"
            )

    for name, payload in data.items():
        decoded = zlib.decompress(base64.b64decode(payload))
        path = ROOT / "data" / name
        if not path.is_file():
            errors.append(f"missing data file: {path.relative_to(ROOT)}")
            continue
        actual = path.read_bytes()
        if decoded != actual:
            errors.append(
                f"embedded/data mismatch: {path.relative_to(ROOT)} "
                f"embedded={sha(decoded)[:12]} data={sha(actual)[:12]}"
            )

    if errors:
        print("FAIL: deterministic Streamlit bundle audit")
        print("\n".join(f"- {e}" for e in errors))
        return 1

    print(f"PASS: {len(modules)} embedded Python modules match src/ byte-for-byte")
    print(f"PASS: {len(data)} embedded data files match data/ byte-for-byte")
    return 0


if __name__ == "__main__":
    sys.exit(main())
