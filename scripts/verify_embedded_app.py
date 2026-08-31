from __future__ import annotations

import ast
import base64
import hashlib
from pathlib import Path
import sys
import zlib

ROOT = Path(__file__).resolve().parents[1]
APP = ROOT / "app.py"
ENTRYPOINTS = [ROOT / "app.py", ROOT / "streamlit_app.py", ROOT / "supercompatibility_r7.py", ROOT / "supercompatibility_final.py"]


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


def expected_module_names() -> set[str]:
    names = set()
    for path in (ROOT / "src").glob("*.py"):
        names.add("src" if path.name == "__init__.py" else f"src.{path.stem}")
    return names


def main() -> int:
    values = assignments()
    modules = values.get("_EMBEDDED_MODULES")
    data = values.get("_EMBEDDED_DATA")
    if not modules or data is None:
        print("FAIL: app.py is not the self-contained build.")
        return 2

    errors: list[str] = []
    expected_modules = expected_module_names()
    if set(modules) != expected_modules:
        errors.append(f"embedded module set mismatch: missing={sorted(expected_modules-set(modules))}, extra={sorted(set(modules)-expected_modules)}")
    expected_data = {p.name for p in (ROOT / "data").iterdir() if p.is_file()}
    if set(data) != expected_data:
        errors.append(f"embedded data set mismatch: missing={sorted(expected_data-set(data))}, extra={sorted(set(data)-expected_data)}")

    for name, payload in modules.items():
        decoded = zlib.decompress(base64.b64decode(payload))
        path = ROOT / "src" / ("__init__.py" if name == "src" else name.split(".")[-1] + ".py")
        if not path.is_file():
            errors.append(f"missing source file: {path.relative_to(ROOT)}")
            continue
        actual = path.read_bytes()
        if decoded != actual:
            errors.append(f"embedded/source mismatch: {path.relative_to(ROOT)} embedded={sha(decoded)[:12]} source={sha(actual)[:12]}")

    for name, payload in data.items():
        decoded = zlib.decompress(base64.b64decode(payload))
        path = ROOT / "data" / name
        if not path.is_file():
            errors.append(f"missing data file: {path.relative_to(ROOT)}")
            continue
        actual = path.read_bytes()
        if decoded != actual:
            errors.append(f"embedded/data mismatch: {path.relative_to(ROOT)} embedded={sha(decoded)[:12]} data={sha(actual)[:12]}")

    if all(p.is_file() for p in ENTRYPOINTS):
        digests = {p.name: sha(p.read_bytes()) for p in ENTRYPOINTS}
        if len(set(digests.values())) != 1:
            errors.append("self-contained entrypoints are not byte-identical: " + str({k: v[:12] for k, v in digests.items()}))
    else:
        errors.append("one or more self-contained entrypoints are missing")

    if errors:
        print("FAIL: deterministic Streamlit bundle audit")
        print("\n".join(f"- {e}" for e in errors))
        return 1

    print(f"PASS: {len(modules)} embedded Python modules match src/ byte-for-byte and cover every src module")
    print(f"PASS: {len(data)} embedded data files match data/ byte-for-byte and cover every data asset")
    print("PASS: app.py, streamlit_app.py, supercompatibility_r7.py and supercompatibility_final.py are byte-identical")
    return 0


if __name__ == "__main__":
    sys.exit(main())
