from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_reconstruction_templates_and_release_script_exist():
    required = [
        ROOT / "data" / "reconstruction_orientations_template.csv",
        ROOT / "data" / "reconstruction_adjacency_template.csv",
        ROOT / "data" / "reconstruction_demo_orientations.csv",
        ROOT / "data" / "reconstruction_demo_adjacency.csv",
        ROOT / "scripts" / "verify_release.py",
    ]
    assert all(p.is_file() and p.stat().st_size > 0 for p in required)


def test_repository_omits_disallowed_author_reference():
    needle = ("ca" + "yron").lower()
    skip = {".git", ".pytest_cache", "__pycache__", ".venv"}
    allowed_suffixes = {".py", ".md", ".txt", ".csv", ".toml", ".yml", ".yaml"}
    hits = []
    for path in ROOT.rglob("*"):
        if not path.is_file() or any(part in skip for part in path.parts) or path.suffix.lower() not in allowed_suffixes:
            continue
        text = path.read_text(encoding="utf-8", errors="ignore").lower()
        if needle in text:
            hits.append(str(path.relative_to(ROOT)))
    assert not hits, hits
