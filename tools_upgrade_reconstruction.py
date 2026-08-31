from pathlib import Path

ROOT = Path(__file__).resolve().parent
APP = ROOT / "app.py"
START = '# ---------- 9 parent / daughter reconstruction ----------'
END = '# ---------- 10 independent compatibility methods ----------'

text = APP.read_text(encoding='utf-8')
i = text.find(START)
j = text.find(END, i)
if i < 0 or j < 0:
    raise SystemExit('Could not find workspace 9/10 markers in app.py; no file was changed.')

block = '''# ---------- 9 parent / daughter reconstruction ----------\nelif workspace.startswith("9"):\n    from src.reconstruction_workbench import render_reconstruction_workbench\n    render_reconstruction_workbench()\n    math_used(["APP-RECON"], title="Reconstruction method provenance")\n\n\n'''
backup = APP.with_suffix('.py.pre_reconstruction_workbench_backup')
if not backup.exists():
    backup.write_text(text, encoding='utf-8')
APP.write_text(text[:i] + block + text[j:], encoding='utf-8')
print('Patched app.py workspace 9 and saved backup:', backup.name)
