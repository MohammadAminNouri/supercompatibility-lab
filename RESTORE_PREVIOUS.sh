#!/usr/bin/env bash
set -euo pipefail
ROOT="/workspaces/supercompatibility-lab"
cd "$ROOT"

python - <<'PY'
from pathlib import Path
import json, shutil

root=Path("/workspaces/supercompatibility-lab")
pointer=root/".last_v4_2_backup.json"
if not pointer.exists():
    raise SystemExit("No V4.2 backup pointer found.")
backup=Path(json.loads(pointer.read_text())["backup_dir"])
entry=backup/"entrypoints"
names=["app.py","streamlit_app.py","supercompatibility_r7.py","supercompatibility_final.py"]
for name in names:
    src=entry/name
    if not src.exists():
        raise SystemExit(f"Missing backup: {src}")
for name in names:
    shutil.copy2(entry/name,root/name)
    print("RESTORED:",name)

# Restore every policy-cleaned stale file as well.
policy=backup/"policy_hits"
if policy.exists():
    for src in policy.rglob("*"):
        if src.is_file() and src.name!="policy_hits_manifest.json":
            rel=src.relative_to(policy)
            dst=root/rel
            dst.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(src,dst)
            print("RESTORED POLICY FILE:",rel)
print("Restoration complete from:",backup)
PY
