from __future__ import annotations
import ast
import base64
import hashlib
import sys
import tempfile
import zlib
from pathlib import Path
import numpy as np

APP = Path("supercompatibility_final.py")
if not APP.exists():
    raise SystemExit("FAIL: supercompatibility_final.py not found in current folder.")

text = APP.read_text(encoding="utf-8")
compile(text, str(APP), "exec")
print("PASS: Python syntax")

tree = ast.parse(text)
embedded = None
for node in tree.body:
    if isinstance(node, ast.Assign):
        if any(isinstance(t, ast.Name) and t.id == "_EMBEDDED_MODULES" for t in node.targets):
            embedded = ast.literal_eval(node.value)
            break
if not embedded:
    raise SystemExit("FAIL: embedded research engine not found.")

root = Path(tempfile.mkdtemp(prefix="supercompat_verify_"))
(root/"src").mkdir()
for name, payload in embedded.items():
    src = zlib.decompress(base64.b64decode(payload)).decode("utf-8")
    path = root/"src/__init__.py" if name == "src" else root/f"src/{name.rsplit('.',1)[-1]}.py"
    path.write_text(src, encoding="utf-8")

sys.path.insert(0, str(root))
from src.symmetry import cubic_point_group, double_cosets
from src.equation_engine import correspondence_left_cosets
from src.presets import PRESETS

G = cubic_point_group()
cosets = correspondence_left_cosets()
dcs = double_cosets()

assert len(G) == 48, len(G)
assert len(cosets) == 12, len(cosets)
assert len(dcs) == 7, len(dcs)
print("PASS: 48 cubic operations / 12 correspondence variants / 7 double-coset classes")

def member(m, mats):
    return any(np.allclose(m, x, atol=1e-10, rtol=0.0) for x in mats)

for ci in cosets:
    gi=np.asarray(ci.representative,float)
    for cj in cosets:
        gj=np.asarray(cj.representative,float)
        rel=np.linalg.inv(gi)@gj
        hits=[d.label for d in dcs if member(rel,d.matrices)]
        assert len(hits)==1, hits
print("PASS: every ordered variant pair maps to exactly one operator class")

inp = PRESETS["Published binary NiTi example"]
a0,am,bm,cm=inp.a_b2,inp.a_b19p,inp.b_b19p,inp.c_b19p
beta=np.deg2rad(inp.beta_deg); s2=np.sqrt(2.0)
b2=np.array([[am/a0,0,0],[0,bm/(s2*a0),-cm/(s2*a0)],[0,bm/(s2*a0),cm/(s2*a0)]])
mono=np.array([[np.sin(beta),0,0],[0,1,0],[np.cos(beta),0,1]])
basis=np.array([[1,0,0],[0,.5,.5],[0,-.5,.5]])
F0=b2@mono@basis
variants=[]
for g in G:
    f=np.linalg.inv(g)@F0@g
    if not any(np.allclose(f,q,atol=1e-9,rtol=0.0) for q in variants):
        variants.append(f)
assert len(variants)==12, len(variants)
print("PASS: reference NiTi distortion generates 12 unique IW variants")

sigma=100.0*np.outer([1.,0.,0.],[1.,0.,0.])
vals=[float(np.sum(sigma*(f-np.eye(3)))) for f in variants]
assert all(np.isfinite(vals))
print("PASS: IW calculation finite for all 12 variants")

required = [
    "Variants, operators, twins & interaction work",
    "interaction_work_variant_table",
    "correspondence_operator_map",
    "operator_composition_table",
    "Download IW evidence record (JSON)",
]
for item in required:
    assert item in text, item
print("PASS: integrated UI/workflow hooks present")

if "use_container_width=" in text:
    raise SystemExit("FAIL: deprecated use_container_width remains.")
print("PASS: deprecated Streamlit width API removed")

print("SHA256:", hashlib.sha256(APP.read_bytes()).hexdigest())
print("ALL UPGRADE CHECKS PASS")
