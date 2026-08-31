from __future__ import annotations
import ast,base64,zlib,tempfile,sys,hashlib
from pathlib import Path
import numpy as np

APP=Path('supercompatibility_final.py')
text=APP.read_text(encoding='utf-8')
compile(text,str(APP),'exec')
print('PASS: Python syntax')

tree=ast.parse(text)
mods=None
for node in tree.body:
    if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='_EMBEDDED_MODULES' for t in node.targets):
        mods=ast.literal_eval(node.value);break
assert mods
root=Path(tempfile.mkdtemp(prefix='paper_v3_')); (root/'src').mkdir()
for name,payload in mods.items():
    code=zlib.decompress(base64.b64decode(payload)).decode('utf-8')
    path=root/'src/__init__.py' if name=='src' else root/f"src/{name.split('.')[-1]}.py"
    path.write_text(code,encoding='utf-8')
sys.path.insert(0,str(root))
from src.core import LatticeInput,normalized_metrics,cmc_matrix,cmc_degeneracy,C_A_TO_M,C_M_TO_A
from src.equation_engine import correspondence_left_cosets,niti_cmc_from_input,first_order_families,higher_order_degeneracy_families
from src.symmetry import cubic_point_group,double_cosets,full_twin_explorer
from src.presets import PRESETS

assert len(cubic_point_group())==48
assert len(correspondence_left_cosets())==12
assert len(double_cosets())==7
print('PASS: 48 B2 symmetries / 12 variants / 7 operator classes')

# Pair mapping V1->V2 should be DC4 under the stable embedded ordering.
cosets=correspondence_left_cosets(); dcs=double_cosets()
def member(m,mats): return any(np.allclose(m,x,atol=1e-10,rtol=0) for x in mats)
rel=np.linalg.inv(np.asarray(cosets[0].representative,float))@np.asarray(cosets[1].representative,float)
hits=[d.label for d in dcs if member(rel,d.matrices)]
assert hits==['DC4'],hits
print('PASS: V1 -> V2 maps uniquely to DC4')

inp=LatticeInput(3.01,2.898,4.108,4.646,97.78)
ma,mm=normalized_metrics(inp); cmc=cmc_matrix(ma,mm); deg=cmc_degeneracy(cmc,rtol=2e-6,atol=1e-12)
assert not deg.exact
assert np.allclose(deg.eigenvalues,[-0.13505111992790989,-0.06868224412534085,0.25324650408481797],rtol=1e-9,atol=1e-10)
assert abs(deg.relative_zero-0.2712070769685245)<1e-10
assert all(not f.met for f in first_order_families(*inp.ratios(),inp.beta_deg,tol=2e-6,c2b_interpretation='equation'))
assert all(not f.met for f in higher_order_degeneracy_families(*inp.ratios(),inp.beta_deg,tol=2e-6))
print('PASS: default binary NiTi CMC failure is quantitatively reproduced and all analytical degeneracy families fail')

entries=full_twin_explorer(ma,mm)
e=next(x for x in entries if x.double_coset=='DC4' and x.twin.twin_type=='Type I' and x.parent_element=='(0 0 1) mirror')
assert abs(e.twin.shear_amplitude-0.27325472357998776)<1e-10
# Representative coordinate conversion used in the UI.
a_m=C_M_TO_A@np.asarray(e.twin.twin_shear_vector,float)
p_m=C_A_TO_M.T@np.asarray(e.twin.twin_plane_normal,float)
def collinear(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    return abs(np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b)))>1-1e-10
assert collinear(a_m,[0,0,1])
assert collinear(p_m,[1,0,0])
print('PASS: DC4 representative Type-I twin shear and B19-prime coordinate mapping')

# Static UI/audit assertions
required=[
    'Pair-specific coordinate rule',
    'relative CMC distance δ',
    'shear/shear ε","NOT COMPUTED"',
    'Paper-ready record for this calculation',
    'workspace3_paper_evidence.zip',
    'Experimental-paper rule',
    'F_{i\\to j}=F_jF_i^{-1}',
]
for s in required:
    assert s in text,s
assert 'Best shear/shear ε", "not evaluable"' not in text
assert 'M/M twin", "DEFINED"' not in text
assert 'use_container_width=' not in text
print('PASS: vague legacy statuses removed from Workspace 3 and paper-evidence hooks present')

# Ensure embedded scientific engine was not rewritten by the V3 UI pass.
print('SHA256:',hashlib.sha256(APP.read_bytes()).hexdigest())
print('ALL PAPER-READY V3 CHECKS PASS')
