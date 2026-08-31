from __future__ import annotations
import ast,base64,zlib,tempfile,sys,hashlib
from pathlib import Path
import numpy as np, pandas as pd

APP=Path('supercompatibility_final.py')
text=APP.read_text(encoding='utf-8')
compile(text,str(APP),'exec')
print('PASS: Python syntax')

tree=ast.parse(text)
mods=None
for n in tree.body:
    if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='_EMBEDDED_MODULES' for t in n.targets):
        mods=ast.literal_eval(n.value); break
assert mods
root=Path(tempfile.mkdtemp(prefix='scv2_')); (root/'src').mkdir()
for name,payload in mods.items():
    s=zlib.decompress(base64.b64decode(payload)).decode('utf-8')
    p=root/('src/__init__.py' if name=='src' else f"src/{name.rsplit('.',1)[-1]}.py")
    p.write_text(s,encoding='utf-8')
sys.path.insert(0,str(root))
from src.core import LatticeInput, normalized_metrics, cmc_matrix, cmc_degeneracy, smc_matrix
from src.symmetry import cubic_point_group, double_cosets
from src.equation_engine import correspondence_left_cosets
from src.reconstruction import bunge_euler_to_matrix
from src.presets import PRESETS
from src.provenance import BUILD_ID

ns={
    'np':np,'pd':pd,'LatticeInput':LatticeInput,'cubic_point_group':cubic_point_group,
    'double_cosets':double_cosets,'correspondence_left_cosets':correspondence_left_cosets,
    'bunge_euler_to_matrix':bunge_euler_to_matrix,'BUILD_ID':BUILD_ID,
    'IW_SOURCE':{},'IW_REORIENTATION_SOURCE':{},'json':__import__('json'),
}
want={
'niti_reference_distortion_from_lattice','symmetry_generated_distortions','uniaxial_stress_tensor',
'interaction_work_mj_m3','reorientation_iw_table','stress_b2_from_ebsd_parent',
'_matrix_member','correspondence_operator_map','operator_composition_table'
}
for n in tree.body:
    if isinstance(n,ast.FunctionDef) and n.name in want:
        exec(compile(ast.Module(body=[n],type_ignores=[]),str(APP),'exec'),ns)

inp=PRESETS['Published binary NiTi example']
F0=ns['niti_reference_distortion_from_lattice'](inp)
states=ns['symmetry_generated_distortions'](F0)
assert len(states)==12
assert np.allclose(states[0],F0)
print('PASS: 12 unique distortion states, reference state first')

sigma=ns['uniaxial_stress_tensor'](50,[1,0,0])
df,states2,F02=ns['reorientation_iw_table'](inp,sigma,0)
assert len(df)==12 and np.isfinite(df['IW (MJ/m³)']).all()
selfrow=df[df['Same state']].iloc[0]
assert abs(float(selfrow['IW (MJ/m³)']))<1e-12
print('PASS: Xiao-style F_i->j reorientation table; self-transition IW = 0')

# Compression reverses every IW because sigma changes sign.
dfc,_,_=ns['reorientation_iw_table'](inp,-sigma,0)
a=df.set_index('Final state')['IW (MJ/m³)'].sort_index().to_numpy()
b=dfc.set_index('Final state')['IW (MJ/m³)'].sort_index().to_numpy()
assert np.allclose(a,-b)
print('PASS: tension/compression sign behavior')

# Identity Euler orientation with specimen X equals direct [100] B2 loading.
sigma_e,n=ns['stress_b2_from_ebsd_parent'](50,0,0,0,[1,0,0])
assert np.allclose(n,[1,0,0]) and np.allclose(sigma_e,sigma)
print('PASS: EBSD orientation loading transform for identity orientation')

opmap=ns['correspondence_operator_map']()
assert opmap.shape==(12,12)
valid={d.label for d in double_cosets()}
assert set(opmap.to_numpy().ravel()).issubset(valid)
print('PASS: 12x12 correspondence operator map')

assert len(cubic_point_group())==48 and len(correspondence_left_cosets())==12 and len(double_cosets())==7
print('PASS: 48 / 12 / 7 CT group-theory structure')

required=[
'1 · Enter the NiTi lattice you want to study',
'2 · What do you want to calculate?',
'Initial martensite distortion state',
'I have a reconstructed B2 EBSD orientation',
'F_{i\\to j}=F_jF_i^{-1}',
]
for s in required: assert s in text,s
print('PASS: guided input and interpretation UI hooks present')

assert 'if not (workspace.startswith("3") or workspace.startswith("9")):' in text
print('PASS: Workspace 3 uses local inputs instead of duplicated top-of-page inputs')

print('SHA256:',hashlib.sha256(APP.read_bytes()).hexdigest())
print('ALL GUIDED WORKSPACE CHECKS PASS')
