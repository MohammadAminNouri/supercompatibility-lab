from pathlib import Path
import ast,base64,hashlib,sys,tempfile,zlib
import numpy as np
p=Path("supercompatibility_final.py")
t=p.read_text(encoding="utf-8")
compile(t,str(p),"exec")
print("PASS: Python syntax")
tr=ast.parse(t)
mods=None
for n in tr.body:
    if isinstance(n,ast.Assign) and any(isinstance(x,ast.Name) and x.id=="_EMBEDDED_MODULES" for x in n.targets):
        mods=ast.literal_eval(n.value)
        break
assert mods and len(mods)==27, len(mods or {})
canon="".join(k+"\0"+mods[k]+"\0" for k in sorted(mods))
assert hashlib.sha256(canon.encode()).hexdigest()=="f6a9bc226d14da563547dfe0754124bd9c62cd5eb7daa790f6cae7581124157d"
print("PASS: 27 embedded scientific modules preserved exactly")
w3=t[t.index("# ---------- 3 variants / operators / twins / IW ----------"):t.index("# ---------- 4 temperature + uncertainty ----------")]
assert hashlib.sha256(w3.encode()).hexdigest()=="3ce2b6dd0007df961475ea3143a5ca7a801a7ad25b02035971db436da32fc8ef"
print("PASS: Workspace 3 preserved exactly")
# Extract embedded engine and run numerical smoke/regression checks.
root=Path(tempfile.mkdtemp(prefix="supercompat_v4_verify_")); (root/"src").mkdir()
for name,payload in mods.items():
    s=zlib.decompress(base64.b64decode(payload)).decode("utf-8")
    f=root/"src/__init__.py" if name=="src" else root/("src/"+name.rsplit(".",1)[-1]+".py")
    f.write_text(s,encoding="utf-8")
sys.path.insert(0,str(root))
from src.core import LatticeInput, normalized_metrics, cmc_matrix, cmc_degeneracy
from src.distances import compatibility_dashboard, all_cofactor_systems
from src.ptmc import stretch_from_lattice
from src.symmetry import double_cosets
from src.equation_engine import correspondence_left_cosets, verify_analytic_cmc_against_general
inp=LatticeInput(3.01,2.898,4.108,4.646,97.78)
ma,mm=normalized_metrics(inp); deg=cmc_degeneracy(cmc_matrix(ma,mm)); d=compatibility_dashboard(inp); st=stretch_from_lattice(inp)
assert len(correspondence_left_cosets())==12
assert len(double_cosets())==7
assert not deg.exact
assert abs(deg.relative_zero-0.27120708)<5e-6, deg.relative_zero
assert abs(st.eigenvalues[1]-d.lambda2)<1e-12
chk=verify_analytic_cmc_against_general(inp)
assert chk["matrix_frobenius_residual"]<1e-10
assert chk["eigenvalue_max_abs_residual"]<1e-10
print("PASS: default NiTi CMC/PTMC/analytical parity regression")
for phrase in ["Compatibility verdict & diagnostics","Classical PTMC / cofactor verification","Temperature & measurement uncertainty","Experimental context & literature","Inverse lattice design","EBSD parent ↔ daughter reconstruction","Independent theorem cross-checks","Applicability, references & methods","Equation & analytical solution library","Manuscript audit & reproducible export","Download COMPLETE manuscript evidence bundle (ZIP)"]:
    assert phrase in t, phrase
print("PASS: final academic UX/workspace hooks present")
assert "use_container_width=" not in t
print("PASS: deprecated Streamlit width API absent")
print("SHA256:",hashlib.sha256(p.read_bytes()).hexdigest())
print("ALL FINAL-ACADEMIC V4 CHECKS PASS")
