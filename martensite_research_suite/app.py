import io, json
from pathlib import Path
import numpy as np, pandas as pd, streamlit as st
import plotly.graph_objects as go

from mrs import __version__, BUILD_ID
from mrs.crystal import cubic_proper_group, monoclinic_2_group, niti_natural_or
from mrs.groupoid import common_subgroup, variants_from_cosets, operator_classes, composition_table, arrow_operator_matrix
from mrs.rotations import bunge_euler, angle_axis
from mrs.iw import small_strain, green_lagrange, uniaxial_stress, interaction_work, variant_strains, max_iw_uniaxial
from mrs.compatibility import lambda2_residual, det_volume_ratio
from mrs.exporting import evidence_zip_bytes

st.set_page_config(page_title="Martensite Research Suite", layout="wide")
REFS=Path("data/references.md").read_text()
EQS=Path("data/equations.md").read_text()

st.title("Martensite Research Suite")
st.caption(f"Clean-room research build · {BUILD_ID} · v{__version__}")
st.info("Purpose: connect crystallographic variants/operators, mechanics (interaction work), compatibility diagnostics, reconstruction evidence and reproducible academic export. Published methods are cited; this GUI/code is an independent implementation.")

with st.sidebar:
    st.header("Study definition")
    system=st.selectbox("Transformation preset",["NiTi B2 → B19′ · natural AQ OR"])
    st.markdown("**Notation guard**  \n`(hkl)` = plane · `[uvw]` = direction · `Vᵢ` = variant · `Oₖ` = operator class.")
    aP=st.number_input("a_B2 (Å) — cubic parent lattice parameter",value=3.010000,format="%.6f")
    aM=st.number_input("a_B19′ (Å)",value=2.898000,format="%.6f")
    bM=st.number_input("b_B19′ (Å)",value=4.108000,format="%.6f")
    cM=st.number_input("c_B19′ (Å)",value=4.646000,format="%.6f")
    beta=st.number_input("β_B19′ (deg) — monoclinic angle between a and c",value=97.780000,format="%.6f")

R0,Bp,Bd=niti_natural_or(aP,aM,bM,cM,beta)
Gp=cubic_proper_group(); Gd=monoclinic_2_group()
H=common_subgroup(Gp,Gd,R0,tol_deg=1e-4)
variants,reps,cosets=variants_from_cosets(Gp,H,R0)
op_sets,index_to_op=operator_classes(Gp,H)
A=arrow_operator_matrix(Gp,reps,index_to_op)
comp=composition_table(Gp,op_sets,index_to_op)

tabs=st.tabs(["Context","Variants","Operators & groupoid","Interaction work","Compatibility","Reconstruction validation","Paper export"])

with tabs[0]:
    st.subheader("What this app is doing")
    st.markdown(r"""
**Crystallography engine:** parent symmetry + daughter symmetry + an orientation relationship define a common subgroup. Cosets give crystallographically distinct variants; double cosets give classes of variant-to-variant relationships.

**Mechanics engine:** a stress tensor and a transformation/deformation strain tensor give interaction work \(IW=\sigma:\varepsilon\). This is kept separate from the groupoid engine.

**Evidence engine:** every table can be exported together with equations, references, input metadata and build ID. Cross-method agreement is not treated as accuracy unless known truth exists.
""")
    st.markdown("### NiTi preset provenance")
    st.markdown(r"""
The preset uses the **natural AQ OR** reported for B2→B19′ NiTi:

\[
(010)_{B19'} \parallel (110)_{B2},\qquad
[101]_{B19'} \parallel [\bar{1}11]_{B2}.
\]

The literature correspondence is not silently replaced by KS/NW/Pitsch/Bain.
""")
    st.metric("Parent proper rotations |Gᴾ|",len(Gp))
    st.metric("Common subgroup |Hᴾ|",len(H))
    st.metric("Generated variants Nᵥ",len(variants))
    if len(variants)!=12:
        st.warning("The current OR/symmetry tolerance did not produce 12 variants. Inspect the subgroup/OR before using the result in a paper.")

with tabs[1]:
    st.subheader("Orientational variants")
    rows=[]
    for i,R in enumerate(variants,1):
        e=bunge_euler(R)
        rows.append({"variant":f"V{i}","phi1_deg":e[0],"Phi_deg":e[1],"phi2_deg":e[2],
                     **{f"R{r+1}{c+1}":R[r,c] for r in range(3) for c in range(3)}})
    vdf=pd.DataFrame(rows)
    st.dataframe(vdf,use_container_width=True)
    st.caption("Bunge Euler angles are a display representation of each rotation matrix; the matrix is the primary stored result.")

    pole=np.array([R@np.array([0.,1.,0.]) for R in variants])
    z=np.where(pole[:,2]<0,-1,1); pole=pole*z[:,None]
    X=pole[:,0]/(1+pole[:,2]); Y=pole[:,1]/(1+pole[:,2])
    fig=go.Figure(go.Scatter(x=X,y=Y,mode="markers+text",text=[f"V{i+1}" for i in range(len(variants))],textposition="top center"))
    fig.update_layout(title="Stereographic projection of daughter [010] directions in parent frame",xaxis_title="X",yaxis_title="Y",yaxis_scaleanchor="x")
    st.plotly_chart(fig,use_container_width=True)

with tabs[2]:
    st.subheader("Operator classes")
    st.write(f"Number of double-coset operator classes: **{len(op_sets)}**")
    op_rows=[]
    for k,S in enumerate(op_sets):
        angles=[]
        for gi in S:
            a,ax=angle_axis(Gp[gi]); angles.append((a,ax))
        a,ax=min(angles,key=lambda x:x[0])
        op_rows.append({"operator":f"O{k}","double_coset_size":len(S),"minimum_parent_rotation_deg":a,
                        "axis_x":ax[0],"axis_y":ax[1],"axis_z":ax[2]})
    odf=pd.DataFrame(op_rows); st.dataframe(odf,use_container_width=True)

    st.markdown("#### Variant → variant operator map")
    adf=pd.DataFrame([[f"O{x}" for x in row] for row in A],
                     index=[f"V{i+1}" for i in range(len(variants))],
                     columns=[f"V{i+1}" for i in range(len(variants))])
    st.dataframe(adf,use_container_width=True)

    st.markdown("#### Reduced multivalued composition table")
    cdf=pd.DataFrame(index=[f"O{i}" for i in range(len(op_sets))],columns=[f"O{i}" for i in range(len(op_sets))])
    for i in range(len(op_sets)):
        for j in range(len(op_sets)):
            cdf.iloc[i,j]="{" + ", ".join(f"O{k}" for k in comp[(i,j)]) + "}"
    st.dataframe(cdf,use_container_width=True)
    st.caption("A cell may contain more than one operator. That is intentional: the reduced operator composition is generally multivalued.")

with tabs[3]:
    st.subheader("Interaction-work laboratory")
    st.markdown(r"""
For declared stress \(\sigma\) and transformation/deformation strain \(\varepsilon\):

\[
IW=\sigma:\varepsilon.
\]

With stress in MPa, the numerical result is MJ/m³. A larger positive value means the loading performs more positive mechanical work on that declared mode under this sign convention.
""")
    st.warning("IW is a mechanical driving-work term, not a complete kinetic law. Comparing deformation modes requires physically defensible strain/deformation tensors and barriers for those modes.")

    st.markdown("**Reference deformation gradient F₀** — editable. Default values are the published NiTi natural-OR example from CT-reference-author's correspondence treatment; do not use unchanged for another alloy.")
    F_default=np.array([[1.0251,-0.0589,0.1297],[-0.0589,1.0251,-0.1297],[-0.0016,0.0016,0.9511]])
    F=np.zeros((3,3))
    cols=st.columns(3)
    for i in range(3):
        for j in range(3):
            F[i,j]=cols[j].number_input(f"F{i+1}{j+1}",value=float(F_default[i,j]),format="%.6f",key=f"F{i}{j}")
    strain_kind=st.radio("Strain used for IW",["small strain: sym(F)-I","Green-Lagrange: (FᵀF-I)/2"],horizontal=True)
    E0=small_strain(F) if strain_kind.startswith("small") else green_lagrange(F)
    allE=variant_strains(E0,Gp)[:len(variants)]

    sigma=st.number_input("Uniaxial stress σ₀ (MPa)",value=100.0)
    n=np.array([st.number_input("loading n_x",value=1.0),st.number_input("loading n_y",value=0.0),st.number_input("loading n_z",value=0.0)])
    if np.linalg.norm(n)==0: n=np.array([1.,0.,0.])
    S=uniaxial_stress(sigma,n)
    iw=np.array([interaction_work(S,E) for E in allE])
    idf=pd.DataFrame({"variant":[f"V{i+1}" for i in range(len(iw))],"IW_MJ_m3":iw}).sort_values("IW_MJ_m3",ascending=False)
    st.dataframe(idf,use_container_width=True)
    st.metric("Highest IW for declared loading (MJ/m³)",f"{idf.IW_MJ_m3.iloc[0]:.6g}")
    st.write("Selected by pure IW ranking:",idf.variant.iloc[0])

    if st.button("Search IWmax over uniaxial loading directions"):
        mx,dirs=max_iw_uniaxial(sigma,allE,n_theta=25,n_phi=49)
        mdf=pd.DataFrame({"variant":[f"V{i+1}" for i in range(len(mx))],"IWmax_MJ_m3":mx,
                          "n_x":dirs[:,0],"n_y":dirs[:,1],"n_z":dirs[:,2]}).sort_values("IWmax_MJ_m3",ascending=False)
        st.dataframe(mdf,use_container_width=True)
        st.session_state["iwmax"]=mdf

with tabs[4]:
    st.subheader("Compatibility diagnostics")
    st.markdown(r"""
The current stable build calculates **diagnostics** from a supplied deformation gradient:

\[
U=\sqrt{F^TF},\qquad r_{\lambda_2}=|\lambda_2(U)-1|.
\]

A small \(r_{\lambda_2}\) is important in martensitic compatibility theory, but **it is not, by itself, proof of cofactor conditions or full CT supercompatibility**.
""")
    # use current/default F if IW tab has not been interacted with
    F=np.array([[1.0251,-0.0589,0.1297],[-0.0589,1.0251,-0.1297],[-0.0016,0.0016,0.9511]])
    res,lams=lambda2_residual(F)
    c1,c2,c3=st.columns(3)
    c1.metric("λ₁",f"{lams[0]:.8f}"); c2.metric("λ₂",f"{lams[1]:.8f}"); c3.metric("λ₃",f"{lams[2]:.8f}")
    st.metric("|λ₂ − 1|",f"{res:.6e}")
    st.metric("det(F) — volume ratio",f"{det_volume_ratio(F):.8f}")
    st.info("Advanced CT CMC/SMC + M/M twin + shear/shear supercompatibility belongs in a separately verified module. This build intentionally does not fabricate those equations.")

with tabs[5]:
    st.subheader("Reconstruction validation")
    st.markdown("Upload a CSV with columns `grain_id,true_parent,pred_parent`. Optional adjacency CSV: `i,j` using zero-based row indices.")
    up=st.file_uploader("Known-truth reconstruction CSV",type="csv",key="truth")
    if up:
        from mrs.validation import partition_metrics,boundary_metrics
        q=pd.read_csv(up)
        need={"grain_id","true_parent","pred_parent"}
        if need<=set(q.columns):
            m=partition_metrics(q.true_parent,q.pred_parent)
            st.dataframe(pd.DataFrame([m]),use_container_width=True)
            st.caption("These are truth-referenced metrics. Method-vs-method ARI/NMI must not be described as accuracy.")
            edges_up=st.file_uploader("Optional adjacency CSV",type="csv",key="edges")
            if edges_up:
                e=pd.read_csv(edges_up)
                bm=boundary_metrics(e[["i","j"]].to_numpy(int),q.true_parent.to_numpy(),q.pred_parent.to_numpy())
                st.dataframe(pd.DataFrame([bm]),use_container_width=True)
        else:
            st.error(f"Missing required columns: {sorted(need-set(q.columns))}")

with tabs[6]:
    st.subheader("Paper-ready evidence export")
    st.markdown("The evidence bundle records tables + exact input metadata + equation ledger + references + build ID. It is designed to make it obvious what was calculated and what was merely interpreted.")
    vrows=[]
    for i,R in enumerate(variants,1):
        e=bunge_euler(R); vrows.append({"variant":i,"phi1_deg":e[0],"Phi_deg":e[1],"phi2_deg":e[2],**{f"R{r+1}{c+1}":R[r,c] for r in range(3) for c in range(3)}})
    vdf=pd.DataFrame(vrows)
    odf=pd.DataFrame([{"operator":k,"double_coset_size":len(S)} for k,S in enumerate(op_sets)])
    adf2=pd.DataFrame(A,columns=[f"to_V{i+1}" for i in range(len(variants))]); adf2.insert(0,"from_variant",[f"V{i+1}" for i in range(len(variants))])
    crows=[{"left_operator":i,"right_operator":j,"result_operators":";".join(map(str,comp[(i,j)]))} for i in range(len(op_sets)) for j in range(len(op_sets))]
    tables={"variants":vdf,"operators":odf,"variant_operator_map":adf2,"operator_composition":pd.DataFrame(crows)}
    meta={"build_id":BUILD_ID,"software_version":__version__,"transformation":system,
          "lattice_parameters_A":{"a_B2":aP,"a_B19p":aM,"b_B19p":bM,"c_B19p":cM,"beta_deg":beta},
          "OR_definition":{"daughter_plane":"(010)","parent_plane":"(110)","daughter_direction":"[101]","parent_direction":"[-1 1 1]"},
          "symmetry":{"parent":"cubic proper rotational group 432 (24)","daughter":"monoclinic proper rotational group 2 (2)"},
          "variant_count":len(variants),"common_subgroup_order":len(H),"operator_class_count":len(op_sets)}
    blob=evidence_zip_bytes(tables,meta,EQS,REFS)
    st.download_button("Download academic evidence ZIP",blob,file_name="martensite_research_evidence.zip",mime="application/zip")
    with st.expander("Equation ledger"): st.markdown(EQS)
    with st.expander("References / attribution"): st.markdown(REFS)
