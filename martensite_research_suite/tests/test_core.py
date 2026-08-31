import numpy as np
from mrs.crystal import cubic_proper_group,monoclinic_2_group,niti_natural_or
from mrs.groupoid import common_subgroup,variants_from_cosets,operator_classes,composition_table,arrow_operator_matrix
from mrs.iw import interaction_work,uniaxial_stress,small_strain
from mrs.compatibility import stretch_eigenvalues

def test_cubic_group():
    G=cubic_proper_group()
    assert len(G)==24
    assert all(abs(np.linalg.det(g)-1)<1e-12 for g in G)

def test_niti_variants_are_12():
    Gp=cubic_proper_group(); Gd=monoclinic_2_group()
    R,_,_=niti_natural_or()
    H=common_subgroup(Gp,Gd,R,tol_deg=1e-4)
    V,reps,C=variants_from_cosets(Gp,H,R)
    assert len(H)==2
    assert len(V)==12

def test_operator_table_internal():
    Gp=cubic_proper_group(); Gd=monoclinic_2_group()
    R,_,_=niti_natural_or()
    H=common_subgroup(Gp,Gd,R,tol_deg=1e-4)
    V,reps,C=variants_from_cosets(Gp,H,R)
    ops,idx=operator_classes(Gp,H)
    A=arrow_operator_matrix(Gp,reps,idx)
    T=composition_table(Gp,ops,idx)
    assert A.shape==(12,12)
    assert all(A[i,i]==0 for i in range(12))
    assert len(T)==len(ops)**2

def test_iw_units_and_sign():
    E=np.diag([0.01,0,0])
    S=uniaxial_stress(100,[1,0,0])
    assert abs(interaction_work(S,E)-1.0)<1e-12

def test_stretch_sorted_positive():
    F=np.diag([.95,1,1.1])
    l=stretch_eigenvalues(F)
    assert np.allclose(l,[.95,1,1.1])
