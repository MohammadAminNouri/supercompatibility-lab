import numpy as np
from .rotations import angle_axis

def mat_key(M, nd=8):
    return tuple(np.round(np.asarray(M,float).ravel(),nd))

def contains(group,M,tol=1e-7):
    return any(np.linalg.norm(Q-M)<tol for Q in group)

def common_subgroup(Gp,Gd,R,tol_deg=1e-5):
    H=[]
    for gp in Gp:
        best=min(angle_axis(gp.T@(R@gd@R.T))[0] for gd in Gd)
        if best<tol_deg:
            H.append(gp)
    return H

def left_cosets(G,H,tol=1e-8):
    unused=list(range(len(G))); cosets=[]
    while unused:
        i=unused[0]; g=G[i]; C=[]
        for h in H:
            gh=g@h
            j=min(range(len(G)),key=lambda k:np.linalg.norm(G[k]-gh))
            if np.linalg.norm(G[j]-gh)<tol: C.append(j)
        C=sorted(set(C)); cosets.append(C); unused=[u for u in unused if u not in C]
    return cosets

def double_coset_indices(G,H,g_idx,tol=1e-8):
    S=set()
    g=G[g_idx]
    for h1 in H:
        for h2 in H:
            x=h1@g@h2
            j=min(range(len(G)),key=lambda k:np.linalg.norm(G[k]-x))
            if np.linalg.norm(G[j]-x)<tol: S.add(j)
    return frozenset(S)

def operator_classes(G,H):
    # Force the identity double-coset to be O0 for stable scientific labeling.
    identity_idx=min(range(len(G)), key=lambda k: np.linalg.norm(G[k]-np.eye(3)))
    identity_set=double_coset_indices(G,H,identity_idx)
    op_sets=[identity_set]
    for i in range(len(G)):
        D=double_coset_indices(G,H,i)
        if D not in op_sets:
            op_sets.append(D)
    index_to_op={}
    for i in range(len(G)):
        D=double_coset_indices(G,H,i)
        index_to_op[i]=op_sets.index(D)
    return op_sets,index_to_op

def group_mul_index(G,A,B):
    X=A@B
    j=min(range(len(G)),key=lambda k:np.linalg.norm(G[k]-X))
    return j

def composition_table(G,op_sets,index_to_op):
    # CT-reference-author operator composition convention: (Om,On) -> Om^{-1} On.
    table={}
    for m,Sm in enumerate(op_sets):
        for n,Sn in enumerate(op_sets):
            outs=set()
            for i in Sm:
                for j in Sn:
                    k=group_mul_index(G,G[i].T,G[j])
                    outs.add(index_to_op[k])
            table[(m,n)]=tuple(sorted(outs))
    return table

def variants_from_cosets(Gp,H,R0):
    cosets=left_cosets(Gp,H)
    reps=[Gp[C[0]] for C in cosets]
    R=[g@R0 for g in reps]
    return R,reps,cosets

def arrow_operator_matrix(G,reps,index_to_op):
    N=len(reps); A=np.zeros((N,N),int)
    for i in range(N):
        for j in range(N):
            k=group_mul_index(G,reps[i].T,reps[j])
            A[i,j]=index_to_op[k]
    return A
