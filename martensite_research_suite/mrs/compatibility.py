import numpy as np

def right_stretch(F):
    C=np.asarray(F,float).T@np.asarray(F,float)
    w,V=np.linalg.eigh(C)
    U=V@np.diag(np.sqrt(np.clip(w,0,None)))@V.T
    return U

def stretch_eigenvalues(F):
    return np.linalg.eigvalsh(right_stretch(F))

def lambda2_residual(F):
    l=stretch_eigenvalues(F)
    return float(abs(l[1]-1.0)), l

def rank_one_residual(A,B):
    # distance of A-B from rank-one matrix, Frobenius norm from best SVD rank-one approximation
    D=np.asarray(A,float)-np.asarray(B,float)
    s=np.linalg.svd(D,compute_uv=False)
    return float(np.sqrt(np.sum(s[1:]**2)))

def det_volume_ratio(F):
    return float(np.linalg.det(np.asarray(F,float)))
