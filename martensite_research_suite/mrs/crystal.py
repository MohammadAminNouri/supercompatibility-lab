import itertools, numpy as np
from .rotations import project_rotation

def cubic_proper_group():
    out=[]
    I=np.eye(3)
    for p in itertools.permutations(range(3)):
        P=I[:,p]
        for s in itertools.product([-1,1], repeat=3):
            R=P@np.diag(s)
            if round(np.linalg.det(R))==1 and not any(np.allclose(R,Q) for Q in out):
                out.append(R)
    assert len(out)==24
    return out

def monoclinic_2_group():
    return [np.eye(3), np.diag([-1.,1.,-1.])]

def lattice_basis(a,b,c,alpha=90,beta=90,gamma=90):
    ar,br,gr=np.radians([alpha,beta,gamma])
    va=np.array([a,0.,0.])
    vb=np.array([b*np.cos(gr), b*np.sin(gr),0.])
    cx=c*np.cos(br)
    cy=c*(np.cos(ar)-np.cos(br)*np.cos(gr))/np.sin(gr)
    cz=np.sqrt(max(c*c-cx*cx-cy*cy,0))
    return np.column_stack([va,vb,np.array([cx,cy,cz])])

def direction_cart(B,uvw):
    v=B@np.asarray(uvw,float); return v/np.linalg.norm(v)

def plane_normal_cart(B,hkl):
    n=np.linalg.inv(B).T@np.asarray(hkl,float); return n/np.linalg.norm(n)

def frame_from_plane_direction(B, hkl, uvw):
    n=plane_normal_cart(B,hkl)
    d=direction_cart(B,uvw)
    d=d-n*np.dot(d,n)
    d=d/np.linalg.norm(d)
    y=np.cross(n,d); y=y/np.linalg.norm(y)
    return np.column_stack([d,y,n])

def orientation_from_parallelisms(Bp,Bd,parent_plane,daughter_plane,parent_dir,daughter_dir):
    Fp=frame_from_plane_direction(Bp,parent_plane,parent_dir)
    Fd=frame_from_plane_direction(Bd,daughter_plane,daughter_dir)
    # maps daughter Cartesian coordinates into parent Cartesian coordinates
    return project_rotation(Fp @ Fd.T)

def niti_natural_or(a_b2=3.01,a_m=2.898,b_m=4.108,c_m=4.646,beta=97.78):
    Bp=lattice_basis(a_b2,a_b2,a_b2)
    Bd=lattice_basis(a_m,b_m,c_m,beta=beta)
    # AQ natural OR: (010)M // (110)P ; [101]M // [-1,1,1]P
    R=orientation_from_parallelisms(
        Bp,Bd,
        parent_plane=(1,1,0), daughter_plane=(0,1,0),
        parent_dir=(-1,1,1), daughter_dir=(1,0,1)
    )
    return R,Bp,Bd
