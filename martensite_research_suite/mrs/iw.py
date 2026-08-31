import numpy as np

def small_strain(F):
    F=np.asarray(F,float)
    return 0.5*(F+F.T)-np.eye(3)

def green_lagrange(F):
    F=np.asarray(F,float)
    return 0.5*(F.T@F-np.eye(3))

def uniaxial_stress(sigma_mpa, direction):
    n=np.asarray(direction,float); n=n/np.linalg.norm(n)
    return float(sigma_mpa)*np.outer(n,n)

def interaction_work(stress_mpa, strain):
    # MPa * dimensionless = MJ/m^3
    return float(np.tensordot(np.asarray(stress_mpa,float),np.asarray(strain,float),axes=2))

def variant_strains(E0, parent_symmetry):
    return [g@E0@g.T for g in parent_symmetry]

def iw_table(stress, strains):
    return np.array([interaction_work(stress,E) for E in strains])

def unit_sphere_directions(n_theta=37,n_phi=73):
    out=[]
    for th in np.linspace(0,np.pi,n_theta):
        for ph in np.linspace(0,2*np.pi,n_phi,endpoint=False):
            out.append([np.sin(th)*np.cos(ph),np.sin(th)*np.sin(ph),np.cos(th)])
    return np.asarray(out)

def max_iw_uniaxial(sigma_mpa,strains,n_theta=37,n_phi=73):
    dirs=unit_sphere_directions(n_theta,n_phi)
    vals=np.empty((len(strains),len(dirs)))
    for i,E in enumerate(strains):
        vals[i]=[interaction_work(uniaxial_stress(sigma_mpa,n),E) for n in dirs]
    imax=np.argmax(vals,axis=1)
    return vals[np.arange(len(strains)),imax], dirs[imax]
