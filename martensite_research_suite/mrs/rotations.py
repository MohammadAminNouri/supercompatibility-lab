import numpy as np

def project_rotation(M):
    U, _, Vt = np.linalg.svd(np.asarray(M, float))
    R = U @ Vt
    if np.linalg.det(R) < 0:
        U[:, -1] *= -1
        R = U @ Vt
    return R

def angle_axis(R, degrees=True):
    R = project_rotation(R)
    c = np.clip((np.trace(R)-1)/2, -1.0, 1.0)
    a = np.arccos(c)
    if abs(a) < 1e-12:
        axis = np.array([1.,0.,0.])
    else:
        axis = np.array([R[2,1]-R[1,2], R[0,2]-R[2,0], R[1,0]-R[0,1]])/(2*np.sin(a))
        n=np.linalg.norm(axis)
        axis = axis/n if n else np.array([1.,0.,0.])
    return (np.degrees(a) if degrees else a), axis

def rot_axis_angle(axis, angle_deg):
    a=np.asarray(axis,float); a=a/np.linalg.norm(a)
    t=np.radians(angle_deg)
    K=np.array([[0,-a[2],a[1]],[a[2],0,-a[0]],[-a[1],a[0],0]])
    return np.eye(3)+np.sin(t)*K+(1-np.cos(t))*(K@K)

def misorientation_angle(R1,R2):
    return angle_axis(np.asarray(R1).T @ np.asarray(R2))[0]

def bunge_euler(R):
    R=project_rotation(R)
    Phi=np.arccos(np.clip(R[2,2],-1,1))
    if abs(np.sin(Phi))>1e-10:
        phi1=np.arctan2(R[2,0],-R[2,1])
        phi2=np.arctan2(R[0,2], R[1,2])
    else:
        phi1=np.arctan2(R[0,1],R[0,0]); phi2=0.0
    return tuple(np.degrees([phi1%(2*np.pi),Phi,phi2%(2*np.pi)]))
