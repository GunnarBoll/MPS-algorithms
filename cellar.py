"""Unused (possible broken) code"""

import numpy as np

#Algorithm implementing Hastings trick (Caution: Incorrect!)
def hastings(B,Lam,Ut,i, chi, d):
    chia = B[i].shape[1]
    chib = B[i+1].shape[2]
    
    thetabar = np.tensordot(B[i],B[i+1],(2,1))
    
    phibar = np.transpose(np.tensordot(Ut[i], thetabar, ([1,3], [0,2])), (0,2,1,3))
    
    phi = np.transpose(np.tensordot(np.diag(Lam[i]),phibar, (1,1)), (1,0,2,3))
    phi = np.reshape(phi, (d*chia,d*chib))
    
    U, S, V = np.linalg.svd(phi, full_matrices = True, compute_uv = True)
    V = V.T

    chic = np.min([np.sum(S>10**-12), chi])
    
    U = np.reshape(U[:d*chia,:chic], (d, chia, chic))
    V = np.transpose(np.reshape(V[:d*chib,:chic], (d, chib, chic)), (0,2,1))
    
    #State update
    B[i+1] = V
    Lam[i+1] = S[:chic]/np.sqrt(sum(S[:chic]**2))
    B[i] = np.tensordot(phibar,np.transpose(np.conj(B[i+1]), (0,2,1)), ([2,3],[0,1]))
    return B, Lam
    
#Function for getting the energy of a state given using Hastings trick
def get_hastings(Lam, B, Hchain, N):
    E = []
    for c in range(0,N-1):
        LB = np.tensordot(np.diag(Lam[c]),B[c], (1,1))
        LBB = np.tensordot(LB, B[c+1], (2,1))
        C = np.tensordot(LBB, np.reshape(Hchain[c], (2,2,2,2)), ([1,2],[2,3]))
        E.append(np.tensordot(np.conj(LBB), C, ([0,3,1,2],[0,1,2,3]) ))
    return E