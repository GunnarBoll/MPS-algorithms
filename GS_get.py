import numpy as np
import os
import importlib as imp

import storage as st

imp.reload(st)

def get_GS(filename):
    
    GS_mat = []
    
    with open(filename, "r") as fr:
        GS_mat += [line.strip("\n") for line in fr]
    
    N = int(GS_mat.pop(0))
    d = int(GS_mat.pop(0))
    err = float(GS_mat.pop(0))
    notation = GS_mat.pop(0)
    chi = int(GS_mat.pop(0))
    
    GS_mat = [float(dat) for dat in GS_mat]
    tens_net = []
    lam_net = []
    
    for i in range(int(N / 2)):
        shape1 = d**i if d**i < chi else chi
        shape2 = d**(i+1) if d**(i+1) < chi else chi
        GS_mat, tens_comps = GS_mat[d*shape1*shape2:], GS_mat[:d*shape1*shape2]
        tens_net += [np.array(tens_comps).reshape(d, shape1, shape2)]
    
    for j in range(int(N / 2), 0, -1):
        shape1 = d**j if d**j < chi else chi
        shape2 = d**(j-1) if d**(j-1) < chi else chi
        GS_mat, tens_comps = GS_mat[d*shape1*shape2:], GS_mat[:d*shape1*shape2]
        tens_net += [np.array(tens_comps).reshape(d, shape1, shape2)]
    
    for k in range(N+1):
        np_dat = min([d**k, d**(N-k), chi])
        GS_mat, lam_comps = GS_mat[np_dat:], GS_mat[:np_dat]
        lam_net += [np.array(lam_comps)]
    
    indata = [err, notation, lam_net, tens_net]
    Psi = st.StateChain(N, d, chi, "tDMRG", load_state=True, load_mat=indata)
    
    return Psi