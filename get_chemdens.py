"""
Main program constructs reference data for the static MF algorithm. The data is
found "exactly" using exact diagonalization. For the reference data to be used
a recompiling algorithm (refdata_comp.py) must be used.
"""
import importlib as imp
import numpy as np
import sys
import time

import ExactDiag as ed
import cwd_storage as store

imp.reload(ed)
imp.reload(store)

# Given the model parameters returns the system particle density
def get_dens(coup1, coup2, N, d, model):
    a = np.array([[0, 0], [1, 0]])
    adag = np.array([[0, 1], [0, 0]])
    ED = ed.ExactD(coup1, coup2, N, d, model)
    ED.exact_GS()
    
    dens = 0
    for k in range(N):
        dens += ED.ED_correl(ED.GS, adag, a, k, k)
    dens /= N
    
    return dens

# Main program
def ED_MF_refdata():
    start = time.time()
    
    # Takes arguments from console
    U = float(sys.argv[1])
    N = int(sys.argv[2])
    tperp = float(sys.argv[3])
    ord_par = float(sys.argv[4])
    mu = float(sys.argv[5])
    
    d = 2
    model = "HCboson"
    g1 = [1., U]
    g2 = [0. , 0.]
    
    g2[1] = 4 * tperp * ord_par
    
    # Gets a density given a chemical potential mu and stores them in a file
    dens = get_dens(g1, [mu, g2[1]], N, d, model)
    data = [dens, mu]
    
    direc_name = "rho_of_mu_tperp=" + str(tperp) + ",U=" + str(U) + "/"
    file_name = "ordpar=" + str(ord_par) + ",mu=" + str(mu) + ".txt"
    store.cwd_store(direc_name, file_name, data)
    end = time.time()
    
    print("Time cost:", end-start)
    return
ED_MF_refdata()