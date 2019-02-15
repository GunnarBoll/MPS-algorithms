import importlib as imp
import numpy as np
import sys
import time

import ExactDiag as ed
import cwd_storage as store

imp.reload(ed)
imp.reload(store)

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

def ED_MF_loop():
    start = time.time()
    U = float(sys.argv[1])
    N = int(sys.argv[2])
    tperp = float(sys.argv[3])
    ord_par = float(sys.argv[4])
    # mu = float(sys.argv[5])
    
    d = 2
    model = "HCboson"
    g1 = [1., U]
    g2 = [0. , 0.]
    mu_list = [-1., 0.5, 0., 0.5, 1.0, 1.5, 2.]
    
    g2[1] = 4 * tperp * ord_par
    
    dens_list = [get_dens(g1, [mu, g2[1]], N, d, model) for mu in mu_list]
    data = dens_list + mu_list
    
    direc_name = "rho_of_mu_tperp=" + str(tperp) + ",U=" + str(U)
    file_name = "ordpar=" + str(ord_par) + ".txt"
    store.cwd_store(direc_name, file_name, data)
    end = time.time()
    
    print("Time cost:", end-start)
    return
ED_MF_loop()