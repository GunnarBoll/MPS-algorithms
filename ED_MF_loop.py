import importlib as imp
import numpy as np
import sys

import ExactDiag as ed
from cwd_storage import cwd_store

imp.reload(ed)

def get_dens(mu, coup1, coup2, N, d, model):
    a = np.array([[0, 0], [1, 0]])
    adag = np.array([[0, 1], [0, 0]])
    ED = ed.ExactD(coup1, coup2, N, d, model)
    ED.exact_GS()
    
    dens = 0
    for k in range(N):
        dens += ED.ED_correl(ED.GS, adag, a, k, k)
    dens /= N
    
    return dens

def ED_MF_loop(U, N, tperp, ord_par):
    U, N, tperp, ord_par = sys.argv[1:5]
    
    
    d = 2
    model = "HCboson"
    g1 = [1., U]
    g2 = [0. , 0.]
    
    a = np.array([[0, 0], [1, 0]])
    adag = np.array([[0, 1], [0, 0]])
    g2[1] = 4 * tperp * ord_par
    mu_list = [-1., -0.5, 0.0, 0.5, 1.0, 1.5, 2.0]
    
    ED = ed.ExactD(g1, g2, N, d, model)
    ED.exact_GS()
    
    new_dens = 0
    for k in range(N):
        new_dens += ED.ED_correl(ED.GS, adag, a, k, k)
    new_dens /= N
    
    
    dens_list = [get_dens(mu) for mu in mu_list]
    data = dens_list + mu_list
    
    direc_name = "rho_of_mu_tperp=" + str(tperp) + ",U=" + str(U)
    file_name = "ordpar=" + str(ord_par) + ".txt"
    cwd_store(direc_name, file_name, data)
    
    return