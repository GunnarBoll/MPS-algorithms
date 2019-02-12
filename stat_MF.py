"""Loop for self-consistent static mean-field calculation."""


import importlib as imp
import datetime
import os
import pathlib
import numpy as np
import sys

import storage as st
from static_MF_loop import SMF_loop

imp.reload(st)

def exp_func(x, a, b, c):
    return a*np.exp(-b*x) + c
    
def extrap_res(g1, g2, N, dt, d, chi, model, order, algo, step_num, a):
    M = st.Measure()
    H = st.Hamiltonian(g1, g2, N, dt, d, chi, model, TO=order,
                           grow_chi=False)
    Psi = st.StateChain(N, d, algo)
    Psi = H.time_evolve(Psi, step_num, algo)
        
    extr_ord = M.expec(Psi, a, int(N / 2))
    
    return extr_ord

def mk_dat(U, N):
    g1 = [1., U]
    g2 = [0., 0.01]
    tperp = 0.05
    T = 30
    chi = 10
    
    # Args: (tperp, g1, g2, N, chi, T)
    ord_pars, Psi = SMF_loop(tperp, g1, g2, N, chi, T)        
    file_name = "N=" + str(N) + ",U=" + str(U)
    return [ord_pars + [file_name], Psi]

def main():
    
    N = int(sys.argv[1])
    U = float(sys.argv[2])
    
    dat_list, Psi_list = mk_dat(U, N)
    
    print(Psi_list.B)
    print("\n", Psi_list.L)
    
    
    
    name = "SMF_" + "N=" + str(N)
    run_number = 1
    filename = dat_list.pop()
    GS_name = "GS" + filename
    while True:
        try:
            direc = os.getcwd() + "/" + name
            pathlib.Path(direc + "_" + str(run_number) + "/").mkdir(
                parents=True, exist_ok=True)
            direc += "_" + str(run_number) + "/"
            with open(direc+filename+".txt", "x") as fw:
                for data_point in dat_list:
                    fw.write(str(data_point) + "\n")
            break
        except FileExistsError:
            run_number += 1
    
        
    return
    
__spec__ = None
main()
