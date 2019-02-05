"""Loop for self-consistent static mean-field calculation."""

import numpy as np
import importlib as imp
import datetime
import os
import pathlib

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

def main():
    name = "SMF_" + str(datetime.date.today())
    direc = os.getcwd() + "/" + name
    run_number = 1
    while True:
        try:
            pathlib.Path(direc + "_" + str(run_number) + "/").mkdir(
                parents=True, exist_ok=False)
            break
        except FileExistsError:
            run_number += 1
    
    direc += "_" + str(run_number) + "/"
    
    N_list = [20]
    U_list = [1.5,2.]
    
    for N in N_list:
        for U in U_list:
            g1 = [1., U]
            g2 = [0., 0.01]
            tperp = 0.05
            T = 100
            chi = 70
            
            # Args: (tperp, g1, g2, N, chi, T)
            ord_pars = SMF_loop(tperp, g1, g2, N, chi, T)        
            file_name = "N=" + str(N) + ",U=" + str(U)
            with open(direc+file_name+".txt", "x") as fw:
                for ord_par in ord_pars:
                    fw.write(str(ord_par) + "\n")
    
    return
    
main()
