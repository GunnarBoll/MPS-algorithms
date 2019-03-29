"""
Program for self-consistent static mean-field calculation. The actual loop is
contained in static_MF_loop.py.
"""

import importlib as imp
import os
import pathlib
import numpy as np
import sys
import time

import storage as st
from static_MF_loop import SMF_loop

imp.reload(st)

# Fittable exponential function
def exp_func(x, a, b, c):
    return a*np.exp(-b*x) + c

# Function for extrapolating the order parameter without continuing the loop
def extrap_res(g1, g2, N, dt, d, chi, model, order, algo, step_num, a):
    M = st.Measure()
    H = st.Hamiltonian(g1, g2, N, dt, d, chi, model, TO=order,
                           grow_chi=False)
    Psi = st.StateChain(N, d, algo)
    Psi = H.time_evolve(Psi, step_num, algo)
        
    extr_ord = M.expec(Psi, a, int(N / 2))
    
    return extr_ord

# Calls the main loop function for parameter U and N. Returns a list with
# the order parameters from the loop and a file name for them, the ground state
# Psi and a list of chemical potentials.
def mk_dat(U, N):
    g1 = [1., U]
    g2 = [0., 0.01]
    tperp = 0.05
    T = 200
    chi = 50
    
    # Args: (tperp, g1, g2, N, chi, T)
    ord_pars, Psi, mus = SMF_loop(tperp, g1, g2, N, chi, T)   
    file_name = "N=" + str(N) + ",U=" + str(U)
    return [ord_pars + [file_name], Psi, mus]

# Main program, takes arguments N and U from console. Calls a data collection
# algorithm and saves the order parameters, chemical potentials and ground
# state in files.
def main():
    start = time.time()
    
    N = int(sys.argv[1])
    U = float(sys.argv[2])
    
    print("U is:", U)
    dat_list, Psi, mu_dat = mk_dat(U, N)
    
    # Data for ground state
    GS_mat = ([coup for coup in Psi.g1] + [coup2 for coup2 in Psi.g2]
              + [Psi.N, Psi.d, Psi.err, Psi.notation, Psi.chi])
    for tens in Psi.B:
        GS_mat += [elem for elem in 
                   tens.reshape(tens.shape[0]*tens.shape[1]*tens.shape[2])]
    
    for lams in Psi.L:
        GS_mat += [lam_elem for lam_elem in lams]
    
    
    name = "SMF_" + "N=" + str(N)
    run_number = 1
    filename = dat_list.pop()
    GS_name = "GS_" + filename
    while True:
        try:
            direc = os.getcwd() + "/" + name
            pathlib.Path(direc + "_" + str(run_number) + "/").mkdir(
                parents=True, exist_ok=True)
            direc += "_" + str(run_number) + "/"
            
            # Create file for order parameters
            with open(direc+filename+".txt", "x") as fw:
                for data_point in dat_list:
                    fw.write(str(data_point) + "\n")
            
            # Create file for chemical potentials
            with open(direc+"mu_"+filename+".txt", "x") as fw:
                for mu in mu_dat:
                    fw.write(str(mu) + "\n")
            
            # Create file for ground state
            with open(direc+GS_name+".txt", "x") as fw:
                for comp in GS_mat:
                    fw.write(str(comp) + "\n")
            break
        except FileExistsError:
            run_number += 1
    
    end = time.time()
    print("Time cost:", end-start)
    return
    
__spec__ = None
main()
