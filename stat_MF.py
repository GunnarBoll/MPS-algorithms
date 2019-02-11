"""Loop for self-consistent static mean-field calculation."""

import numpy as np
import importlib as imp
import datetime
import os
import pathlib
import multiprocessing as mp

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
    T = 150
    chi = 100
    
    # Args: (tperp, g1, g2, N, chi, T)
    ord_pars = SMF_loop(tperp, g1, g2, N, chi, T)        
    file_name = "N=" + str(N) + ",U=" + str(U)
    return ord_pars + [file_name]

def main():
    
    if __name__ == "__main__":
        N_list = [32]
        n_core = 5
        U_list = [0., 0.5, 1.0, 1.5, 2.]
        
        for N in N_list:
            
            
            pool = mp.Pool(n_core)
            res = pool.starmap_async(mk_dat, [(U, N) for U in U_list])
            pool.close()
            pool.join()
            
            dat_list = [dat for dat in res.get()]
            
            name = "SMF_" + "N=" + str(N)
            direc = os.getcwd() + "/" + name
            run_number = 1
            while True:
                try:
                    pathlib.Path(direc + "_" + str(run_number) + "/").mkdir(
                        parents=True, exist_ok=False)
                    direc += "_" + str(run_number) + "/"
                    break
                except FileExistsError:
                    run_number += 1
            
            for data in dat_list:
                filename = data.pop()
                with open(direc+filename+".txt", "x") as fw:
                    for data_point in data:
                        fw.write(str(data_point) + "\n")
    
    return
    
__spec__ = None
main()
