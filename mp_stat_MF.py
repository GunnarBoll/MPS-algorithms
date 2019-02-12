"""Loop for self-consistent static mean-field calculation."""


import importlib as imp
import datetime
import multiprocessing as mp
import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["VECLIB_NUM_THREADS"] = "1"
import pathlib
import numpy as np

import storage as st
from static_MF_loop import SMF_loop

imp.reload(st)

os.system("taskset -p 0xff %d" % os.getpid())

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
    
    if __name__ == "__main__":
        N_list = [4]
        n_core = 2
        U_list = [0.,0.5]
        
        for N in N_list:
            
            
            pool = mp.Pool(n_core)
            res = pool.starmap_async(mk_dat, [(U, N) for U in U_list])
            pool.close()
            pool.join()
            
            dat_list = [dat[0] for dat in res.get()]
            Psi_list = [dat[1] for dat in res.get()]
            
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
