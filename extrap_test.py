import numpy as np
import imp

import storage as st

imp.reload(st)

def main():
    
    g1 = [1, 1]
    g2 = [1, 0.01]
    dt = 0.1
    N = 20
    chi_list = [10, 20]
    d = 2
    T = 100
    
    date = "2019-01-10"
    run_num = "_run#1/"
    direc = ("C:/Users/Gunnar/Documents/Ph.D/Learning/DMRG/Tryout code/output/"
             + date + run_num)
    
    with open(direc + "ED.txt", "r") as fed:
        E_GS = float(fed.readline())
    
    E_algo = []
    Err_algo = []
    for chi in chi_list:
        H = st.Hamiltonian(g1, g2, N, dt, d, chi, "HCboson", "fourth",
                           grow_chi=False)
        Psi = st.StateChain(N, d, "tDMRG")
        Psi = H.time_evolve(Psi, int(T / dt), "tDMRG")
        E_algo.append(sum(Psi.get_ener(H.Hchain)))
        Err_algo.append(Psi.err)
    
    p = np.polyfit(Err_algo, E_algo, 1)
    
    print(p[1])
    print(E_algo[1])
    print(E_GS)
    
    return
    
    
main()