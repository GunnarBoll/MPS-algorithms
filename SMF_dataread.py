from matplotlib import pyplot as plt
import importlib as imp
import numpy as np

import storage as st
from GS_get import get_GS

imp.reload(st)

def tail(filename):
    
    with open(filename, 'r') as fr:
        lines = fr.read().splitlines()
        last_line = lines[-1]
    
    return float(last_line)

def main():
    data_direc = "C:/Users/Gunnar/Documents/Ph.D/Data/Static_MF/"
    N_list = [20, 24, 32]
    for N in N_list:
        name = "SMF_"+"N="+str(N)+"_1/"
        direc = data_direc + name
        U_list = [0., 0.5, 1., 1.5, 2.]
        
        files = [direc+"N="+str(N)+",U="+str(U)+".txt" for U in U_list]
        ord_par = [tail(file) for file in files]
                
        plt.plot(U_list, ord_par)
    plt.show()
    
    N = 32
    dt = 0.1
    d = 2
    chi = 100
    alf = ord_par[0]*4*0.05
    g1 = [1.0, 0.]
    g2 = [0.0, alf]
    
    H = st.Hamiltonian(g1, g2, N, dt, d, chi, "tDMRG", "fourth")
    
    Psi = get_GS(data_direc+"SMF_N=32_1/GS_N=32,U=2.0.txt")
    
    a = np.array([[0, 0], [1, 0]])
    adag = np.array([[0, 1], [0, 0]])
    num_op = np.matmul(adag, a)
    
    M = st.Measure()
    
    dens = 0
    for k in range(N):
        dens += M.expec(Psi, num_op, k)
    print(dens)
    return

main()