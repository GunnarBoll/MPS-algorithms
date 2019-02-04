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
    U_list = [0]
    
    for N in N_list:
        for U in U_list:
            g1 = [1., U]
            g2 = [0., 0.01]
            tperp = 0.05
            T = 100
            chi = 70
            
            # dt = 0.1
            # step_num = int(T / dt)
            # d = 2
            # model = "HCboson"
            # order = "fourth"
            # algo = "tDMRG"
            
            #Args: (tperp, g1, g2, N, chi, T)
            ord_pars = SMF_loop(tperp, g1, g2, N, chi, T)
            # M = st.Measure()
            # a = np.array([[0, 0], [1, 0]])
            # adag = np.array([[0, 1], [0, 0]])
            # num_op = np.matmul(adag, a)
            # a_exp_guess = 1 / np.sqrt(N / 2)
            # i = 0
            # err = 1
            # g2[1] = 4 * tperp * a_exp_guess
            # ord_pars = [a_exp_guess]
            # mu_list = [g2[0]]
            # 
            # while i < 10 and err > 10 ** -4:
            #     H = st.Hamiltonian(g1, g2, N, dt, d, chi, model, TO=order,
            #                     grow_chi=False)
            #     Psi = st.StateChain(N, d, algo)
            #     Psi = H.time_evolve(Psi, step_num, algo)
            #     
            #     if i == 0:
            #         dens = 0
            #         for k in range(N):
            #             dens += M.expec(Psi, num_op, k)
            #     else:
            #         new_dens = 0
            #         for k in range(N):
            #             new_dens += M.expec(Psi, num_op, k)
            #         g2[0] = g2[0]*dens / new_dens
            #         dens = new_dens
            #         mu_list.append(g2[0])
            #         
            #     new_ord_par = M.expec(Psi, a, int(N / 2))
            #     err = abs((abs(ord_pars[i]) - abs(new_ord_par)) / abs(ord_pars[i]))
            #     g2[1] = 4*new_ord_par*tperp
            #     i += 1
            #     ord_pars.append(abs(new_ord_par))
            
            file_name = "N=" + str(N) + ",U=" + str(U)
            with open(direc+file_name+".txt", "x") as fw:
                for ord_par in ord_pars:
                    fw.write(str(ord_par) + "\n")
    
    
    
    
    
    return
    
main()