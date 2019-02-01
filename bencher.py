"""
Program creates data for benchmarking the MPS algorithm written and stored in
storage.py.
"""
import numpy as np
import time
import datetime
import pathlib
import imp
import os

import storage as st
import ExactDiag as ed

imp.reload(st)
imp.reload(ed)

def run_algo(g1, g2, N, dt, d, chi_max, model, order, T, algo, bis_err,
             trunc_err_check):
    H = st.Hamiltonian(g1, g2, N, dt, d, chi_max, model, order,
                       grow_chi=trunc_err_check)
    Psi = st.StateChain(N, d, algo, bis_err)
    step_num = int(T / dt)
    H.time_evolve(Psi, step_num, algo)
    return Psi, H

def filewrite(Psi, H, direc, T):
    a = np.array([[0, 0], [1, 0]])
    adag = np.array([[0, 1], [0, 0]])
    E_GS = sum(Psi.get_ener(H.Hchain))
    meas = st.Measure()
    corr_mat = meas.corr_mat(Psi, adag, a)
    
    name = ("chi=" + str(H.chi_max) + ",T=" + str(T) + ",dt=" + str(H.dt)
            + ",BE=" + str(Psi.bis_err))
    with open(direc+name+".txt",'x') as fw:
        fw.write(str(E_GS) + "\n")
        fw.write(str(Psi.err) + "\n")
        for correlations in corr_mat.reshape(Psi.N ** 2):
            fw.write(str(correlations) + "\n")
    fw.close()

def main():
    model = "HCboson"
    date = str(datetime.date.today())
    run_number = 1
    direc = (os.getcwd() + "/" + date)
    while True:
        try:
            pathlib.Path(direc + "_run#" + str(run_number) + "/").mkdir(
                parents=True, exist_ok=False)
            break
        except FileExistsError:
            run_number += 1
    direc = direc + "_run#" + str(run_number) + "/"
    
    g1 = [1.0, 1.0]
    g2 = [1.0, 0.01]
    dt_list = [0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1]
    T_list = [20, 24, 28, 30, 34, 40, 65, 100]
    T_list.reverse()
    chi_max_list = [10, 12, 14, 16, 18, 20, 30, 50, 70]
    chi_max_list.reverse()
    bis_err_list = [10**-7, 5*10**-8, 10**-8, 5*10**-9, 10**-9, 5*10**-10,
                    10**-10, 5*10**-11, 10**-11]
    bis_err_list.reverse()
    order = "fourth"
    algo = "tDMRG"
    d = 2
    N = 4
    
    start = time.process_time()
    start2 = time.time()
    ed_file = direc+"ED.txt"
    with open(ed_file, "x") as fed:
        a = np.array([[0, 0], [1, 0]])
        adag = np.array([[0, 1], [0, 0]])
        ED = ed.ExactD(g1, g2, N, d, model, order)
        write_mat = []
        ED.exact_GS()
        write_mat.append(ED.E_GS)
        for ind1 in range(N):
            for ind2 in range(N):
                write_mat.append(ED.ED_correl(ED.GS, adag, a, ind1, ind2))
        for data in write_mat:
            fed.write(str(data) + "\n")
        fed.close()
    
    iter_lists = [dt_list, T_list[1:], chi_max_list[1:], bis_err_list[1:]]
    trunc_err_check = False
    for ind in range(3):
        if ind == 3:
            trunc_err_check = True
        params = [dt_list[0], T_list[0], chi_max_list[0], bis_err_list[0]]
        param_list = iter_lists[ind]
        for params[ind] in param_list:
            dt = params[0]
            T = params[1]
            chi_max = params[2]
            bis_err = params[3]
            Psi, H = run_algo(g1, g2, N, dt, d, chi_max, model, order, T,
                                  algo, bis_err, trunc_err_check)
            filewrite(Psi, H, direc, T)
        
    end = time.process_time()
    end2 = time.time()
    print("Time:", end-start)
    print("Real time:", end2-start2)
    return 0

main()