
import numpy as np
import time
import datetime
import pathlib
import imp

import storage as st
import ExactDiag as ed

imp.reload(st)
imp.reload(ed)

def run_algo(g1, g2, N, dt, d, chi_max, model, order, T, algo, bis_err):
    H = st.Hamiltonian(g1, g2, N, dt, d, chi_max, model, order)
    Psi = st.StateChain(N, d, algo, bis_err)
    step_num = int(T / dt)
    H.time_evolve(Psi, step_num, algo)
    ED = ed.ExactD(g1, g2, N, dt, d, model)
    ED.exact_GS()
    return Psi, H, ED

def filewrite(Psi, H, ED, direc, T):
    writing_mat = []
    a = np.array([[0, 1], [0, 0]])
    adag = np.array([[1, 0], [0, 0]])
    E_GS = sum(Psi.get_ener(H.Hchain))
    writing_mat.append(E_GS)
    writing_mat.append(ED.E_GS)
    corr_algo = []
    ED_corr = []
    meas = st.Measure()
    for ind1 in range(H.N):
        for ind2 in range(H.N):
            corr_algo.append(meas.correl(Psi, adag, a, ind1, ind2))
            ED_corr.append(ED.ED_correl(ED.GS, adag, a, ind1, ind2))
    writing_mat += corr_algo + ED_corr
    
    name = ("chi=" + str(H.chi_max) + ",T=" + str(T) + ",dt=" + str(H.dt)
            + ",BE=" + str(Psi.bis_err))
    with open(direc+name+".txt",'x') as fw:
        for data in writing_mat:
            fw.write(str(data)+"\n")
    fw.close()

def main():
    model = "HCboson"
    date = str(datetime.date.today())
    run_number = 1
    direc = ("C:/Users/Gunnar/Documents/Ph.D/Learning/DMRG/Tryout code/output/"
             + date)
    while True:
        try:
            pathlib.Path(direc + "_run#" + str(run_number) + "/").mkdir(
                parents=True, exist_ok=False)
            break
        except FileExistsError:
            run_number += 1
    direc = direc + "_run#" + str(run_number) + "/"
    
    g1 = [1.0, 0.0]
    g2 = [0.0, 1.0]
    dt_list = [0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2]
    T_list = [20, 30, 40, 50, 60, 70, 80, 90, 100]
    T_list.reverse()
    chi_max_list = [20, 25, 30, 35, 40, 45, 50]
    chi_max_list.reverse()
    bis_err_list = [10**-6, 10**-7, 10**-8, 10**-9, 10**-10]
    bis_err_list.reverse()
    order = "fourth"
    algo = "tDMRG"
    d = 2
    N = 4
    
    start = time.process_time()
    iter_lists = [dt_list, T_list[1:], chi_max_list[1:], bis_err_list[1:]]
    for ind in range(4):
        params = [dt_list[0], T_list[0], chi_max_list[0], bis_err_list[0]]
        param_list = iter_lists[ind]
        for params[ind] in param_list:
            dt = params[0]
            T = params[1]
            chi_max = params[2]
            bis_err = params[3]
            Psi, H, ED = run_algo(g1, g2, N, dt, d, chi_max, model, order, T,
                                  algo, bis_err)
            filewrite(Psi, H, ED, direc, T)
        
    end = time.process_time()
    print("Time:", end-start)
    return 0

main()