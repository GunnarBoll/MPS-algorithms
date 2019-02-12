
import time
import datetime
import pathlib
import imp
import os
# os.environ["OMP_NUM_THREADS"] = "1"
import multiprocessing as mp
import numpy as np

import storage as st
import ExactDiag as ed

imp.reload(np)
imp.reload(st)
imp.reload(ed)

def process_func(g1, g2, N, d, model, order, algo, parlist, ED_calc):
    if ED_calc:
        a = np.array([[0, 0], [1, 0]])
        adag = np.array([[0, 1], [0, 0]])
        ED = ed.ExactD(g1, g2, N, d, model, order)
        ED.exact_GS()
        data = [ED.E_GS]
        data += [ED.ED_correl(ED.GS, adag, a, ind1, ind2) for ind1 in range(N)
                 for ind2 in range(N)]
    else:
        Psi, H = run_algo(g1, g2, N, parlist[0], d, parlist[2], model, order,
                        parlist[1], algo)
        data = mk_dat(Psi, H, parlist[1])
    return data

def run_algo(g1, g2, N, dt, d, chi_max, model, order, T, algo):
    H = st.Hamiltonian(g1, g2, N, dt, d, chi_max, model, order,
                       grow_chi=False)
    Psi = st.StateChain(N, d, chi_max, algo)
    step_num = int(T / dt)
    H.time_evolve(Psi, step_num, algo)
    return Psi, H

def mk_dat(Psi, H, T):
    a = np.array([[0, 0], [1, 0]])
    adag = np.array([[0, 1], [0, 0]])
    E_GS = sum(Psi.get_ener(H.Hchain))
    meas = st.Measure()
    corr_mat = meas.corr_mat(Psi, adag, a)
    data = [E_GS, Psi.err]
    data += [corrs for corrs in corr_mat.reshape(H.N ** 2)]
    
    data += ["chi=" + str(H.chi_max) + ",T=" + str(T) + ",dt=" + str(H.dt)]
    return data
    
    
def main():
    
    
    if __name__ == '__main__':
        model = "HCboson"
        n_core = 2
        date = str(datetime.date.today())
        run_number = 1
        g1 = [1.0, 1.0]
        g2 = [1.0, 0.0]
        
        N = 6
        d = 2
        order = "fourth"
        algo = "tDMRG"
        dt_list = [0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1]
        T_list = [20, 24, 28, 30, 34, 40, 65, 100]
        T_list.reverse()
        chi_max_list = [10, 12, 14, 16, 18, 20, 30, 50, 70]
        chi_max_list.reverse()
        
        params = [[dt, T_list[0], chi_max_list[0]] for dt in dt_list]
        params += [[dt_list[0], T, chi_max_list[0]] for T in T_list[1:]]
        params += [[dt_list[0], T_list[0], chi_max] for chi_max in
                chi_max_list[1:]]
        arglist = [(g1, g2, N, d, model, order, algo, parl, False) for parl in
                params]        
        
        start = time.time()
        pool = mp.Pool(n_core)
        res = pool.starmap_async(process_func, arglist)
        pool.close()
        pool.join()
        end = time.time()
        
        print("Time cost:", end-start)
        
        dat_mat = [dat_list for dat_list in res.get()]
        
        date = str(datetime.date.today())
        run_number = 1
        direc = (os.getcwd() + "/" + "alf=" + str(g2[1]) + "_"  + date)
        while True:
            try:
                pathlib.Path(direc + "_" + str(run_number) + "/").mkdir(
                            parents=True, exist_ok=False)
                direc = direc + "_" + str(run_number) + "/"
                break
            except FileExistsError:
                run_number += 1
        
        for data in dat_mat:
            filename = data.pop()
            with open(direc+filename+".txt", "x") as fw:
                for dat_point in data:
                    fw.write(str(dat_point) + "\n")
        
        
    
    return

__spec__ = None
main()