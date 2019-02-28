import numpy as np
import importlib as imp
import os

import storage as st
import ExactDiag as ed

imp.reload(st)
imp.reload(ed)

def new_mu(coup1, coup2, N, dt, d, chi, T, num_op, old_mu, start_dens,
           rho_maxerr=1e-4):
    goal_dens = 1 / 2
    Meas = st.Measure()
    order = "fourth"
    model = "HCboson"
    if start_dens < goal_dens:
        mu0 = old_mu
        dens0 = start_dens
    else:
        mu1 = old_mu
        dens1 = start_dens
    for ind in range(20):
        mu_ham = st.Hamiltonian(coup1, coup2, N, dt, d, chi, model, TO=order,
                                grow_chi=False)
        mu_psi = st.StateChain(N, d, chi, "tDMRG")
        
        steps = int(T / dt)
        #for ind2 in range(5):
        mu_psi = mu_ham.time_evolve(mu_psi, steps, "tDMRG", fast_run=True)
        
        mu_dens = 0
        for m in range(N):
            mu_dens += Meas.expec(mu_psi, num_op, m)
        mu_dens /= N
        
        dens_err = (mu_dens - goal_dens) / goal_dens
        
        if mu_dens > goal_dens:
            mu1 = coup2[0]
            dens1 = mu_dens
        else:
            mu0 = coup2[0]
            dens0 = mu_dens
        
        if abs(dens_err) < rho_maxerr or ind == 19:
            mu = coup2[0]
            if ind == 19:
                print("Density error:", dens_err)
            break
        else:
            yp = ([mu0, mu1] if mu0 < mu1 else [mu1, mu0])
            xp = ([dens0, dens1] if mu0 < mu1 else [dens1, dens0])

            coup2[0] = np.interp(goal_dens, xp, yp)
                    
    return mu

def guess_mu(ord_par, U, tperp, over, run_nr=1):
    
    dname = (os.path.expanduser("~") + "/Data/" + "/Density_ref/"
             + "rho_of_mu_tperp=" + str(tperp) + ",U=" + str(U) + "_"
             + str(run_nr) + "/")
    op_poss = int(ord_par*10)/10
    eff_op = op_poss if ord_par-op_poss<0.05 else op_poss+0.1
        
    fetchname = "ord_par=" + str(eff_op)
    
    data = []
    with open(dname+fetchname+".txt", "r") as fr:
        for line in fr:
            data.append(float(line.strip("\n")))
    listlen = int(len(data)/2)
    den_list = data[:listlen]
    mu_list = data[listlen:]
    
    HF_flag = False
    if over:
        den_list.reverse()
        mu_list.reverse()
        for i in range(len(den_list)):
            if den_list[i] < 0.49:
                mug = mu_list[i]
                HF_flag = True
                break
    else:
        for i in range(len(den_list)):
            if den_list[i] > 0.51:
                mug = mu_list[i]
                HF_flag = True
                break
        
    if not HF_flag:
        print("No half-filling")
    
    return mug

def SMF_loop(tperp, g1, g2, N, chi, T, rho_maxerr=1e-4, orp_maxerr=1e-5):
    dt = 0.1
    model = "HCboson"
    order = "fourth"
    d = 2
    algo = "tDMRG"
    
    step_num = int(T / dt)
    
    M = st.Measure()
    a = np.array([[0, 0], [1, 0]])
    adag = np.array([[0, 1], [0, 0]])
    num_op = np.matmul(adag, a)
    a_exp_guess = np.sqrt(1 / 2)
    i = 0
    err = 1
    g2[1] = 4 * tperp * a_exp_guess
    ord_pars = [a_exp_guess]
    mu_list = [g2[0]]
    dens = 1 / 2
    
    while i < 100 and err > orp_maxerr:
        H = st.Hamiltonian(g1, g2, N, dt, d, chi, model, TO=order,
                           grow_chi=False)
        Psi = st.StateChain(N, d, chi, algo)
        Psi = H.time_evolve(Psi, step_num, algo, fast_run=True)
        
        new_dens = 0
        for k in range(N):
            new_dens += M.expec(Psi, num_op, k)
        new_dens /= N
        
        if abs(new_dens - dens)/dens > rho_maxerr:
            if new_dens > dens:
                over = True
            else:
                over = False
            mu_guess = guess_mu(g2[1], g1[1], tperp, over)
            g2[0] = new_mu(g1, [mu_guess, g2[1]], N, dt, d, chi, T, num_op,
                           g2[0], new_dens)
        
        mu_list.append(g2[0])
            
        new_ord_par = M.expec(Psi, a, int(N / 2))
        print(abs(new_ord_par))
        err = abs((abs(ord_pars[i]) - abs(new_ord_par)) / abs(ord_pars[i]))
        g2[1] = abs(4*new_ord_par*tperp)
        i += 1
        ord_pars.append(abs(new_ord_par))
        if i == 20 or i ==40:
            print("Error in order parameter is:", err)
            print("Truncation error:", Psi.err)
    print("Error in order parameter is:", err)
    print("Truncation error:", Psi.err)
    return ord_pars, Psi