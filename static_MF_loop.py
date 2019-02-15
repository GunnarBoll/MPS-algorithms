import numpy as np
import importlib as imp

import storage as st
import ExactDiag as ed

imp.reload(st)
imp.reload(ed)

def new_mu(coup1, coup2, N, dt, d, chi, T, num_op, old_mu, start_dens):
    goal_dens = N / 2
    Meas = st.Measure()
    order = "fourth"
    model = "HCboson"
    mu0 = old_mu
    dens0 = start_dens
    
    for ind in range(5):
        mu_ham = st.Hamiltonian(coup1, coup2, N, dt, d, chi, model, TO=order,
                                grow_chi=False)
        mu_psi = st.StateChain(N, d, chi, "tDMRG", fast_run=True)
        
        steps = int(T / dt)
        #for ind2 in range(5):
        mu_psi = mu_ham.time_evolve(mu_psi, steps, "tDMRG")
        
        mu_dens = 0
        for m in range(N):
            mu_dens += Meas.expec(mu_psi, num_op, m)
        mu_dens /= N
        
        dens_err = (mu_dens - goal_dens) / goal_dens
        
        print(mu_dens)
        mu1 = coup2[0]
        dens1 = mu_dens
        
        if abs(dens_err) < 10**-3 or ind == 4:
            mu = coup2[0]
            break
        else:
            den_slope = (dens1-dens0) / (mu1-mu0)
            insec = (dens0*mu1 - dens1*mu0) / (mu1-mu0)
            # coup2[0] = np.interp(goal_dens, [un_dens, ov_dens],
            #                      [un_mu, ov_mu])
            
            
            coup2[0] = (goal_dens - insec) / den_slope
        
    return mu

def SMF_loop(tperp, g1, g2, N, chi, T):
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
    
    while i < 20 and err > 10 ** -4:
        H = st.Hamiltonian(g1, g2, N, dt, d, chi, model, TO=order,
                           grow_chi=False)
        Psi = st.StateChain(N, d, chi, algo)
        Psi = H.time_evolve(Psi, step_num, algo, fast_run=True)
        
        new_dens = 0
        for k in range(N):
            new_dens += M.expec(Psi, num_op, k)
        new_dens /= N
        
        print("Start dens", new_dens)
        print("Mu", g2[0])
        if new_dens > dens:
            mu_guess = 0.9*g2[0]
        else:
            mu_guess = 1.1*g2[0]
        
        if abs(new_dens - dens)/dens > 10**-3:
            g2[0] = new_mu(g1, [mu_guess, g2[1]], N, dt, d, chi, T, num_op,
                           g2[0], new_dens)
        
        mu_list.append(g2[0])
            
        new_ord_par = M.expec(Psi, a, int(N / 2))
        err = abs((abs(ord_pars[i]) - abs(new_ord_par)) / abs(ord_pars[i]))
        g2[1] = 4*new_ord_par*tperp
        i += 1
        ord_pars.append(abs(new_ord_par))
    return ord_pars, Psi