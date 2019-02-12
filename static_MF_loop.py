import numpy as np

import storage as st

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
    a_exp_guess = 1 / np.sqrt(N / 2)
    i = 0
    err = 1
    g2[1] = 4 * tperp * a_exp_guess
    ord_pars = [a_exp_guess]
    mu_list = [g2[0]]
    
    while i < 20 and err > 10 ** -4:
        H = st.Hamiltonian(g1, g2, N, dt, d, chi, model, TO=order,
                           grow_chi=False)
        Psi = st.StateChain(N, d, chi, algo)
        Psi = H.time_evolve(Psi, step_num, algo)
        
        if i == 0:
            dens = 0
            for k in range(N):
                dens += M.expec(Psi, num_op, k)
        else:
            new_dens = 0
            for k in range(N):
                new_dens += M.expec(Psi, num_op, k)
            g2[0] = g2[0]*dens / new_dens
            dens = new_dens
            mu_list.append(g2[0])
            
        new_ord_par = M.expec(Psi, a, int(N / 2))
        err = abs((abs(ord_pars[i]) - abs(new_ord_par)) / abs(ord_pars[i]))
        g2[1] = 4*new_ord_par*tperp
        i += 1
        ord_pars.append(abs(new_ord_par))
    return ord_pars, Psi