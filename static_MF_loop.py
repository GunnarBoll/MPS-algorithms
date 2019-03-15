"""
Self-consistency loop for the order parameter. Attempts to keep the density at
half-filling such that superfluidity arises.
"""
import numpy as np
import importlib as imp
import os

import storage as st
import ExactDiag as ed

imp.reload(st)
imp.reload(ed)

# Function which tries to find the correct chemical potential for half-filling
# to a precision rho_maxerr via interpolation. This function requires one
# chemical potential correspondent to a density below goal_dens and one above.
def new_mu(coup1, coup2, N, dt, d, chi, T, num_op, old_mu, start_dens,
           rho_maxerr):
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
    for ind in range(40):
        mu_ham = st.Hamiltonian(coup1, coup2, N, dt, d, chi, model, TO=order,
                                grow_chi=False)
        mu_psi = st.StateChain(coup1, coup2, N, d, chi, "tDMRG")
        
        steps = int(T / dt)
        #for ind2 in range(5):
        mu_psi = mu_ham.time_evolve(mu_psi, steps, "tDMRG", fast_run=True)
        
        mu_dens = 0
        for m in range(N):
            mu_dens += Meas.expec(mu_psi, num_op, m)
        mu_dens /= N
        
        dens_err = abs((mu_dens - goal_dens)) / mu_dens
        
        if mu_dens > goal_dens:
            mu1 = coup2[0]
            dens1 = mu_dens
        else:
            mu0 = coup2[0]
            dens0 = mu_dens
        
        if abs(dens_err) < rho_maxerr or ind == 39:
            mu = coup2[0]
            if ind == 39:
                print("Density error:", dens_err)
            break
        else:
            yp = ([mu0, mu1] if mu0 < mu1 else [mu1, mu0])
            xp = ([dens0, dens1] if mu0 < mu1 else [dens1, dens0])

            coup2[0] = np.interp(goal_dens, xp, yp)
                    
    return mu

# Function which provides a guess at a chemical potential close and on the 
# opposite side of the chemical potential given by the main loop. This is done
# with reference data produced by get_chemdens.py
def guess_mu(ord_par, U, tperp, over, dens, run_nr=1):
    
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
        print("Density is:", dens)
        if over:
            mug = -10.0
        else:
            mug = 10.0
    
    return mug

# Main self-consistency loop. Given maximum allowed errors in density and
# order parameters calculates an order parameter by looping the algorithm until
# input order parameter and output order parameter agree.
def SMF_loop(tperp, g1, g2, N, chi, T, rho_maxerr=1e-4, orp_maxerr=1e-6):
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
    
    # Loops for at most 100 iterations. Some runs do not converge even at this
    # point.
    while i < 150 and err > orp_maxerr and ord_pars[-1] > 1e-8:
        H = st.Hamiltonian(g1, g2, N, dt, d, chi, model, TO=order,
                           grow_chi=False)
        Psi = st.StateChain(g1, g2, N, d, chi, algo)
        Psi = H.time_evolve(Psi, step_num, algo, fast_run=True)
        
        new_dens = 0
        for k in range(N):
            new_dens += M.expec(Psi, num_op, k)
        new_dens /= N
        
        # Only finds a new chemical potential if density deviates too much from
        # the goal density. Finds a guess which should lie on the other side
        # of goal density.
        if abs(new_dens - dens)/dens > rho_maxerr:
            if new_dens > dens:
                over = True
            else:
                over = False
            mu_guess = guess_mu(g2[1], g1[1], tperp, over, new_dens)
            g2[0] = new_mu(g1, [mu_guess, g2[1]], N, dt, d, chi, T, num_op,
                           g2[0], new_dens, rho_maxerr)
        
        mu_list.append(g2[0])
            
        new_ord_par = M.expec(Psi, a, int(N / 2))
        print(abs(new_ord_par))
        err = abs((abs(ord_pars[i]) - abs(new_ord_par)) / abs(ord_pars[i]))
        g2[1] = abs(4*new_ord_par*tperp)
        i += 1
        ord_pars.append(abs(new_ord_par))
        
        # Some print statements to see if the convergence is coming along
        if i%20 == 0:
            print("Error in order parameter is:", err)
            print("Truncation error:", Psi.err)
    print("Error in order parameter is:", err)
    print("Truncation error:", Psi.err)
    
    
#    Psi = H.time_evolve(Psi, 1000, algo, fast_run=False)
    
#    finaldens = 0
#    for i in range(N):
#        finaldens += M.expec(Psi, num_op, i)
#    finaldens /= N
#    final_err = 1e-8
#    if abs(finaldens - dens)/dens > final_err:
#        if finaldens > dens:
#            over = True
#        else:
#            over = False
#        mu_guess = guess_mu(g2[1], g1[1], tperp, over)
#        g2[0] = new_mu(g1, [mu_guess, g2[1]], N, dt, d, chi, T, num_op, g2[0],
#                       finaldens, final_err)
#        mu_list.append(g2[0])
        
#    print(abs(M.expec(Psi, a, int(N / 2))))
    return ord_pars, Psi, mu_list