"""Loop for self-consistent static mean-field calculation."""

import numpy as np
import imp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from benchreader import figcloser

import storage as st

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
    
    g1 = [1., 0.]
    g2 = [0., 0.01]
    tperp = 0.05
    dt = 0.1
    T = 100
    step_num = int(T / dt)
    chi = 70
    d = 2
    N = 20
    model = "HCboson"
    order = "fourth"
    algo = "tDMRG"
    
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
    
    while i < 10 and err > 10 ** -4:
        H = st.Hamiltonian(g1, g2, N, dt, d, chi, model, TO=order,
                           grow_chi=False)
        Psi = st.StateChain(N, d, algo)
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
    
    
    figcloser()
    plt.figure(1)
    plt.plot(range(1, len(ord_pars) + 1), ord_pars)
    # popt, pcov = curve_fit(exp_func, range(1, 4), ord_pars[0 : 3],
    #                             p0=(-1, 1, 2))
    # guess = 4 * tperp * popt[2]
    # ex_res = extrap_res(g1, [mu_list[-1], guess], N, dt, d, chi, model, order,
    #                     algo, step_num, a)
    # print(abs(ex_res), ord_pars[-1], ord_pars[3])
    # 
    # xx = np.linspace(1, len(ord_pars))
    # yy = exp_func(xx, *popt)
    # plt.plot(xx,yy)
    
    plt.figure(2)
    plt.plot(range(1, len(mu_list) + 1), mu_list)
    plt.show()
    
    return
    
main()