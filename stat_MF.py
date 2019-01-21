"""Loop for self-consistent static mean-field calculation."""

import numpy as np
import imp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import storage as st

imp.reload(st)

def exp_func(x, a, b, c):
    return a*np.exp(-b*x) + c

def main():
    
    g1 = [2., 2.]
    g2 = [2., 1]
    dt = 0.1
    T = 100
    step_num = int(T / dt)
    chi = 64
    d = 2
    N = 4
    model = "HCboson"
    order = "fourth"
    algo = "tDMRG"
    
    M = st.Measure()
    a = np.array([[0, 0], [1, 0]])
    adag = np.array([[0, 1], [0, 0]])
    num_op = np.matmul(adag, a)
    i = 0
    err = 1
    alp_list = [g2[1]]
    
    while i < 10:
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
            g2[0] = dens / new_dens
            
        new_alp = 4 * M.expec(Psi, a, int(N/2))
        err = (abs(g2[1]) - abs(new_alp))/abs(g2[1])
        g2[1] = new_alp
        i += 1
        alp_list.append(abs(g2[1]))
    
    
    plt.plot(range(1, len(alp_list) + 1), alp_list)
    popt, pcov = curve_fit(exp_func, range(1, 4), alp_list[:3],
                               p0=(-1, 1, 2))
    xx = np.linspace(1, len(alp_list))
    yy = exp_func(xx, *popt)
    plt.plot(xx,yy)
    plt.show()
    
    return
    
main()