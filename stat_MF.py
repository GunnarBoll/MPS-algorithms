"""Loop for self-consistent static mean-field calculation."""

import numpy as np

import storage as st

def main():
    
    g1 = [1., 0.]
    g2 = [1., 1.]
    dt = 0.1
    T = 50
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
    
    while i < 10:
        H = st.Hamiltonian(g1, g2, N, dt, d, chi, model, order)
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
        err = abs(g2[1] - new_alp)
        g2[1] = new_alp
        i += 1
        print(g2[1])
    
    print(err)
    print(g2[1])
    
    return
    
main()