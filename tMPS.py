"""
tMPS algorithm solving an either a Heisenberg or Hardcore boson model.
Written by: Gunnar Bollmark
"""
import numpy as np
import time as t
import scipy.sparse as sp
import importlib as imp
import os

import storage as st
import ExactDiag as ed

imp.reload(st)
imp.reload(ed)
imp.reload(np)

# Calculates and prints some interesting output if called. For small systems
# it makes a double check with a full Hilbert space groundstate.
def algo_output(Psi, H):
    adag = np.array([[0, 1], [0, 0]])
    a = np.array([[0, 0], [1, 0]])
    num_op = np.matmul(adag,a)
    print("\nMPS algorithm results")
    E = Psi.get_ener(H.Hchain)
    print("GS energy:", sum(E))
    print("Bond energies:", E)
        
    if H.N < 12:
        B = Psi.B
        psi = B[0]
        for i in range(len(B) - 1):
            psi = np.tensordot(psi, B[i + 1], (i+2, 1))
        psi = np.reshape(psi, (H.d ** H.N))
        print("Energy from product state: ", sum(ed.ExactD.get_ener(H, psi)))
        
    M = st.Measure()
    dens = 0
    for i in range(H.N):
        dens += M.expec(Psi, num_op, i)
    print("Expec algo:", dens/H.N)

# Given an Exact Diagonalization object prints some interesting quantities
# for comparison with the algorithm results.
def exact_output(ED):
    adag = np.array([[0, 1], [0, 0]])
    a = np.array([[0, 0], [1, 0]])
    num_op = np.matmul(adag,a)
    print("\nExact diagonalization results")
    print("Exact groundstate energy:", ED.E_GS)
    print("Exact energies:",ED.elist)
    dens = 0
    for i in range(ED.N):
        dens += ED.ED_correl(ED.GS, adag, a, i, i)
    print("Exact dens", dens/ED.N)
    
# Checks the variance and expectation value of an observable op. This is useful
# for finding states of definite value in some observable ( i.e. conserved
# particle number)
def obser_test(ED, op):
    corr_list = []
    maglist = []
    t1 = t.process_time()
    for k in range(len(ED.elist)):
        ED_corr = 0
        ED_mag = []
        for i in range(ED.N):
            mag = ED.ED_correl(ED.GSl[: ED.d**ED.N, k], op, sp.eye(2), i,
                               ED.N-1)
            ED_mag.append(mag)
            for j in range(ED.N):
                ED_corr += ED.ED_correl(ED.GSl[: ED.d**ED.N, k], op, op, i, j)
        ED_corr = ED_corr - sum(ED_mag)**2
        corr_list.append(ED_corr)
        maglist.append(sum(ED_mag))
    t2 = t.process_time()
    print("\nVariance of observable")
    print(corr_list)
    print("Expectation value of observable")
    print(maglist)
    print("Time taken for measurement (ED):", t2-t1)    

# Main program
def main():
    # model = input("Enter model (Heisen, HCboson): ")
    model = "HCboson"
    
    # Fetches parameters from a file in current directory with name model.txt
    directory = (os.path.dirname(os.path.realpath(__file__)) + "/" + model
                 + ".txt")
    f = open(directory, "r")
    inputlist = []
    keylist = {}
    for line in f:
        if line.startswith('#'):
            continue
        d = []
        inputlist.append(line.rstrip('\n').split(" = "))
        x = inputlist[-1][0]
        yl = inputlist[-1][1].split(',')
        for y in yl:
            try:
                d.append(int(y))
            except ValueError:
                try:
                    d.append(float(y))
                except ValueError:
                    d.append(y)
        if len(d) == 1:
            d = d[0]
        keylist[x] = d
    f.close()
    
    g1 = keylist['g1']
    g2 = keylist['g2']
    dt = keylist['dt']          #Time step
    d = keylist['d']            #One-particle Hilbert space dim
    chi = keylist['chi']        #Maximum MPS dim      
    N = keylist['N']            #Site number
    T = keylist['T']            #Total time evolved
    step_number = int(T / dt)   #Number of time steps
    order = keylist['order']    #Trotter order
    algo = keylist['algo']      #Which algorithm to use
    
    # H belongs to class Hamiltonian which features time evolution and model
    # construction.
    H = st.Hamiltonian(g1, g2, N, dt, d, chi, model, TO=order, grow_chi=False)
    
    # Build the initial state
    Psi = st.StateChain(g1, g2, N, d, chi, algo, bis_err=10**-11)
        
    # Time evolve the state (imaginary time)
    start = t.process_time()
    actstart = t.time()
    Psi = H.time_evolve(Psi, step_number, algo, fast_run=False)
    end = t.process_time()
    actend = t.time()
    print("Process time taken for algorithm is:", end-start)
    print("Actual time taken:", actend-actstart)
    
    # Output data to terminal
    algo_output(Psi, H)

    # Exact diagonalization trotter
    if N <= 20:
        start = t.process_time()
        ED = ed.ExactD(g1, g2, N, d, model, TO=order)
        end1 = t.process_time()
        ED.exact_GS()
        end2 = t.process_time()
        exact_output(ED)
        print("\nProcess time taken for full Hamiltonian:", end1-start)
        print("\nProcess time taken for ED is:", end2-end1)
    
    if model == "HCboson":
        Free = st.FreeFerm(g1, g2, N)
        print("\nFree fermion result:", Free.E_GS)
    
    # ener = 0
    # M = st.Measure()
    # for i in range(N):
    #     ener += M.correl()
    
    print(Psi.err)
    return Psi
psi2 = main()