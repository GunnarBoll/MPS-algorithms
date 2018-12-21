"""tMPS algorithm solving an antiferromagnetic transverse Ising chain
Written by: Gunnar Bollmark
"""
import numpy as np
import time as t
import scipy.sparse as sp
import datetime
import imp
import pathlib

import storage as st
import ExactDiag as ed

imp.reload(st)
imp.reload(ed)

def algo_output(Psi, H):
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
            
def exact_output(ED):
    print("\nExact diagonalization results")
    print("Exact groundstate energy:", ED.E_GS)
    print("Exact energies:",ED.elist)    

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

def filewrite(Psi, H, ED, T, keylist):
    date = str(datetime.date.today())
    run_number = 1
    direc = ("C:/Users/Gunnar/Documents/Ph.D/Learning/DMRG/Tryout code/output/"
             + date)
    while True:
        try:
            pathlib.Path(direc + "_" + str(run_number) + "/").mkdir(
                parents=True, exist_ok=False)
            break
        except FileExistsError:
            run_number += 1
    direc = direc + "_" + str(run_number) + "/"
    name = "chi=" + str(H.chi) + ",T=" + str(T) + ",dt=" + str(H.dt)
    fw = open(direc+name+".txt", "x")
    
    adag = np.array([[0, 1], [0, 0]])
    a = np.array([[0, 0], [1, 0]])
    fw.write("# Algorithm data\n")
    fw.write("# GS Energy\n")
    algo_GSE = sum(Psi.get_ener(H.Hchain))
    fw.write(str(algo_GSE) + "\n")
    meas = st.Measure()
    fw.write("# Correlators\n")
    for ind1 in range(H.N):
        for ind2 in range(H.N):
            corr = meas.correl(Psi, adag, a, ind1, ind2)
            fw.write(str(corr) + "\n")
    fw.write("# ED data\n")
    fw.write("# GS Energy\n")
    fw.write(str(ED.E_GS) + "\n")
    fw.write("# Correlators\n")
    for ind1 in range(ED.N):
        for ind2 in range(ED.N):
            corr = ED.ED_correl(ED.GS, adag, a, ind1, ind2)
            fw.write(str(corr) + "\n")
    fw.close()
    

# Main program
def main():
    model = input("Enter model (Heisen, HCboson): ")
    directory = ("C:/Users/Gunnar/Documents/Ph.D/Learning/DMRG/Tryout code/"
                 + "models/" + model + ".txt")
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
    H = st.Hamiltonian(g1, g2, N, dt, d, chi, model, TO=order)
    
    # Build the initial state
    Psi = st.StateChain(N, d, algo)
        
    # Time evolve the state (imaginary time)
    start = t.process_time()
    Psi = H.time_evolve(Psi, step_number, algo)
    end = t.process_time()
    print("Time taken for algorithm is:", end-start)
    
    # Output data to terminal
    algo_output(Psi, H)

    # Exact diagonalization trotter
    if N <= 22:
        start = t.process_time()
        ED = ed.ExactD(g1, g2, N, dt, d, model, TO=order)
        end1 = t.process_time()
        ED.exact_GS()
        end2 = t.process_time()
        # exact_output(ED)
        print("\nTime taken for full Hamiltonian:", end1-start)
        print("\nTime taken for ED is:", end2-end1)
    
    if model == "HCboson":
        Free = st.FreeFerm(g1, g2, N)
        # print("\nFree fermion result:", Free.E_GS)
    
    filewrite(Psi, H, ED, T, keylist)
    
    return 0
main()
