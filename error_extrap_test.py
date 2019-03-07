# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 13:43:10 2019

@author: Gunnar
"""

import numpy as np
import scipy as sci
import importlib as imp
import os
import matplotlib.pyplot as plt

import storage as st

imp.reload(st)

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
    
    chi_list = [28, 30, 32, 34]
    err = []
    ener = []
    
    for chi in chi_list:
        H = st.Hamiltonian(g1, g2, N, dt, d, chi, model, TO=order, 
                           grow_chi=False)
        Psi = st.StateChain(g1, g2, N, d, chi, algo)
        Psi = H.time_evolve(Psi, step_number, algo, fast_run=False)
        
        ener += [sum(Psi.get_ener(H.Hchain))]
        err += [Psi.err]
    
    
    plt.figure(1)
    plt.plot(err, ener)
    p = sci.polyfit(err, ener, 1)
    xran = np.linspace(0, max(err), 20)
    fitfunc = lambda x: p[0]*x + p[-1]
    plt.plot(xran, fitfunc(xran))
    Free = st.FreeFerm(g1, g2, N)
    gs_list = Free.E_GS*np.ones(len(xran))
    plt.plot(xran, gs_list)
    
    plt.show()
    
    print("Algo energies:", ener)
    print("Errors:", err)
    
    
    
    print("Extrap energy", p[-1])
    
    print("\nFree fermion result:", Free.E_GS)
main()