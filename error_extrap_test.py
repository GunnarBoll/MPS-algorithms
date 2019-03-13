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
from GS_get import get_GS

imp.reload(st)

def main():
    model = "HCboson"
    
    direc = "C:/Users/Gunnar/Documents/Ph.D/Data/Static_MF/"
    file = "SMF_N=40_1/GS_N=40,U=3.0.txt"
    
    Psi = get_GS(direc+file)
    
    g1 = Psi.g1
    g2 = Psi.g2
    dt = 0.1
    d = Psi.d
    chi_init = Psi.chi
    N = Psi.N
    T = 100
    step_number = int(T / dt)
    order = "fourth"
    algo = "tDMRG"
    
    H = st.Hamiltonian(g1, g2, N, dt, d, chi_init, model, TO=order, 
                       grow_chi=False)
    
    print(Psi.err)
    
    M = st.Measure()
        
    chi_list = [Psi.chi-2, Psi.chi-4]
    err = [Psi.err]
    a = np.array([[0, 0], [1, 0]])
    measurement = [abs(M.expec(Psi, a, int(Psi.N / 2)))]
    
    for chi in chi_list:
        H.chi = chi
        H.chi_max = chi
        Psi.chi = chi
        Psi = H.time_evolve(Psi, 10, algo, fast_run=False)
        
        measurement += [abs(M.expec(Psi, a, int(Psi.N / 2)))]
        err += [Psi.err]
    
    
    plt.figure(1)
    plt.plot(err, measurement, "ro")
    p = sci.polyfit(err, measurement, 1)
    xran = np.linspace(0, max(err), 20)
    fitfunc = lambda x: p[0]*x + p[-1]
    plt.plot(xran, fitfunc(xran))
#    Free = st.FreeFerm(g1, g2, N)
#    gs_list = Free.E_GS*np.ones(len(xran))
#    plt.plot(xran, gs_list)
    
    plt.show()
    
    print("Algo measure:", measurement)
    print("Errors:", err)
    
    
    
    print("Extrap measure", p[-1])
    
#    print("\nFree fermion result:", Free.E_GS)
main()