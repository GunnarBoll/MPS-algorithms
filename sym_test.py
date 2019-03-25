# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 11:13:00 2019

@author: Gunnar
"""
import importlib as imp
import numpy as np

from GS_get import get_GS
from symmetry_func import symmetry_check
import storage as st

imp.reload(st)

def main():
    direc = "C:/Users/Gunnar/Documents/Ph.D/Data/Static_MF/"
    name = "SMF_N=20_1/GS_N=20,U=0.0.txt"
    Psi = get_GS(direc+name)
    
    g1 = Psi.g1
    g2 = Psi.g2
    dt = 0.1
    d = Psi.d
    chi_init = Psi.chi
    N = Psi.N
    order = "fourth"
    model = "HCboson"
    
    H = st.Hamiltonian(g1, g2, N, dt, d, chi_init, model, TO=order, 
                       grow_chi=False)
    
    a = np.array([[0, 0], [1, 0]])
    adag = np.array([[0, 1], [0, 0]])
    num_op = np.matmul(adag, a)
    symmetry_check(Psi, H, a, "Ord. param.", plot=True)
    symmetry_check(Psi, H, num_op, "Density", plot=True)
main()