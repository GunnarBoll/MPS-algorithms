# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 13:15:50 2019

@author: Gunnar
"""
import sys
import importlib as imp
import datetime
import numpy as np

import storage as st
import ExactDiag as ed
from cwd_storage import cwd_store

imp.reload(st)
imp.reload(ed)

def main():
    date = str(datetime.date.today())
    N = int(sys.argv[1])
    dt = float(sys.argv[2])
    T = int(sys.argv[3])
    chi = int(sys.argv[4])
    alf = float(sys.argv[5])
    step_num = int(T / dt)
    g1 = [1., 1.]
    g2 = [1., alf]
    d = 2
    algo = "tDMRG"
    
    H = st.Hamiltonian(g1, g2, N, dt, d, chi, model="HCboson", TO="fourth",
                       grow_chi=False)
    Psi = st.StateChain(g1, g2, N, d, chi, algo)
    
    Psi = H.time_evolve(Psi, step_num, algo)
    
    M = st.Measure()
    
    adag = np.array([[0., 1.], [0, 0]])
    a = np.array([[0, 0], [1, 0]])
    corr_matrix = M.corr_mat(Psi, adag, a)
    E_GS = sum(Psi.get_ener(H.Hchain))
    
    data = [E_GS]
    corr_matrix = np.reshape(corr_matrix, (corr_matrix.shape[0] ** 2))
    data += [corr for corr in corr_matrix]
    
    direc_name = "alf=" + str(alf) + date
    file_name = "chi=" + str(chi) + ",T=" + str(T) + ",dt=" + str(dt) + ".txt"
    
    cwd_store(direc_name, file_name, data)
    
main()