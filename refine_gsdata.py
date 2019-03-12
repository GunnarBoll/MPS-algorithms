# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 08:48:54 2019

@author: Gunnar
Takes old Ground state data files and updates them to for longer time 
evolution.
"""
import importlib as imp
import sys
import os

from GS_get import get_GS
from GS_store import store_GS
from cwd_storage import cwd_store
import storage as st

imp.reload(st)

def refine():
    model = "HCboson"
    
    N = int(sys.argv[1])
    U = float(sys.argv[2])
    
    dt = 0.1
    hom = os.path.expanduser("~")
    dname = hom + "/Data/Static_MF/SMF_N=" + str(N) + "_1/"
    fname = "GS_N=" + str(N) + ",U=" + str(U) + ".txt"
    
    Psi = get_GS(dname + fname)
    
    H = st.Hamiltonian(Psi.g1, Psi.g2, N, dt, Psi.d, Psi.chi, model,
                       TO="fourth")
    
    Psi = H.time_evolve(Psi, 1000, "tDMRG", fast_run=False)
    
    gs_mat = store_GS(Psi)
    
    cwd_folder = "refined_gsdata"
    
    cwd_store(cwd_folder, fname, gs_mat)
    