# -*- coding: utf-8 -*-
"""
Created on Fri May  3 08:52:22 2019

@author: Gunnar
"""
import sys
import numpy as np
import importlib as imp

import mptk_class as mp
from proj_storage import proj_store

imp.reload(mp)

def measure_tree_builder():
    proj_direc = "/proj/snic2019-8-26/mptk_states/"
    
    tperp = eval(sys.argv[1])
    N = eval(sys.argv[2])
    numb = eval(sys.argv[3])
    U = eval(sys.argv[4])
    chi = eval(sys.argv[5])
    obser = str(sys.argv[6])
    
    state_direc = (proj_direc + "tperp=" + str(tperp) + "/N=" + str(N) + "/n="
                   + str(numb) + "/U=" + str(U) + "/chi=" + str(chi) + "/")
    
    sol = mp.MPTKState(state_direc)
    
    if obser == "Energy":
        handle = sol.get_ener
    elif obser == "OrderPar":
        handle = sol.get_orp
    elif obser == "Density":
        handle = sol.get_dens
    
    measure = handle()
    
    save_direc = ("measurements/tperp=" + str(tperp) + "/N=" + str(N) + "/n="
                  + str(numb) + "/U=" + str(U) + "/chi=" + str(chi) + "/")
    
    proj_store(save_direc, obser, measure, replace=True)
    
measure_tree_builder()