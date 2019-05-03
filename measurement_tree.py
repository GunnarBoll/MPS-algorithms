# -*- coding: utf-8 -*-
"""
Created on Fri May  3 08:52:22 2019

@author: Gunnar
"""
import sys
import importlib as imp

import mptk_class as mp
from proj_storage import proj_store

imp.reload(mp)

def measure_tree_builder():    
    tperp = eval(sys.argv[1])
    N = eval(sys.argv[2])
    numb = eval(sys.argv[3])
    U = eval(sys.argv[4])
    chi = eval(sys.argv[5])
    obser = str(sys.argv[6])
    
    state_direc = ("mptk_states/tperp=" + str(tperp) + "/N=" + str(N) + "/n="
                   + str(numb) + "/U=" + str(U) + "/chi=" + str(chi) + "/")
    
    sol = mp.MPTKState(state_direc)
    
    if obser == "Energy":
        measure = [sol.get_ener()]
    elif obser == "OrderPar":
        measure = [sol.get_orp()]
    elif obser == "Density":
        measure = [sol.get_dens()]
    elif obser == "Trunc_err":
        measure = [sol.get_trunc_err()]
        
    save_direc = ("measurements/tperp=" + str(tperp) + "/N=" + str(N) + "/n="
                  + str(numb) + "/U=" + str(U) + "/chi=" + str(chi) + "/")
    
    proj_store(save_direc, obser, measure, replace=True)
    
measure_tree_builder()