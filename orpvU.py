# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 11:48:14 2019

@author: Gunnar
"""
import importlib as imp
import sys

import mptk_class as mp
from avorp import get_orp
from proj_storage import proj_store

imp.reload(mp)

def get_dat(sol_name, op, N):
    sol = mp.MPTKState(sol_name)
    orp = get_orp(op, sol, N)
    trunc_err = sol.get_trunc_err()
    return orp, trunc_err

def get_orp_vs_U():
    orpl = []
    truncl = []
    
    N = eval(sys.argv[1])
    chi = eval(sys.argv[2])
    tperp = eval(sys.argv[3])
    U_list = [i/4 for i in range(21)]
    
    for U in U_list:
        loc = ("mptk_states/N=" + str(N) + ",tperp=" + str(tperp)
               + "/U=" + str(U) + "/chi=" + str(chi))
        orp, trunc_err = get_dat(loc, "B", N)
        orpl.append(orp)
        truncl.append(trunc_err)
    
    dname = ("orp_vs_U/tperp=" + str(tperp) + "/N=" + str(N) + "/chi=" 
             + str(chi) + "/")
    proj_store(dname, "order_param", orpl)    
    proj_store(dname, "trunc_err", truncl)
    
get_orp_vs_U()