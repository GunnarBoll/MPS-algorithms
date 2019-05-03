# -*- coding: utf-8 -*-
"""
Created on Thu May  2 11:28:19 2019

@author: Gunnar
"""

import importlib as imp
import sys

import mptk_class as mp
from proj_storage import proj_store

imp.reload(mp)

def get_dat(sol_name, obser):
    sol = mp.MPTKState(sol_name)
    if obser == "orp":
        fhandle = sol.get_orp
    elif obser == "Energy":
        fhandle = sol.get_ener
    elif obser == "Num":
        fhandle = sol.get_dens
    
    obsdat = fhandle()
    trunc_err = sol.get_trunc_err()
    return obsdat, trunc_err

def get_obser():
    observables = []
    trunc_errs = []
    
    N = eval(sys.argv[1])
    chi = eval(sys.argv[2])
    tperp = eval(sys.argv[3])
    relb = eval(sys.argv[4])
    U_start = eval(sys.argv[5])
    U_end = eval(sys.argv[6])
    num_points = eval(sys.argv[7])
    measure = str(sys.argv[8])
    U_list = [U_start + i*(U_end-U_start)/num_points
              for i in range(num_points)]
    
    num_boson = N + relb
    
    
    for U in U_list:
        loc = ("mptk_states/N=" + str(N) + ",tperp=" + str(tperp)
               + "/U=" + str(U) + "/chi=" + str(chi))
        obs, trunc_err = get_dat(loc, measure)
        observables.append(obs)
        trunc_errs.append(trunc_err)
    
    
    
    dname = ("mptk_states/N=" + str(N) + ",n=" + str(num_boson) + ",tperp="
             + str(tperp) + "/")
    
    proj_store(dname, measure + "_chi=" + str(chi), observables)    
    proj_store(dname, "trunc_err_chi=" + str(chi), trunc_errs)
    
get_obser()