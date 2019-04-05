# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 12:26:26 2019

@author: Gunnar
"""
import importlib as imp
import sys
import scipy as sci

import mptk_class as mp
from avorp import get_orp

imp.reload(mp)

def get_dat(sol_name, op, N):
    
    sol = mp.MPTKState(sol_name)
    orp = get_orp(op, sol, N)
    trunc_err = sol.get_trunc_err()
    return orp, trunc_err

def lin_extr(xdat, ydat):
    p = sci.polyfit(xdat, ydat, deg=1)
    return p[-1]

def save_data():
    pass

def main():
    extr_orp = []
    
    U = float(sys.argv[1])
    tperp = float(sys.argv[2])
    chi_max = int(sys.argv[3])
    
    for N in [80, 60, 50, 40, 30, 20]:
        orpl = []
        trunc_l = []
        for chi in range(chi_max, chi_max-9, -2):
            solname = ("mptk_states/N=" + str(N) + ",tperp=" + str(tperp)
                       + "/U=" + str(U) + "/chi=" + str(chi))
            orp, trunc_err = get_dat(solname, "B", N)
            orpl.append(orp)
            trunc_l.append(trunc_err)
        extr_orp.append(lin_extr(trunc_l, orpl))
    
    save_data()
    
    return