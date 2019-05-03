# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 15:51:38 2019

@author: Gunnar
"""
import sys
import scipy.optimize as sci

from proj_storage import proj_store

def get_extr_dat(N, tperp, ind):
    name = ("/proj/snic2019-8-26/extr_orp_vs_U/tperp=" + str(tperp) + "/N=" 
            + str(N) + "/order_param")
    with open(name, 'r') as fr:
        orps = fr.read().splitlines()
        orp = eval(orps[ind])
    return orp

def powlaw_extr(xdat, ydat):
    fitfunc = lambda x, a, b, c: a*x**(b) + c
    p, cov = sci.curve_fit(fitfunc, xdat, ydat)
    return p[-1]

def main():
    
    isize_orp = []
    tperp = eval(sys.argv[1])
    
    N_list = [n*10 for n in range(2,7)] + [80]
    inv_N = [1/N for N in N_list]
    inv_N.reverse()
    
    U_list = [i/4 for i in range(21)]
    
    for ind, U in enumerate(U_list):
        orp_v_N = []
        for N in N_list:    
            orp_v_N.append(get_extr_dat(N, tperp, ind))
        orp_v_N.reverse()
        orp = powlaw_extr(inv_N, orp_v_N)
        isize_orp.append(orp)
    
    dname = "isize_orp/tperp=" + str(tperp) + "/"
    fname = "order_param"
    
    proj_store(dname, fname, isize_orp)

main()