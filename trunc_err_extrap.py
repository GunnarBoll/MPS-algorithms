# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 12:26:26 2019

@author: Gunnar
"""
import sys
import scipy as sci

from proj_storage import proj_store

def get_measure(direc, ind):
    forp = "order_param"
    ferr = "trunc_err"
    
    with open(direc+forp, 'r') as fr:
        orpvU = fr.read().splitlines()
        orp = eval(orpvU[ind])
    with open(direc+ferr, 'r') as fr:
        truncs = fr.read().splitlines()
        trunc_err = eval(truncs[ind])
    
    return orp, trunc_err

def lin_extr(xdat, ydat):
    p = sci.polyfit(xdat, ydat, deg=1)
    return p[-1]

def main():
    extr_orp = []
    
    N = int(sys.argv[1])
    tperp = float(sys.argv[2])
    chi_max = int(sys.argv[3])
    U_list = [i/4 for i in range(21)]
    
    for ind, U in enumerate(U_list):
        orpl = []
        trunc_l = []
        for chi in range(chi_max, chi_max-9, -2):
            data_direc = ("/proj/snic2019-8-26/orp_vs_U/tperp=" + str(tperp)
                          + "/N=" + str(N) + "/chi=" + str(chi) + "/")
            
#            data_direc = ("C:/Users/Gunnar/Documents/N=80/chi=" + str(chi) 
#                          + "/")
            
            orp, trunc_err = get_measure(data_direc, ind)
            orpl.append(orp)
            trunc_l.append(trunc_err)
        extr_orp.append(lin_extr(trunc_l, orpl))
    
    save_dir = ("extr_orp_vs_U/tperp=" + str(tperp)+ "/N=" + str(N) + "/")
    
    save_file = "order_param"
    
    proj_store(save_dir, save_file, extr_orp)
    
    return

main()