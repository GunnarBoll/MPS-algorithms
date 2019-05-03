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
    extr_meas = []
    
    N = eval(sys.argv[1])
    chi_max = eval(sys.argv[2])
    tperp = eval(sys.argv[3])
    relb = eval(sys.argv[4])
    U_start = eval(sys.argv[5])
    U_end = eval(sys.argv[6])
    num_points = eval(sys.argv[7])
    measure = str(sys.argv[8])
    
    U_list = [U_start + i*(U_end-U_start)/num_points
              for i in range(num_points)]
    
    num_boson = N + relb
    
    for ind, U in enumerate(U_list):
        meas_list = []
        trunc_list = []
        for chi in range(chi_max, chi_max-9, -2):
            data_direc = ("/proj/snic2019-8-26/mptk_states/N=" + str(N)
                          + ",n=" + str(num_boson) + "/tperp=" + str(tperp)
                          + "/U=" + str(U) + "/chi=" + str(chi) + "/")
            
#            data_direc = ("C:/Users/Gunnar/Documents/N=80/chi=" + str(chi) 
#                          + "/")
            
            meas, trunc_err = get_measure(data_direc, ind)
            meas_list.append(meas)
            trunc_list.append(trunc_err)
        extr_meas.append(lin_extr(trunc_list, meas_list))
    
    
    save_dir = ("/proj/snic2019-8-26/mptk_states/N=" + str(N)
                + "/n=" + str(num_boson) + "/tperp=" + str(tperp)
                + "/U=" + str(U) + "/")
    
    save_dir = ("/proj/snic2019-8-26/measurements/tperp=" + str(tperp) 
                + "/N=" + str(N) + "/")
    
    save_file = measure + "_vs_U" + ".dat"
    
    proj_store(save_dir, save_file, extr_meas)
    
    return

main()