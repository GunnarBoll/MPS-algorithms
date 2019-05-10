# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 08:42:44 2019

@author: Gunnar
"""
import os

def guess_mu(ord_par, U, tperp, over, dens, goal_dens=0.5, run_nr=1):
    
    dname = (os.path.expanduser("~") + "/Data/" + "/Density_ref/"
             + "rho_of_mu_tperp=" + str(tperp) + ",U=" + str(U) + "_"
             + str(run_nr) + "/")
    op_poss = int(ord_par*10)/10
    eff_op = op_poss if ord_par-op_poss<0.05 else op_poss+0.1
        
    fetchname = "ord_par=" + str(eff_op)
    
    data = []
    with open(dname+fetchname+".txt", "r") as fr:
        for line in fr:
            data.append(float(line.strip("\n")))
    listlen = int(len(data)/2)
    den_list = data[:listlen]
    mu_list = data[listlen:]
    
    HF_flag = False
    if over:
        den_list.reverse()
        mu_list.reverse()
        for i in range(len(den_list)):
            if den_list[i] < goal_dens-0.02:
                mug = mu_list[i]
                HF_flag = True
                break
    else:
        for i in range(len(den_list)):
            if den_list[i] > goal_dens+0.02:
                mug = mu_list[i]
                HF_flag = True
                break
        
    if not HF_flag:
        print("Density is:", dens)
        if over:
            mug = -10.0
        else:
            mug = 20.0
    
    return mug