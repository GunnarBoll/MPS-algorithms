# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 11:41:36 2019

@author: Gunnar
"""
import sys
import math
import os
import importlib as imp

import mptk_class as mp
from guess_mu import guess_mu
from cwd_storage import cwd_store

imp.reload(mp)

def new_mu(mu_gu, new_dens, goal_dens, mpsol, rtol):
    
    err = abs(new_dens - goal_dens) / abs(new_dens)
    
    mu_old = mpsol.mu
    alp = mpsol.alp
    N = mpsol.N
    U = mpsol.U
    chi = mpsol.chi
    cl_flag = mpsol.cluster
    
    
    if new_dens < goal_dens:
        mu0 = mu_old
        dens0 = new_dens
    else:
        mu1 = mu_old
        dens1 = new_dens
    
    ind = 0
    direc = lambda i: mpsol.loc + "/mu_" + str(i)
    
    
    while ind < 40 and err > rtol:

        musol = mp.MPTKState(direc(ind), mpsol.run_script, 
                             [N, U, mu_gu, alp, chi], model="SMF",
                             cluster=cl_flag)
        musol.mptk_run()
        
        mu_dens = get_av_expec(musol, "N", 1, N+1)
        
        dens_err = abs(mu_dens - goal_dens) / abs(mu_dens)
        
        if mu_dens < new_dens:
            mu0 = mu_gu
            dens0 = mu_dens
        else:
            mu1 = mu_gu
            dens1 = mu_dens
        
        if abs(dens_err) < rtol or ind == 39:
            mu = mu_gu
            if ind==39:
                print("Density error is:", dens_err)
            break
        else:
            yp = ( [mu0, mu1] if mu0 < mu1 else [mu1, mu0])
            xp = ([dens0, dens1 if mu0 < mu1 else [dens1, dens0]])
            
            mu_gu = man_interp(goal_dens, xp, yp)
            
        
        ind += 1
    
    return mu

def man_interp(x, xdat, ydat):
    k = (ydat[1] - ydat[0]) / (xdat[1] - xdat[0])
    m = ydat[0] - k*xdat[0]
    assert abs(m - (ydat[1] - k*xdat[1]))/abs(m) < 1e-6
    interp_val = k*x + m
    return interp_val

def get_av_expec(mpsol, oper, start, end):    
    avexpec = 0
    av_num_sites = 0
    op_func = lambda ind: oper + "(" + str(ind) + ")"
        
    for orpind in range(start, end):
        av_num_sites += 1
        avexpec += mpsol.expec(op_func(orpind))
    
    avexpec /= av_num_sites
    
    return avexpec

def mptk_SMF():
    rho_max_err = 1e-5
    orp_max_err = 1e-6
    
    N = eval(sys.argv[1])
    U = eval(sys.argv[2])
    tperp = eval(sys.argv[3])
    chi = eval(sys.argv[4])
    num_bos = eval(sys.argv[5])
    cl_flag = bool(int(sys.argv[6]))
    
    mu = 0
    orp_guess = 1 / math.sqrt(2)
    alp = 4 * tperp * orp_guess
    dens = num_bos / N
    
    mu_list = [0]
    orp_list = [orp_guess]
    err = 1
    i = 0
    
    run_script = "bin/mptk_script.sh"
    dname = lambda p1, p2, p3, p4, p5: ("mptk_states/tperp=" + str(p1) + "/N="
                                        + str(p2) + "/n=" + str(p3) + "/U="
                                        + str(p4) + "/chi=" + str(p5))
    
    while True:
        arg_list = [N, U, mu, alp, chi]
        
        mpsol = mp.MPTKState(dname(tperp, N, num_bos, U, chi), run_script,
                             arg_list, model="SMF", cluster=cl_flag)
        
        is_finished = mpsol.mptk_run()
        if is_finished:
            return
        
        new_dens = get_av_expec(mpsol, "N", 1, N+1)
        density_err = abs(new_dens - dens)/dens
        
        if density_err > rho_max_err:
            if new_dens > dens:
                over = True
            else:
                over = False
            mu_guess = guess_mu(alp, U, tperp, over, new_dens, goal_dens=dens)
            mu = new_mu(mu_guess, new_dens, dens, mpsol, rho_max_err)
        
        
        av_start = int(N/2)-int(N/4) + 1
        av_end = int(N/2) + int(N/4) + 1 if N%2==0 else int(N/2) + int(N/4) + 2
                
        new_orp = abs(get_av_expec(mpsol, "B", av_start, av_end))
        err = abs(new_orp - orp_list[-1]) / abs(new_orp)
        print(new_orp)
        
        orp_list.append(new_orp)
        mu_list.append(mu)
        alp = 4 * tperp * new_orp
        
        if i%20 == 0:
            print("Error in order param. is", err)
        
        if ( (err < orp_max_err and density_err < rho_max_err) or i == 500
            or orp_list[-1] < 1e-8 ):
            break
        else:
            mpsol.delete_solution()
                
        i += 1
    print("Final error in order param. is", err)
    
    if cl_flag and not is_finished:
#        Do NOT add chi at the end: it will create a chi=$chi/chi=$chi 
#        structure
        mpsol.copy_solution("/mptk_states/tperp=" + str(tperp) + "/N=" + str(N)
                            + "/n=" + str(num_bos) + "/U=" + str(U) + "/")
    
    fol = "transf_n=" + str(num_bos) + "/SMF_N="+str(N)
    fnam = "N=" + str(N) + ",U=" + str(U) + ".txt"
    cwd_store(fol, fnam, orp_list)
    cwd_store(fol, "mu_" + fnam, mu_list)
    
    return orp_list, mu_list
mptk_SMF()