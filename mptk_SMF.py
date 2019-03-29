# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 11:41:36 2019

@author: Gunnar
"""
import sys
import math
import subprocess as subp

from guess_mu import guess_mu

def new_mu(U, mu_gu, alp, mu_old, new_dens, bargs, rtol):
    goal_dens = 1 / 2
    err = abs(new_dens - goal_dens) / abs(new_dens)
    
    if new_dens < goal_dens:
        mu0 = mu_old
        dens0 = new_dens
    else:
        mu1 = mu_old
        dens1 = new_dens
    
    ind = 0
    
    
    while ind < 40 and err > rtol:
        args_mu = [bargs[0] + "/mu_" + str(ind), bargs[1], str(U), str(mu_gu),
                   str(alp)]
        out = bash_call("toolkit_runs/mptk_script.sh", args_mu)
        
        mu_dens = get_av_expec(args_mu[0], "N", 1, int(args_mu[1])+1)
        
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

def bash_call(scr_name, argv):
    arg_list = ["/home/gunbo249/" + scr_name] + [str(arg) for arg in argv]
    
    res = subp.check_output(arg_list)
    
    return res

def get_av_expec(dname, oper, start, end):    
    expec = 0
    av_num_sites = 0
    
    scr = "bin/mp-expectation"
    arg0 = dname + "/GS_file.psi.12"
    arg1 = lambda ind: dname + "/lattice:" + str(oper) + "(" + str(ind) + ")"
        
    for orpind in range(start, end):
        av_num_sites += 1
        res = bash_call(scr, [arg0, arg1(orpind)])
        expec += eval(res)[0]
    
    expec /= av_num_sites
    
    return expec

def mptk_SMF():
    rho_max_err = 1e-5
    orp_max_err = 1e-6
    
    dname = sys.argv[1]
    N = int(sys.argv[2])
    U = float(sys.argv[3])
    tperp = float(sys.argv[4])
    mu = 0
    orp_guess = 1 / math.sqrt(2)
    alp = 4 * tperp * orp_guess
    dens = 1 / 2
    
    mu_list = [0]
    orp_list = [orp_guess]
    err = 1
    i = 0
    
    direc = lambda ind: dname + "_" + str(ind)
    
    
    while err > orp_max_err and i < 150:
        arg_list = [direc(i), str(N), str(U), str(mu), str(alp)]
        
        stdo = bash_call("toolkit_runs/mptk_script.sh", arg_list)
        
        new_dens = get_av_expec(arg_list[0], "N", 1, N+1)
        
        if abs(new_dens - dens)/dens > rho_max_err:
            if new_dens > dens:
                over = True
            else:
                over = False
            mu_guess = guess_mu(alp, U, tperp, over, new_dens)
            mu = new_mu(U, mu_guess, alp, mu, new_dens, arg_list, rho_max_err)
        
        
        av_start = int(N/2)-int(N/4) + 1
        av_end = int(N/2) + int(N/4) + 1 if N%2==0 else int(N/2) + int(N/4) + 2
                
        new_orp = abs(get_av_expec(arg_list[0], "B", av_start, av_end))
        err = abs(new_orp - orp_list[-1]) / abs(new_orp)
        print(new_orp)
        
        orp_list.append(new_orp)
        mu_list.append(mu)
        alp = 4 * tperp * new_orp
        
        if i%20 == 0:
            print("Error in order param. is", err)
        
        i += 1
    print("Final error in order param. is", err)
    
    return orp_list, mu_list
mptk_SMF()