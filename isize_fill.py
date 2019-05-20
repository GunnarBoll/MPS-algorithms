# -*- coding: utf-8 -*-
"""
Created on Mon May  6 14:37:14 2019

@author: Gunnar
"""
import sys
import scipy.optimize as sciop
import scipy as sci
from numpy import inf

from data_retrieval import get_plot_data
from proj_storage import proj_store
from cwd_storage import cwd_store

def powlaw_extr(xdat, ydat):
    fitfunc = lambda x, a, b, c: a*x**(b) + c
    p, cov = sciop.curve_fit(fitfunc, xdat, ydat)
    return p[-1]

def quad_extr(xdat, ydat):
    p = sci.polyfit(xdat, ydat, deg=2)
    return p[-1]

def obs_treatment(inv_N, obsdat, obs):
    if obs == "OrderPar":
        obser_dict = dict(zip(obsdat, inv_N))
        for obs in obsdat:
            if abs(obs) < 1e-6:
                obser_dict.pop(obs)
            elif obser_dict[obs] == 1/20 or obser_dict[obs] == 1/30:
                obser_dict.pop(obs)
        inv_N = list(obser_dict.keys())
        obsdat = list(obser_dict.values())
        print(inv_N)
        print(obsdat)
        extr_func = powlaw_extr
    elif obs == "ChargeGap":
        extr_func = quad_extr
        
    return inv_N, obsdat, extr_func

def treefill_isize(*args):
    if args == ():
        args = sys.argv[1:]
        
    obser = str(args[0])
    tperp = eval(args[1])
    n = eval(args[2])
    U = eval(args[3])
    chi = eval(args[4])
    
    N_dat1, obser_dat = get_plot_data(obser, 'N', tperp, n, U, chi)
    
    if inf in N_dat1: 
        obser_dat.pop()
        N_dat1.pop()
    
    inv_N = [1/size for size in N_dat1]
    inv_N.reverse()
    obser_dat.reverse()
    
    
    inv_N, obser_dat, extrap_function = obs_treatment(inv_N, obser_dat, obser)
    
    try:
        isize_obser = extrap_function(inv_N, obser_dat)
    except RuntimeError:
        print("Unable to fit to suggested function. U is ", U)
        isize_obser = obser_dat[0]
    
    folder = ('/measurements/tperp=' + str(tperp) + '/N=inf' + '/n='
              + str(n) + '/U=' + str(U) + '/chi=' + str(chi) + '/')
    
    proj_store(folder, obser + '.dat', isize_obser, replace=True)
    
treefill_isize()