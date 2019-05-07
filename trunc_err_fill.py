# -*- coding: utf-8 -*-
"""
Created on Mon May  6 13:55:03 2019

@author: Gunnar
"""
import sys
import scipy as sci
from numpy import inf

from data_retrieval import get_plot_data
from proj_storage import proj_store

def lin_extrap(xdat, ydat):
    p = sci.polyfit(xdat, ydat, deg=1)
    return p[-1]

def treefill_trunc_extrap(*args, **kwargs):
    if args == ():
        args = sys.argv[1:]
    
    obser = str(args[0])
    tperp = eval(args[1])
    N = eval(args[2])
    n = eval(args[3])
    U = eval(args[4])
    
    chidat1, trunc_err = get_plot_data("Trunc_err", "chi", tperp, N, n, U)
    chidat2, obser_dat = get_plot_data(obser, "chi", tperp, N, n, U)
    
    if inf in chidat1: trunc_err.pop()
    if inf in chidat2: obser_dat.pop()
    
    trunc_err.reverse()
    obser_dat.reverse()
    
    extrap_val = lin_extrap(trunc_err, obser_dat)
    
    folder = ('/measurements/tperp=' + str(tperp) + '/N=' + str(N) + '/n='
              + str(n) + '/U=' + str(U) + "/chi=inf/")
    
    proj_store(folder, obser + '.dat', extrap_val, replace=True)
    proj_store(folder, 'Trunc_err.dat', 0, replace=True)
treefill_trunc_extrap()