# -*- coding: utf-8 -*-
"""
Created on Fri May 17 16:35:42 2019

@author: Gunnar
"""
import sys
from numpy import inf

from data_retrieval import get_plot_data
from proj_storage import proj_store

def charge_gap_fill(*args):
    
    if args == ():
        args = sys.argv[1:]
    
    tperp = eval(args[0])
    N = eval(args[1])
    U = eval(args[2])
    chi = eval(args[3])
    
    print(args)
    energies, fillings = get_plot_data('Energy', 'n', *args)
    
    gap = energies[0] + energies[2] - 2*energies[1]
    
    folder = ('/measurements/tperp=' + str(tperp) + '/N=' + str(N) + '/n='
              + str(int(N/2)) + '/U=' + str(U) + '/chi=' + str(chi) + '/')
    
    proj_store(folder, 'ChargeGap.dat', [gap], replace=True)

charge_gap_fill()