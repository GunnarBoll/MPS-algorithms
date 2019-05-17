# -*- coding: utf-8 -*-
"""
Created on Fri May 17 12:27:25 2019

@author: Gunnar
"""
import sys
import matplotlib.pyplot as plt

from data_retrieval import get_plot_data

def plotdat_2D(*args):
    if args == ():
        args = sys.argv[1:]
    
    xdat, ydat = get_plot_data(*args)
    
    plt.plot(xdat, ydat)
    
    return

plotdat_2D('Energy', 'U', 0.05, 100, 49, 'inf')
plotdat_2D('OrderPar', 'U', 0.05, 'inf', 0, 'inf')