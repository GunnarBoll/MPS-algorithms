# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 10:37:49 2019

@author: Gunnar
"""
import importlib as imp

import matplotlib.pyplot as plt

import storage as st

imp.reload(st)

def symmetry_check(Psi, H, oper, tag, plot=False):
    M = st.Measure()
    
    measures = [abs(M.expec(Psi, oper, ind)) for ind in range(Psi.N)]
    
    residues = [abs(measures[ind]-measures[Psi.N-ind-1])/abs(measures[ind])
                for ind in range(Psi.N)]
    
    if plot:
        
        fig, axes = plt.subplots(nrows=2,ncols=1, sharex="row")
        plt.subplots_adjust(hspace=0.3)
        
        axes[0].plot(range(1, Psi.N+1), measures)
        axes[0].set_ylabel("Measure " + tag)
        axes[1].plot(range(1, Psi.N+1), residues)
        axes[1].set_xlabel("Site number")
        axes[1].set_ylabel("Residue")
        
        plt.show()
        