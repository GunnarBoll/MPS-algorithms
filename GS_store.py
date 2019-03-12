# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 09:04:47 2019

@author: Gunnar
"""


def store_GS(Psi):
    GS_mat = ([coup for coup in Psi.g1] + [coup2 for coup2 in Psi.g2]
              + [Psi.N, Psi.d, Psi.err, Psi.notation, Psi.chi])
    for tens in Psi.B:
        GS_mat += [elem for elem in 
                   tens.reshape(tens.shape[0]*tens.shape[1]*tens.shape[2])]
    
    for lams in Psi.L:
        GS_mat += [lam_elem for lam_elem in lams]
    
    return GS_mat