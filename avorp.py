# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 15:40:39 2019

@author: Gunnar
"""

def get_orp(oper, mpsol, N):
    start = int(N/2)-int(N/4) + 1
    end = int(N/2) + int(N/4) + 1 if N%2==0 else int(N/2) + int(N/4) + 2
    
    orp = 0
    av_num_sites = 0
        
    for orpind in range(start, end):
        av_num_sites += 1
        orp += abs(mpsol.expec(oper, orpind))
    
    orp /= av_num_sites
    
    return orp