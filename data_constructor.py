# -*- coding: utf-8 -*-
"""
Created on Fri May  3 13:59:06 2019

@author: Gunnar
"""
import sys

def data_constructor():
    data_loc = "/proj/snic2019-8-26/measurements/"
    
    params = ["tperp", "N", "n", "U", "chi"]
    
    xdat = str(sys.argv[1])
    ydat = str(sys.argv[2])
    params.remove(xdat)
    params.remove(ydat)