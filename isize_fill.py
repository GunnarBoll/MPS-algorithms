# -*- coding: utf-8 -*-
"""
Created on Mon May  6 14:37:14 2019

@author: Gunnar
"""
import sys

def treefill_isize(*args):
    if args == ():
        args = sys.argv[1:]
        
    obser = str(args[0])
    tperp = eval(args[1])
    n = eval(args[2])
    U = eval(args[3])
    try:
        chi = eval(args[4])
    except NameError:
        chi = str(args[4])
    