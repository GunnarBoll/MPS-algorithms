# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 13:57:34 2019

@author: Gunnar
"""
import sys
import importlib as imp
import numpy as np
import datetime

import ExactDiag as ed
from cwd_storage import cwd_store

imp.reload(ed)

def main():
    date = str(datetime.date.today())
    N = int(sys.argv[1])
    alf = float(sys.argv[2])
    g1 = [1., 1.]
    g2 = [1., alf]
    model = "HCboson"
    d = 2
    
    ED = ed.ExactD(g1, g2, N, d, model)
    
    ED.exact_GS()
    
    adag = np.array([[0., 1.], [0, 0]])
    a = np.array([[0, 0], [1, 0]])
    corr_mat = [ED.ED_correl(ED.GS, adag, a, i, j) for i in range(N) for j in
                range(N)]
    
    data = [ED.E_GS] + corr_mat
    
    file_name = "bench_ED.txt"
    direc_name = "alf=" + str(alf) + date
    cwd_store(direc_name, file_name, data)

main()