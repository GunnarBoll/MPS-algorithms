# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 07:58:51 2019

@author: Gunnar
"""
import numpy as np

def make_plot_dat(alf, chi, T, dt, ED_EGS, ED_corr, i, j):
    direc_name = ("C:/Users/Gunnar/Documents/Ph.D/Data/Benchmarking/alf=" 
                  + str(alf) + "_1")
    file_name = ("chi=" + str(chi) + ",T=" + str(T) + ",dt="
                 + str(dt) + ".txt")
    with open(direc_name+file_name, "r") as fr:
        file_data = [float(line) for line in fr]
    E = file_data.pop(0)
    N = int(np.sqrt(len(file_data)))
    corr_mat = np.array(file_data).reshape(N, N)
    corr = corr_mat[i, j]
    
    diff_E = abs(E - ED_EGS) / abs(E)
    diff_corr = abs(corr - ED_corr) / abs(corr)
    
    return diff_E, diff_corr

def benchmark_dataread():
    
        
    
    dt_list = [0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1]
    T_list = [20, 24, 28, 30, 34, 40, 65]
    chi_max_list = [10, 12, 14, 16, 18, 20, 30, 50, 70]
    
    alp_list = [0, 0.01, 0.1, 1]
    
    dt_mp = 0.1
    chi_mp = 70
    T_mp = 100
    
    params = [[chi, T_mp, dt_mp] for chi in chi_max_list]
    params += [[chi_mp, T, dt_mp] for T in T_list]
    params += [[chi_mp, T_mp, dt] for dt in dt_list]
    
    ci1 = 2
    ci2 = 7
    
    direc = "C:/Users/Gunnar/Documents/Ph.D/Data/Benchmarking/alf=0.0_1/"
    ED_file = "bench_ED.txt"
    with open(direc + ED_file, "r") as fr:
        ed_dat = [float(line) for line in fr]
        ED_EGS = ed_dat.pop(0)
        N = int(np.sqrt(len(ed.dat)))
        ED_corr_mat = np.array(ed_dat).reshape(N, N)
        ED_corr = ED_corr_mat[ci1, ci2]
    
    for alf in alp_list:
        dt_data = [
            [make_plot_dat(alf, chi_mp, T_mp, dt, ED_EGS, ED_corr, ci1, ci2)] 
            for dt in dt_list
            ]
        chi_data = [
            [make_plot_dat(alf, chi, T_mp, dt_mp, ED_EGS, ED_corr, ci1, ci2)] 
            for chi in chi_max_list
            ]
        T_data = [
            [make_plot_dat(alf, chi_mp, T, dt_mp, ED_EGS, ED_corr, ci1, ci2)] 
            for T in T_list
            ]
        
            
benchmark_dataread()