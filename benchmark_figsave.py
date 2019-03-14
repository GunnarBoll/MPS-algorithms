# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 07:58:51 2019

@author: Gunnar
"""
import numpy as np
import matplotlib.pyplot as plt
import importlib as imp
import os

import ExactDiag as ed

imp.reload(ed)

def make_plot_dat(alf, chi, T, dt, ED_EGS, ED_corr, i, j):
    direc_name = ("C:/Users/Gunnar/Documents/Ph.D/Data/Benchmarking/alf=" 
                  + str(alf) + "_1/")
    file_name = ("chi=" + str(chi) + ",T=" + str(T) + ",dt="
                 + str(dt) + ".txt")
    with open(direc_name+file_name, "r") as fr:
        file_data = [float(line) for line in fr]
    E = file_data.pop(0)
    err = file_data.pop(0)
    N = int(np.sqrt(len(file_data)))
    corr_mat = np.array(file_data).reshape(N, N)
    corr = corr_mat[i, j]

    diff_E = abs(E - ED_EGS) / abs(E)
    diff_corr = abs(corr - ED_corr) / abs(corr)
    
    return [diff_E, diff_corr]

def repack(data):
    rep_dat1 = [data[i][0] for i in range(len(data))] 
    rep_dat2 = [data[j][1] for j in range(len(data))]
    return rep_dat1, rep_dat2

def make_fig(data, xdat, name, label):
    corr_data = [data[i][1] for i in range(len(data))]
    E_data = [data[i][0] for i in range(len(data))]
    
    form = "eps"
    
    Efig = plt.figure()
    plt.loglog(xdat, E_data, xdat, E_data, "ro")
    plt.xlabel(label)
    plt.ylabel("Ground state energy, Relative error")
    plt.savefig(name + "_E." + form, format=form)
    corrfig = plt.figure()
    plt.loglog(xdat, corr_data, xdat, corr_data, "ro")
    plt.xlabel(label)
    plt.ylabel("Correlator, Relative error")
    plt.savefig(name + "_Corr." + form, format=form)
    
    return Efig, corrfig

def make_subplot(alf_data, xdat, xinfo, alf_list, fignam):
    n_alf = len(alf_data)
    alf_names = [r"$\alpha$=" + str(alf) for alf in alf_list]
    fig, axes = plt.subplots(nrows=n_alf, ncols=1, sharex="all")
    plt.subplots_adjust(hspace = 0.3)
    
    form = "eps"
    
    for i in range(n_alf):
        data = alf_data[i]
        
        axes[i].set_title(str(alf_names[i]))
        axes[i].loglog(xdat, data, xdat, data, "ro")

        axes[i].set_ylabel("Rel. error")
    
    axes[i].set_xlabel(xinfo)
    fig.savefig(fignam + "." + str(form), format=form)

                
def benchmark_dataread():
    
    dt_list = [0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1]
    T_list = [20, 24, 28, 30, 34, 40, 65, 100]
    chi_max_list = [10, 12, 14, 16, 18, 20, 30, 50, 70]
    
    # alpha = 1.0 behaves very strangely
    alp_list = [0.0, 0.01, 0.1]
    
    dt_mp = 0.1
    chi_mp = 70
    T_mp = 100
    
    ci1 = 2
    ci2 = 7
    
    fig_direc = os.getcwd() + "/figs/"
    dt_name = fig_direc + "dt"
    chi_name = fig_direc + "chi"
    T_name = fig_direc + "T"
    
    dt_alf_E = []
    dt_alf_C = []
    chi_alf_E = []
    chi_alf_C = []
    T_alf_E = []
    T_alf_C = []
    
    for alf in alp_list:
        direc = ("C:/Users/Gunnar/Documents/Ph.D/Data/Benchmarking/alf="
                 + str(alf) + "_1/")
        ED_file = "bench_ED.txt"
        with open(direc + ED_file, "r") as fr:
            ed_dat = [float(line) for line in fr]
            ED_EGS = ed_dat.pop(0)
            N = int(np.sqrt(len(ed_dat)))
            ED_corr_mat = np.array(ed_dat).reshape(N, N)
            ED_corr = ED_corr_mat[ci1, ci2]
        
        
        
        dt_data = [
            make_plot_dat(alf, chi_mp, T_mp, dt, ED_EGS, ED_corr, ci1, ci2) 
            for dt in dt_list
            ]
        dtE_data, dtC_data = repack(dt_data)
        chi_data = [
            make_plot_dat(alf, chi, T_mp, dt_mp, ED_EGS, ED_corr, ci1, ci2) 
            for chi in chi_max_list
            ]
        chiE_data, chiC_data = repack(chi_data)
        T_data = [
            make_plot_dat(alf, chi_mp, T, dt_mp, ED_EGS, ED_corr, ci1, ci2) 
            for T in T_list
            ]
        TE_data, TC_data = repack(T_data)
        
        dt_alf_E += [dtE_data]
        dt_alf_C += [dtC_data]
        T_alf_E += [TE_data]
        T_alf_C += [TC_data]
        chi_alf_E += [chiE_data]
        chi_alf_C += [chiC_data]
    
    make_subplot(dt_alf_E, dt_list, r"Trotter error, dt", alp_list, 
                 dt_name+"_E")
    make_subplot(dt_alf_C, dt_list, r"Trotter error, dt", alp_list, 
                 dt_name+"_corr")
    make_subplot(T_alf_E, T_list, r"Total time evolved, $\tau$", alp_list, 
                 T_name+"_E")
    make_subplot(T_alf_C, T_list, r"Total time evolved, $\tau$", alp_list, 
                 T_name+"_corr")
    make_subplot(chi_alf_E, chi_max_list, r"Bond dimension, $\chi$", 
                 alp_list, chi_name+"_E")
    make_subplot(chi_alf_C, chi_max_list, r"Bond dimension, $\chi$", 
                 alp_list, chi_name+"_corr")
benchmark_dataread()