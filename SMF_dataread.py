"""
Reads files produced with stat_MF.py and plots the order parameter versus the
interaction strength, U. The program reads files with combinations of N and U.
"""
from matplotlib import pyplot as plt
import importlib as imp
import scipy as sci
import numpy as np

import storage as st
from GS_get import get_GS

imp.reload(st)

# Function fetches the last line of a file
def tail(filename):
    with open(filename, 'r') as fr:
        lines = fr.read().splitlines()
        last_line = lines[-1]
    return float(last_line)

def get_err(filename):
    with open(filename, 'r') as fr:
        err = [line for i, line in enumerate(fr) if i == 3]
    return float(err[0])

def truncerr_extrap(orp, gs_file):
    Psi = get_GS(gs_file)
    trunc_errs = [Psi.err]
    inc = int((Psi.chi-20) / 4)
    chi_list = range(Psi.chi, 20, -inc)
    trunc_errs += [err_gen(chi, Psi) for chi in chi_list]
    
    print(trunc_errs)
    return trunc_errs

def err_gen(chi, Psi):
    Psi.chi = chi
    step_num = 100
    H = st.Hamiltonian(Psi.g1, Psi.g2, Psi.N, 0.1, Psi.d, Psi.chi, 
                       model="HCboson", TO="fourth")
    Psi = H.time_evolve(Psi, step_num, "tDMRG", fast_run=True)
    
    return Psi.err

def isize_orp(N_list, orp_list, i, U):
    orp_list.reverse()
    N_list.reverse()
    
    inv_N_list = [1/N for N in N_list]
    p = sci.polyfit(inv_N_list, orp_list, 2)
    
    fit = lambda x: p[0]*x**2 + p[1]*x +p[2]
    x_list = np.linspace(0, 0.05, 20)
    
    fig = plt.figure(i)
    dat_plot, = plt.plot(inv_N_list, orp_list, "ro", label="Data points")
    quad_plot, = plt.plot(x_list, fit(x_list), label="Quadratic fit")
    plt.ylabel("<a>")
    plt.xlabel("1/L")
    plt.legend([dat_plot, quad_plot])
    # fig.savefig("U="+str(U)+".png")
    
    N_list.reverse()
    orp_list.reverse()
    return p[-1]

def finsiz_extrap():
    pass

def namer(N, U):
    name = "N=" + str(N) + ",U=" + str(U) + ".txt"
    return name
# Main program,for N in N_list and U in U_list calls files and retrieves the
# order parameter for each U. Plots the order parameter versus U.
def main():
    data_direc = "C:/Users/Gunnar/Documents/Ph.D/Data/Static_MF/"
    N_list = [30, 40, 50, 60]
    U_list = [0., 0.5, 1.0, 1.5, 2., 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
    orps_vs_U = []
    for N in N_list:
        name = "SMF_"+"N="+str(N)+"_1/"
        direc = data_direc + name
        
        files = [direc+namer(N,U) for U in U_list]
        GS_files = [direc+"GS_"+namer(N,U) for U in U_list]
        
        # trunc_errs = [get_err(file) for file in GS_files]
        
        ord_par = [tail(file) for file in files]
        
#        new_orp = [truncerr_extrap(ord_par[i], GS_files[i]) for i in
#                   range(len(ord_par))]
        
        new_orp = ord_par
        
        orps_vs_U.append(new_orp)
        plt.figure(19)
        plt.plot(U_list, ord_par)
#        plt.figure(2)
#        plt.plot(U_list, new_orp)
    orps_vs_N = [[orps_vs_U[j][i] for j in range(len(N_list))] for i in 
                  range(len(U_list))]
    isize_orp_vs_U = [isize_orp(N_list, orps_vs_N[i], i, U_list[i]) for i in 
                      range(len(U_list))]
    print(isize_orp_vs_U)
    plt.figure(20)
    plt.plot(U_list, isize_orp_vs_U)
    plt.ylabel("<a>")
    plt.xlabel("U")
    
    
    plt.show()
    
    return

main()