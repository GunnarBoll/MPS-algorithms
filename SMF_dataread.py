"""
Reads files produced with stat_MF.py and plots the order parameter versus the
interaction strength, U. The program reads files with combinations of N and U.
"""
from matplotlib import pyplot as plt
import importlib as imp

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
    chi_list = range(20, Psi.chi, 10)
    trunc_errs += [err_gen(chi, Psi) for chi in chi_list]
    
    pass

def err_gen():
    pass

def finsiz_extrap():
    pass

# Main program,for N in N_list and U in U_list calls files and retrieves the
# order parameter for each U. Plots the order parameter versus U.
def main():
    data_direc = "C:/Users/Gunnar/Documents/Ph.D/Data/Static_MF/"
    N_list = [20]
    for N in N_list:
        name = "SMF_"+"N="+str(N)+"_1/"
        direc = data_direc + name
        U_list = [0., 0.5, 1., 1.5, 2., 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]

        
        files = [direc+"N="+str(N)+",U="+str(U)+".txt" for U in U_list]
        GS_files = ["GS_"+name for name in files]
        
        # trunc_errs = [get_err(file) for file in GS_files]
        
        ord_par = [tail(file) for file in files]
        
        new_orp = [truncerr_extrap(ord_par[i], GS_files[i]) for i in
                   range(len(ord_par))]
        plt.figure(1)
        plt.plot(U_list, ord_par)
        plt.figure(2)
        plt.plot(U_list, new_orp)
    plt.show()
    
    return

main()