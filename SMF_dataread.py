"""
Reads files produced with stat_MF.py and plots the order parameter versus the
interaction strength, U. The program reads files with combinations of N and U.
"""
from matplotlib import pyplot as plt
import importlib as imp

import storage as st

imp.reload(st)

# Function fetches the last line of a file
def tail(filename):
    with open(filename, 'r') as fr:
        lines = fr.read().splitlines()
        last_line = lines[-1]
    return float(last_line)

# Main program,for N in N_list and U in U_list calls files and retrieves the
# order parameter for each U. Plots the order parameter versus U.
def main():
    data_direc = "C:/Users/Gunnar/Documents/Ph.D/Data/Static_MF/"
    N_list = [20, 30, 40, 50, 60]
    for N in N_list:
        name = "SMF_"+"N="+str(N)+"_1/"
        direc = data_direc + name
        U_list = [0., 0.5, 1., 1.5, 2., 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]

        
        files = [direc+"N="+str(N)+",U="+str(U)+".txt" for U in U_list]
        ord_par = [tail(file) for file in files]
        
        plt.plot(U_list, ord_par)
    plt.show()
    
    return

main()