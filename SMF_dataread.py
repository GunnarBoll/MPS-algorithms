from matplotlib import pyplot as plt
import importlib as imp
import numpy as np

import storage as st
from GS_get import get_GS

imp.reload(st)

def tail(filename):
    
    with open(filename, 'r') as fr:
        lines = fr.read().splitlines()
        last_line = lines[-1]
    
    return float(last_line)

def main():
    data_direc = "C:/Users/Gunnar/Documents/Ph.D/Data/Static_MF/"
    N_list = [20]
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