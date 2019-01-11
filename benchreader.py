
import numpy as np
import datetime
from matplotlib import pyplot as plt

def figcloser():
    for fignr in range(9, 0, -1):
        plt.close(fignr)

def main():
    date = "2019-01-10"
    run_num = "_run#1/"
    direc = ("D:/Documents/Ph.D/Learning/DMRG/Tryout code/output/"
             + date + run_num)
    
    g1 = [1.0, 1.0]
    g2 = [1.0, 0.01]
    dt_list = [0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1]
    T_list = [20, 24, 28, 30, 34, 40, 65, 100]
    T_list.reverse()
    chi_max_list = [10, 12, 14, 16, 18, 20, 30, 50, 70]
    chi_max_list.reverse()
    bis_err_list = [10**-7, 5*10**-8, 10**-8, 5*10**-9, 10**-9, 5*10**-10,
                    10**-10, 5*10**-11, 10**-11]
    bis_err_list.reverse()
    N = 20
    
    E_lists = [[], [], [], []]
    trunc_err_lists = [[], [], [], []]
    corr_lists = []
    
    with open(direc + "ED.txt", "r") as fr:
        read_mat = []
        for line in fr:
            read_mat.append(float(line.strip("\n")))
        ED_GSE = read_mat[0]
        ED_corr = np.array(read_mat[1 :]).reshape(N,N)
        fr.close()
    
    iter_lists = [dt_list, T_list, chi_max_list, bis_err_list]
    for ind in range(4):
        params = [dt_list[0], T_list[0], chi_max_list[0], bis_err_list[0]]
        param_list = iter_lists[ind]
        E = E_lists[ind]
        trunc_err = trunc_err_lists[ind]
        corr = []
        for params[ind] in param_list:
            dt = params[0]
            T = params[1]
            chi_max = params[2]
            bis_err = params[3]
            name = ("chi=" + str(chi_max) + ",T=" + str(T) + ",dt="
                    + str(dt) + ",BE=" + str(bis_err))
            with open(direc+name+".txt", "r") as fr:
                reading_mat = []
                for line in fr:
                    reading_mat.append(float(line.strip("\n")))
                E.append(reading_mat[0])
                trunc_err.append(reading_mat[1])
                corr.append(np.array(reading_mat[2 :]).reshape(N,N))
                fr.close()
        corr_lists.append(np.array(corr))
    
    figcloser()
    for i in range(4):
        ener_err = abs(np.array(E_lists[i]) - np.array(ED_GSE))
        plt.figure(i+1)
        if i != 3:
            plt.plot(iter_lists[i], ener_err)
        else:
            plt.loglog(iter_lists[i], ener_err)
    
    for j in range(4):
        corr_err = abs(corr_lists[j][:, 10, 15] - ED_corr[10, 15])
        plt.figure(j+5)
        if j != 3:
            plt.plot(iter_lists[j], corr_err)
        else:
            plt.loglog(iter_lists[j], corr_err)
    
    plt.figure(9)
    plt.plot(trunc_err_lists[2][:len(E_lists[2])-2],
             E_lists[2][:len(E_lists[2])-2],
             trunc_err_lists[2][:len(E_lists[2])-2],
             E_lists[2][:len(E_lists[2])-2],
             "ko")
    plt.xlabel(r"$\chi$")
    plt.ylabel(r"$E_{GS}$")
    plt.show()
    
    return 0
    
main()