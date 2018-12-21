
import numpy as np
import datetime
import matplotlib.pyplot as plt

def main():
    date = "2018-12-20"
    direc = ("C:/Users/Gunnar/Documents/Ph.D/Learning/DMRG/Tryout code/output/"
             + date + "_run#1/")
    
    g1 = [1.0, 0.0]
    g2 = [0.0, 1.0]
    dt_list = [0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2]
    T_list = [20, 30, 40, 50, 60, 70, 80, 90, 100]
    T_list.reverse()
    chi_max_list = [20, 25, 30, 35, 40, 45, 50]
    chi_max_list.reverse()
    bis_err_list = [10**-6, 10**-7, 10**-8, 10**-9, 10**-10]
    bis_err_list.reverse()
    N = 4
    
    E_lists = [[], [], [], []]
    EDE_lists = [[], [], [], []]
    corr_lists = [[], [], [], []]
    EDcorr_lists = [[], [], [], []]
    
    iter_lists = [dt_list, T_list, chi_max_list, bis_err_list]
    for ind in range(4):
        params = [dt_list[0], T_list[0], chi_max_list[0], bis_err_list[0]]
        param_list = iter_lists[ind]
        E = E_lists[ind]
        EDE = EDE_lists[ind]
        corr = corr_lists[ind]
        EDcorr = EDcorr_lists[ind]
        for params[ind] in param_list:
            dt = params[0]
            T = params[1]
            chi_max = params[2]
            bis_err = params[3]
            name = ("chi=" + str(chi_max) + ",T=" + str(T) + ",dt="
                    + str(dt) + ",BE=" + str(bis_err))
            with open(direc+name+".txt", 'r') as fr:
                reading_mat = []
                for line in fr:
                    reading_mat.append(float(line.strip('\n')))
                E.append(reading_mat[0])
                EDE.append(reading_mat[1])
                corr = reading_mat[2 : (N**2)+2]
                EDcorr = reading_mat[(N**2)+2 :]
                fr.close()
        corr = np.array(corr).reshape(N,N)
        EDcorr = np.array(EDcorr).reshape(N,N)
    
    plt.close
    for i in range(4):
        plotstuff = np.array(E_lists[i]) - np.array(EDE_lists[i])
        plt.figure(i+1)
        plt.plot(iter_lists[i], plotstuff)
        plt.show()
    return 0
    
main()