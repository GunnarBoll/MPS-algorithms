
import numpy as np
import matplotlib.pyplot as plt

direc = "C:/Users/Gunnar/Documents/Ph.D/Learning/DMRG/Tryout code/output/"
f = open(direc + "t_perp=0.0/2018-12-19_1.txt", "r")
E_GS = []
corr = []

for line in f:
    if line.startswith("#"):
        com = line.strip("# ")
        com = com.strip("\n")
        continue

    if com == "GS Energy":
        E_GS.append(float(line))
    elif com == "Correlators":
        corr.append(float(line))

f.close()

N = int(np.sqrt(len(corr) / 2))
corr_algo = np.reshape(np.array(corr[: N**2]), (N, N))
corr_ED = np.reshape(np.array(corr[N**2 :]), (N, N))
E_algo = E_GS[0]
E_ED = E_GS[1]
E_err = np.abs(E_algo - E_ED)

corr_err = abs(corr_algo - corr_ED)
print(corr_err)
print(E_err)