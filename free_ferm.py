
import numpy as np


N = 30
mu = 0
J = -1.0/2


H = np.zeros([N,N])

for i in range(N-1):
    H[i,i] = mu
    H[i, i + 1] = J
    H[i + 1, i] = J
H[N-1,N-1] = mu

H += mu*np.eye(N)/2

e, v = np.linalg.eigh(H)

E_GS = 0
for ener in e:
    if not ener < 0:
        break
    E_GS += ener


print(E_GS)