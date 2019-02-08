import numpy as np
import math

A = np.arange(12).reshape(3,4)

ei1, evec1 = np.linalg.eigh(np.matmul(A.T.conj(), A))
ei2, evec2 = np.linalg.eigh(np.matmul(A, A.T.conj()))
u, s, v = np.linalg.svd(A)

ei1 = np.flip(ei1, 0)
ei2 = np.flip(ei2, 0)

val = [ei1[i] for i in range(len(ei2)) if math.isclose(ei2[i], ei1[i])]

s = s[:(np.sum(s > 10**-10))]

print(ei1, ei2, s, np.sqrt(val))

print(u, np.flip(evec2, 1))
print(v.T, np.flip(evec1, 1))
