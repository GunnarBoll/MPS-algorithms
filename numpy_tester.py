import numpy as np
import time as t
import imp

imp.reload(np)


N = 5000
A = np.random.rand(N,N)
for i in range(N):
    for j in range(N):
        A[i, j] = A[j, i]

start = t.time()
res = np.linalg.eigh(A)
end= t.time()

print("Time:", end-start)