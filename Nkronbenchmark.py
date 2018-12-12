import numpy as np
import time as t

N = 20000
a = np.eye(2)

start = t.process_time()
c = np.zeros([2*N,2*N])
for i in range(0,2*N,2):
    c[i:i+2,i:i+2] = a
end = t.process_time()
print(end-start)