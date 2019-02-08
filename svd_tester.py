import numpy as np
import time
import importlib as imp

MKL_NUM_THREADS = '1'
OMP_NUM_THREADS = '1'
imp.reload(np)

def main():
    
    N = 100
    A = np.ones([N, N])
    start = time.time()
    U, S, V = np.linalg.svd(A)
    end = time.time()
    print("Time cost for SVD:", end-start)

main()