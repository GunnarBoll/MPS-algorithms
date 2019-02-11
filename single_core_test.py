
import multiprocessing as mp
import time as t
import threading as th
import itertools
import numpy as np
import importlib as imp
import os

imp.reload(np)

thr_num = '1'
os.environ['MKL_NUM_THREADS'] = thr_num
os.environ['OPENBLAS_NUM_THREADS'] = thr_num

def f(x):
    return x*x

def thing(L):
    res = 0
    for i in range(len(L)):
        for j in range(len(L)):
            res += L[i] * L[j]
    return res
    
def main():

    N = 100
    A = np.random.rand(N,N)
        
    if __name__ == '__main__':
        begin = t.time()
        pool = mp.Pool(2)
        result = pool.map_async(thing, list(itertools.repeat(A, 20)))
        pool.close()
        pool.join()
        
        dat = [stuff for stuff in result.get()]
        print(dat)
        finish = t.time()
        print("\nCore test total time:", finish-begin)
    
__spec__ = None
main()

# threads = []
# th_num = 2
# begin = t.time()
# for i in range(th_num):
#     threads.append(th.Thread(target=thing, args=([listy])))
# for i in range(th_num):
#     threads[i].start()
# 
# for i in range(th_num):
#     threads[i].join()
# 
# finish = t.time()
# 
# print("Threading time:", finish-begin)