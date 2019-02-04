
import multiprocessing as mp
import time as t
import threading as th
import itertools

def f(x):
    return x*x

def thing(L):
    res = 0
    for i in range(len(L)):
        for j in range(len(L)):
            res += L[i] * L[j]
    return res

listy = list(itertools.repeat(1, 1000))
if __name__ == '__main__':
    __spec__ = None
    begin = t.time()
    
    pool = mp.Pool(2)
    result = pool.map_async(thing, [listy, listy])
    pool.close()
    pool.join()
    
    finish = t.time()
    
    print("\nCore test total time:", finish-begin)

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