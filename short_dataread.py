# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 18:23:13 2019

@author: Gunnar
"""
import matplotlib.pyplot as plt
import scipy as sci
import scipy.optimize as sciop

def exp_fitfunc(x, a, b, c):
    return a*sci.exp(-b*x) + c

def main():
    file ="C:/Users/Gunnar/Documents/SMF_N=80/N=80,U=3.22.txt"
    with open(file, 'r') as fr:
        lines = [eval(line.strip('\n')) for line in fr]
    
    N = int(len(lines))
    plt.close("all")
    old_lines = lines[int(N/2):]
    old_x = [i for i in range(25,276)]
    lines = lines[250:]
    xdata = [ i for i in range(len(lines))]
    p, cov = sciop.curve_fit(exp_fitfunc, xdata, lines)
    
    xx = sci.linspace(25, 275, 1000)
    plt.plot(xx, exp_fitfunc(xx, p[0], p[1], p[2]))
    plt.plot(old_x, old_lines)
    
    print(p[2])
    
    return
main()