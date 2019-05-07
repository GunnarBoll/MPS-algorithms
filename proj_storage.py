# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 10:42:33 2019

@author: Gunnar
"""

import pathlib

def sum_string(strl):
    summed = ''
    for item in strl:
        summed += item + '/'
    return summed

def proj_store(dname, fname, data, replace=False):
    if type(data) != list:
        data = [data]
    
    folders = dname.split("/")
    root = folders[0]
    rest_name = sum_string(folders[1:])
    
    proj = "/proj/snic2019-8-26/"
    direc = proj + dname
    run_nr = 1
    while True:
        if run_nr > 20:
            break
        try:        
            pathlib.Path(direc).mkdir(parents=True, exist_ok=True)
            with open(direc+fname, "x") as fw:
                for mat in data:
                    fw.write(str(mat) + "\n")
            break
        except FileExistsError:
            if replace:
                with open(direc+fname, 'w') as fw:
                    for mat in data:
                        fw.write(str(mat) + "\n")
                break
            else:
                direc = proj + root + "_" + str(run_nr) + "/" + rest_name
    
    return