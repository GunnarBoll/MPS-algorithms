# -*- coding: utf-8 -*-
"""
Created on Fri May  3 15:29:17 2019

@author: Gunnar
"""
import sys
import os
from numpy import inf

def concatenate_string(stringlist):
    fullword = ''
    for word in stringlist:
        fullword += word
    return fullword

def get_plot_data(*args, **kwargs):
    meas_folder = "/proj/snic2019-8-26/measurements/"
    folder_order = ["tperp", "N", "n", "U", "chi"]
    
    if args == ():
        args = sys.argv[1:]
    
    obser = str(args[0])
    var_param = str(args[1])
    fix_params = [eval(par) if type(par) == str else par for par in args[2:6]]
    
    # If N is the varied parameter the number of bosons usually should be
    # varied with it and is thus treated specially
    if var_param == "N":
        n_ind = folder_order.find('n')
        folder_order.remove('n')
        n = fix_params.pop(n_ind-1)
    
    # Create left/right folder partition list
    var_param_index = folder_order.index(var_param)
    left_folders = folder_order[:var_param_index]
    right_folders = folder_order[var_param_index+1:]
    
    left_params = fix_params[:var_param_index]
    right_params = fix_params[var_param_index:]
    
    # Create the actual left/partition of the folder name as a list e.g.
    # Lp = ['C:/', 'Users/', 'Gunnar/']
    Lp = ([fixpar + "=" + str(left_params[i]) + "/" for i, fixpar in
          enumerate(left_folders)])
    Rp = ['/'] + [fixpar + "=" + str(right_params[j]) + "/"
          for j, fixpar in enumerate(right_folders)]
    
    # Concatenate lists of strings such that the form becomes 
    # Lp = 'C:/Users/Gunnar/'
    
    Lp = concatenate_string(Lp)
    Rp = concatenate_string(Rp)
    print(Lp, Rp)
    
    # os.walk is an iterator which contains a 3-tuple 
    # (top directory, directories, files) at each iteration. It
    # starts in cwd and works its way downward. If a folder fails to contain
    # the relevant observable it is removed from the considered varied
    # parameters.
    for ind, contents in enumerate(os.walk(meas_folder + Lp)):
        if ind == 0:
            var_par_folder_contents = contents[1]
            num_dirs = len(var_par_folder_contents)
        elif ind > num_dirs:
            break
        else:
            if obser + '.dat' not in contents[2]:
                var_par_folder_contents.pop(ind-1)
    
    var_param_vals = [eval(var_par_str.replace(var_param+"=",''))
                      for var_par_str in var_par_folder_contents]
    var_param_vals.sort()
    
    obser_vals = []
    if var_param == 'N':
        file_loc = (lambda p: meas_folder + Lp + 'N=' + str(p) + 'n='
                    + str(p/2 + n) + Rp)
    else:
        file_loc = (lambda p: meas_folder + Lp + folder_order[var_param_index] 
                    + '=' + str(p) + Rp)
    
    # Collect the data specified in variable "obser"
    for val in var_param_vals:
        file = file_loc(val) + obser + ".dat"
        with open(file, 'r') as fr:
            obser_vals.append(eval(fr.readline()))
    
    return var_param_vals, obser_vals
