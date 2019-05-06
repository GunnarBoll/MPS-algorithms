# -*- coding: utf-8 -*-
"""
Created on Fri May  3 15:29:17 2019

@author: Gunnar
"""
import sys
import os

from cwd_storage import cwd_store

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
    fix_params = [eval(par) for par in args[2:6]]
    
    # Create left/right folder partition list
    var_param_index = folder_order.index(var_param)
    left_folders = folder_order[:var_param_index]
    right_folders = folder_order[var_param_index+1:]
    
    left_params = fix_params[:var_param_index]
    right_params = fix_params[var_param_index:]
    
    # Create the actual left/partition of the folder name
    Lp = ([fixpar + "=" + str(left_params[i]) + "/" for i, fixpar in
          enumerate(left_folders)])
    Rp = ['/'] + [fixpar + "=" + str(right_params[j]) + "/"
          for j, fixpar in enumerate(right_folders)]
    
    Lp = concatenate_string(Lp)
    Rp = concatenate_string(Rp)
    
    var_par_folder_contents = os.listdir(meas_folder + Lp)
    var_param_vals = [eval(var_par_str.replace(var_param+"=",''))
                      for var_par_str in var_par_folder_contents]
    var_param_vals.sort()
    
    obser_vals = []
    file_loc = (lambda p: meas_folder + Lp + folder_order[var_param_index] 
                + '=' + str(p) + Rp)
    
    # Collect the data specified in variable "obser"
    for val in var_param_vals:
        file = file_loc(val) + obser + ".dat"
        with open(file, 'r') as fr:
            obser_vals.append(eval(fr.readline()))
    
    cwd_store('plotdat', 'ydat', obser_vals)
    cwd_store('plotdat', 'xdat', var_param_vals)
    
    return var_param_vals, obser_vals

get_plot_data()