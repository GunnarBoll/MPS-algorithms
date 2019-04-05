# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 08:59:24 2019

@author: Gunnar
"""

import subprocess as subp
import os
import pathlib

class MPTKState:
    def __init__(self, dname, run_script=None, params=None, model=None,
                 cluster=False):
        self.loc = dname + "/"
        self.home = os.environ.get("HOME") + "/"
        self.proj = "/proj/snic2019-8-26/"
        self.tmp_dir = os.environ.get("SNIC_TMP") + "/"
        self.cluster = cluster
        
        if cluster:
            self.write_direc = self.tmp_dir
        else:
            self.write_direc = self.home
        
        if os.path.exists(self.proj + self.loc):
            self.finish = True
            self.direc = self.proj + self.loc
        else:
            self.finish = False
            self.direc = self.write_direc + self.loc
            self.call_string = [self.direc]
        
        self.params = params
        self.run_script =  run_script
        self.model = model
        
        if model == "SMF":
            self.N = params[0]
            self.U = params[1]
            self.mu = params[2]
            self.alp = params[3]
            self.chi = params[4]
            self.call_string += [self.N, self.U, self.mu, self.alp, self.chi]
                
        return
    
    def mptk_run(self):
        if not self.finish:
            scr_name = self.run_script
            stdo = self.bash_call(scr_name, self.call_string)
            self.finish = True
        else:
            print("There already exists an MPTK folder")
    
    def expec(self, oper, loc):
        state = self.direc + "GS_file.psi." + "1"
        oper_call = self.direc + "lattice:" + oper + "(" + str(loc) + ")"
        expec_scr = "bin/mp-expectation"        
        exval = eval(self.bash_call(expec_scr, [state, oper_call]))[0]
        
        return exval
    
    def correl(oper1, oper2, loc1, loc2):
        pass
    
    def copy_solution(self, new_name):
        cmd = "cp"
        loc = self.proj + new_name
        
        if not os.path.exists(loc):
            pathlib.Path(loc).mkdir(parents=True, exist_ok=False)
        
        self.bash_call(cmd, ["-r", self.direc, loc], with_home=False)
    
    def delete_solution(self):
        self.bash_call("rm", ["-r", self.write_direc+self.loc])
    
    def get_trunc_err(self):
        file = self.direc + "GS_file.sweep"
        with open(file, "r") as fr:
            lines = fr.read().splitlines()
            last = lines[-1]
            trunc_err = last.split(" ")[6]
        return trunc_err
        
    def bash_call(self, scr_name, argv, with_home=True):
        if with_home:
            call_name = self.home + scr_name
        else:
            call_name = scr_name
        
        arg_list = [call_name] + [str(arg) for arg in argv]
        
        res = subp.check_output(arg_list)
        
        return res