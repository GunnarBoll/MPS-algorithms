# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 08:59:24 2019

@author: Gunnar
"""

import subprocess as subp
import os

class MPTKState:
    def __init__(self, dname, params, model, cluster):
        self.finish = False
        
        self.home = os.environ.get("HOME") + "/"
        
        if cluster:
            self.write_direc = os.environ.get("SNIC_TMP") + "/"
        else:
            self.write_direc = self.home + "mptk_states/"
        
        self.direc = self.write_direc + dname
        self.model = model
        
        if model == "SMF":
            self.N = params[0]
            self.U = params[1]
            self.mu = params[2]
            self.alp = params[3]
            self.call_string = [self.direc, self.N, self.U, self.mu, self.alp]
        else:
            raise ValueError("Not a valid model")

    def mptk_run(self):
        if not self.finish:
            scr_name = "bin/mptk_script.sh"
            stdo = self.bash_call(scr_name, self.call_string)
            self.finish = True
        else:
            print("There already exists an MPTK folder")
    
    def expec(self, oper, loc):
        state = self.direc + "/GS_file.psi." + "12"
        oper_call = self.direc + "/lattice:" + oper + "(" + str(loc) + ")"
        expec_scr = "bin/mp-expectation"        
        exval = eval(self.bash_call(expec_scr, [state, oper_call]))[0]
        
        return exval
    
    def correl(oper1, oper2, loc1, loc2):
        pass
    
    def copy_solution(self, new_name):
        cmd = "cp"
        loc = self.home + new_name
        
        self.bash_call(cmd, ["-r", self.direc, loc], with_home=False)
    
    def delete_solution():
        pass

    def bash_call(self, scr_name, argv, with_home=True):
        if with_home:
            call_name= self.home + scr_name
        else:
            call_name = scr_name
        
        arg_list = [call_name] + [str(arg) for arg in argv]
        
        res = subp.check_output(arg_list)
        
        return res