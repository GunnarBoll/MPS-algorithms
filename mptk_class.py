# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 08:59:24 2019

@author: Gunnar
"""

import subprocess as subp
import os

class MPTKState:
    def __init__(self, dname, params, model):
        self.finish = False
        self.direc = "/home/gunbo249/mptk_states/" + dname
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
            stdo = bash_call(scr_name, self.call_string)
            self.finish = True
        else:
            print("There already exists an MPTK folder")
    
    def expec(self, oper, loc):
        oper_call = self.direc + "/lattice:" + oper + "(" + str(loc) + ")"
        expec_scr = "bin/mp-expectation"        
        exval = eval(bash_call(expec_scr, [oper_call])[0])
        
        return exval
    
    def correl(oper1, oper2, loc1, loc2):
        pass
    
    def delete_solution():
        pass

def bash_call(scr_name, argv):
    arg_list = ["/home/gunbo249/" + scr_name] + [str(arg) for arg in argv]
    
    res = subp.check_output(arg_list)
    
    return res