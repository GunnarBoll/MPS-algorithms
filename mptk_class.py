# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 08:59:24 2019

@author: Gunnar
"""

import subprocess as subp
import os
import pathlib
import re

# Class that features functions for ease of use with the MPToolkit
class MPTKState:
    def __init__(self, dname, run_script=None, params=None, model=None,
                 cluster=False):
        self.loc = dname
        self.home = os.environ.get("HOME") + "/"
#        self.proj = "/proj/snic2019-8-26/"
        self.proj = "/proj/uppstore2019070/"
        self.tmp_dir = os.environ.get("SNIC_TMP") + "/"
        self.cluster = cluster
        
        # Operators (to be made more modal)
        self.orp = "B"
        self.num = "N"
        
        if cluster:
            self.write_direc = self.tmp_dir
        else:
            self.write_direc = self.home
        
        # Check if the GS already exists
        if os.path.exists(self.proj + self.loc):
            self.finish = True
            self.direc = self.proj + self.loc + "/"
            
        else:
            self.finish = False
            self.direc = self.write_direc + self.loc + "/"
            
        self.call_string = [self.direc]
        self.params = params
        self.run_script =  run_script
        self.model = model
        
        # Either retreive GS setup for creation
        if self.finish:
            self.get_hamil()
            self.get_params()
        else:
            if model == "SMF":
                self.N = params[0]
                self.U = params[1]
                self.mu = params[2]
                self.alp = params[3]
                self.chi = params[4]
                self.ener_op = ("H+" + str(self.U) + "*H_V-" + str(self.mu)
                                + "*H_mu+" + str(self.alp) + "*H_SMF")
                self.call_string += [self.N, self.U, self.mu, self.alp,
                                     self.chi]
                
        return
    
    # Retrieves the model parameters
    def get_params(self):
        split_direc = re.split("/|=", self.loc)
        
        N_ind = split_direc.index("N") + 1
        self.N = eval(split_direc[N_ind])
        
        chi_ind = split_direc.index("chi") + 1
        self.chi = eval(split_direc[chi_ind])
        
        
    
    # Retrieves the exact Hamiltonian operator (string) of the ground state
    def get_hamil(self):
        hamils = ["H_V", "H_mu", "H_SMF"]
        params = []
        with open(self.direc+"GS_file.params", 'r') as fr:
            call = fr.read()
            ind = call.find(":H") + 2
            call = call[ind:]
            for ham in hamils:
                ind = call.find(ham)
                params.append(eval(call[:ind-1]))
                call = call[ind+len(ham)+1 :]
            
            self.U = params[0]
            self.mu = params[1]
            self.alp = params[2]
        
        self.ener_op = ("H+" + str(self.U) + "*H_V-" + str(self.mu)
                        + "*H_mu+" + str(self.alp) + "*H_SMF")
        
    
    # Retrieves the state energy with the string "ener_op"
    def get_ener(self):
        E = self.expec(self.ener_op)
        return E
    
    # Retrieves the order parameter by averaging over the expectation value
    # of the annihilation operator over the central sites
    def get_orp(self):
        N=self.N
        start = int(N/2)-int(N/4) + 1
        end = int(N/2) + int(N/4) + 1 if N%2==0 else int(N/2) + int(N/4) + 2
        
        orp = 0
        av_num_sites = 0
        op_func = lambda ind: self.orp + "(" + str(ind) + ")"
            
        for orpind in range(start, end):
            av_num_sites += 1
            orp += abs(self.expec(op_func(orpind)))
        
        orp /= av_num_sites
        return orp
    
    def get_dens(self):
        oper = lambda ind: self.num + "(" + str(ind) + ")"
        dens = 0
        for densind in range(1, self.N+1):
            dens += self.expec(oper(densind))
        dens /= self.N
        return dens

    # Runs the algorithm with given parameters. Uses a list of arguments
    # required in the bash function call
    def mptk_run(self):
        if not self.finish:
            scr_name = self.run_script
            stdo = self.bash_call(scr_name, self.call_string)
            self.finish = True
            return False
        else:
            print("There already exists an MPTK folder")
            return True
    
    # Returns the expectation value of an operator (string).
    def expec(self, oper):
        state = self.direc + "GS_file.psi." + "1" #1 is the name of the state
        oper_call = self.direc + "lattice:" + oper
        expec_scr = "bin/mp-expectation"        
        exval = eval(self.bash_call(expec_scr, [state, oper_call]))[0]
        
        return exval
    
    # Would return correlators
    def correl(oper1, oper2, loc1, loc2):
        pass
    
    # Copies the solution from the working directory into a directory of given
    # name.
    def copy_solution(self, new_name):
        cmd = "cp"
        loc = self.proj + new_name
        
        if not os.path.exists(loc):
            pathlib.Path(loc).mkdir(parents=True, exist_ok=False)
        
        self.bash_call(cmd, ["-r", self.direc, loc], with_home=False)
    
    # Deletes the solution in the local directory
    def delete_solution(self):
        self.bash_call("rm", ["-r", self.direc], with_home=False)
    
    # Retrieves truncation error from the MPToolkit files
    def get_trunc_err(self):
        file = self.direc + "GS_file.sweep"
        with open(file, "r") as fr:
            lines = fr.read().splitlines()
            last = lines[-1]
            trunc_err = last.split(" ")[6]
        
        trunc_err.replace('\n', '')
        trunc_err = eval(trunc_err)
        return trunc_err
    
    # Makes a call to bash terminal with function in scr_name and options in
    # argv.
    def bash_call(self, scr_name, argv, with_home=True):
        if with_home:
            call_name = self.home + scr_name
        else:
            call_name = scr_name
        
        arg_list = [call_name] + [str(arg) for arg in argv]
        
        res = subp.check_output(arg_list)
        
        return res