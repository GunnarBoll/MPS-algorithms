"""
Takes reference data produced bt get_chemdens.py and restructures into files
with reference to order parameter sizes. Needs to be run once to obtain usable
reference data in stat_MF.py.
"""
import os

# Recompiles data in a directory given by U, tperp and run_nr, recompiles it
# into a single file and labels it with ord_par
def cpden_compile(ord_par, U, tperp, run_nr):
    dname = "C:/Users/Gunnar/Documents/Ph.D/Data/Density_ref/"
    dname += ("rho_of_mu_tperp=" + str(tperp) + ",U=" + str(U) + "_"
              + str(run_nr) + "/")
    mu_l = [-1., -0.5, 0., 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
    
    file_names = ["ordpar=" + str(ord_par) + ",mu=" + str(mu) for mu in mu_l]
    dens_list = []
    
    for filename in file_names:
        with open(dname+filename + ".txt", "r") as fr:
            dens_list += [float( fr.readline().strip("\n") )]
        
    
    new_file = "ord_par=" + str(ord_par)
    
    data = dens_list + mu_l
    
    with open(dname+new_file + ".txt", "x") as fw:
        for dps in data:
            fw.write(str(dps) + "\n")
    
    for filename in file_names:
        os.remove(dname+filename+".txt")
    
# Main program, for a certain tperp and run_nr recompiles data for a set of U
def compile_all(tperp, run_nr):
    orp_list = range(11)
    U_list = [3.0, 3.5, 4.0, 4.5, 5.0]
    
    for U in U_list:
        for orp in orp_list:
            cpden_compile(orp/10, U, tperp, run_nr)



compile_all(0.05, 1)
