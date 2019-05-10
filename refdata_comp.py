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
    mu_l = [-2.0, -1.5, -1., -0.5, 0., 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0,
            4.5, 5.0, 5.5, 6.0, 6.5, 7.0]
    
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
    U_list = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75,
              3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 3.02, 3.04,
              3.06, 3.08, 3.1, 3.12, 3.14, 3.16, 3.18, 3.20, 3.22, 3.24]
    # U_list = [3.01 + i/50 for i in range(12)]
    
    for U in U_list:
        for orp in orp_list:
            cpden_compile(orp/10, U, tperp, run_nr)



compile_all(0.05, 1)
