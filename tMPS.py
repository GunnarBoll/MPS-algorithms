"""tMPS algorithm solving an antiferromagnetic transverse Ising chain
Written by: Gunnar Bollmark
"""
import numpy as np
import storage
import imp
import ED_AFTIC as ed
imp.reload(storage)
imp.reload(ed)

### MAIN PROGRAM #######################################
def main():
    #Model parameters Typical input: 1.0 1.0 0.01 2 8 4 second tDMRG
    
    f = open("C:/Users/Gunnar/Documents/Ph.D/Learning/DMRG/Tryout code/TEBD/params.txt","r")
    inputlist = []
    for line in f:
        inputlist.append(line.rstrip('\n'))
    f.close()
    params = []
    for i in range(len(inputlist)):
        if i<3:
            params.append(float(inputlist[i]))
        elif i>2 and i < 6:
            params.append(int(inputlist[i]))
        else:
            params.append(inputlist[i])
    J = params[0]        #Spin-spin coupling
    h = params[1]        #External field coupling
    dt = params[2]       #Time step
    d = params[3]          #One-particle Hilbert space dim
    chi = params[4]        #Maximum MPS dim      
    N = params[5]          #Site number
    step_number = int(30/dt) #Number of time steps
    order = params[6]    #Trotter order
    algo = params[7]     #Which algorithm to use
    
    #H belongs to class Hamiltonian which features time evolution and model construction
    H = storage.Hamiltonian(J,h,N,dt,d,chi,which= "AFTIC",TO = order)
    #Energy operators on each bond
    Hchain = H.Hchain
    
    # Build the initial state #
    Psi = storage.StateChain(N, chi, d, algo)
        
    #Time evolve the state
    Psi = H.time_evolve(Psi, step_number, algo)

    
    #Energy calculation
    #Using a two-site Hamiltonian at each bond and find total energy per site
    E = Psi.get_ener(Hchain)
    print("GS energies")
    print("GS energy from algorithm",algo,": ",sum(E))
    
    #Exact diagonalization trotter
    ED = ed.ExactD(J,h,N,dt,d, TO = order)
    ED.trotter_time(step_number)
    E_ed = ED.get_ener(ED.trot_state)
    #Exact GS energy
    ED.exact_GS()
    print("GS energy from ED Trotter: ", sum(E_ed))
    print("GS energy: ", ED.E_GS)
    
    print("\nBond energies")
    print("Algorithm bond energies: ",E)
    print("\nED trotter bond energies: ", E_ed)
    
    #Build product state from MPS
    B = Psi.B
    psi = B[0]
    for i in range(len(B)-1):
        psi = np.tensordot(psi,B[i+1], (i+2,1))
    
    psi = np.reshape(psi, (d**N))
    print("Algorithm energy from product state: ",sum(ED.get_ener(psi)))
    
    #Comparison of MPS and ED trotter GS
    print("\nNon-zero components of product states")
    for i in range(len(psi)):
        if abs(psi[i])>10**-8 or abs(ED.trot_state[i]) > 10**-8:
            
            print("Algorithm: ", psi[i]," ED trotter: ",ED.trot_state[i])
            
    print("\nOther output")
    print("Lowest energies ",ED.elist)
    ED.exact_time(step_number)
    print("Exact time evolution energy ",ED.E_evol)
    print("Truncation error: ",Psi.err)
    print("Final maximum bond dimension: ",H.chi)
    
    
    return 0
    
main()
