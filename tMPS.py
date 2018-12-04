"""tMPS algorithm solving an antiferromagnetic transverse Ising chain
Written by: Gunnar Bollmark
"""
import numpy as np
import storage as st
import imp
import ExactDiag as ed
imp.reload(st)
imp.reload(ed)

### MAIN PROGRAM #######################################
def main():
    #Model parameters Typical input: 1.0 1.0 0.01 2 8 4 second tDMRG
    
    f = open("C:/Users/Gunnar/Documents/Ph.D/Learning/DMRG/Tryout code/TEBD/params.txt","r")
    inputlist = []
    keylist = {}
    for line in f:
        if line.startswith('#'):
            continue
        inputlist.append(line.rstrip('\n').split(" = "))
        x = inputlist[-1][0]
        y = inputlist[-1][1]
        try:
            d = int(y)
        except ValueError:
            try:
                d = float(y)
            except ValueError:
                d = y
        finally:
            keylist[x] = d
    f.close()
    
    Jx = keylist['Jx']
    Jy = keylist['Jy']
    Jz = keylist['Jz']
    J = [Jx, Jy, Jz]        #Spin-spin coupling
    h = keylist['h']        #External field coupling
    dt = keylist['dt']       #Time step
    d = keylist['d']          #One-particle Hilbert space dim
    chi = keylist['chi']        #Maximum MPS dim      
    N = keylist['N']          #Site number
    T = keylist['T']         #Total time evolved
    step_number = int(T/dt) #Number of time steps
    order = keylist['order']    #Trotter order
    algo = keylist['algo']     #Which algorithm to use
    
    #H belongs to class Hamiltonian which features time evolution and model construction
    H = st.Hamiltonian(J,h,N,dt,d,chi,which= "Heisen",TO = order)
    
    #Energy operators on each bond
    Hchain = H.Hchain
    
    # Build the initial state #
    Psi = st.StateChain(N, chi, d, algo)
        
    #Time evolve the state (imaginary time)
    Psi = H.time_evolve(Psi, step_number, algo)

    #Energy calculation
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
    
    #Energies at each bond
    print("\nBond energies")
    print("Algorithm bond energies: ",E)
    print("\nED trotter bond energies: ", E_ed)
    
    #Build product state from MPS
    # B = Psi.B
    # psi = B[0]
    # for i in range(len(B)-1):
    #     psi = np.tensordot(psi,B[i+1], (i+2,1))
    # 
    # psi = np.reshape(psi, (d**N))
    # print("Algorithm energy from product state: ",sum(ED.get_ener(psi)))
    # print(psi)
    
    #Comparison of MPS and ED trotter GS
    # print("\nNon-zero components of product states")
    # for i in range(len(psi)):
    #     if abs(psi[i])>10**-8 or abs(ED.trot_state[i]) > 10**-8:
    #         
    #         print("Algorithm: ", psi[i]," ED trotter: ",ED.trot_state[i])
    
    #Other output of possible interest
    print("\nOther output")
    print("Lowest energies ",ED.elist)
    ED.exact_time(step_number)
    print("Exact time evolution energy ",ED.E_evol)
    print("Truncation error: ",Psi.err)
    print("Final maximum bond dimension: ",H.chi)
    
    #Correlator testing
    # sz = np.array([[1.,0.],[0.,-1.]])
    # sy = np.array([[0,complex(0,-1)],[complex(0,1),0]])
    # sx = np.array([[0.,1.],[1.,0.]])
    # S = [sx,sy,sz]
    # site1 = 1
    # site2 = 2
    # C=0
    # ED_corr = 0
    # M = st.Measure()
    # for op in S:
    #     C += M.correl(Psi,op,op,site1,site2)
    #     ED_corr += ED.ED_correl(ED.GS[:d**N,0],op,op,site1,site2)
    # print("Spin-spin correlator between site",site1,"and",site2,":",C)
    # print("ED correlator:", ED_corr)
        
    return 0
    
main()
