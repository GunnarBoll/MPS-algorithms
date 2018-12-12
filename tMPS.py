"""tMPS algorithm solving an antiferromagnetic transverse Ising chain
Written by: Gunnar Bollmark
"""
import numpy as np
import storage as st
import imp
import ExactDiag as ed
import time as t
imp.reload(st)
imp.reload(ed)

def algo_output(Psi,H):
    print("\nMPS algorithm results")
    E = Psi.get_ener(H.Hchain)
    print("GS energy:", sum(E))
    print("Bond energies:", E)
    
    if H.N < 12:
        B = Psi.B
        psi = B[0]
        for i in range(len(B)-1):
            psi = np.tensordot(psi,B[i+1], (i+2,1))
        
        psi = np.reshape(psi, (H.d**H.N))
        print("Energy from product state: ",sum(ed.ExactD.get_ener(H,psi)))
    
def exact_output(ED):
    
    print("\nExact diagonalization results")
    print("Exact groundstate energy:", ED.E_GS)
    print("Exact energies:",ED.elist)    

def obser_test(ED, op, N):
    corr_list = []
    maglist = []
    for k in range(len(ED.elist)):
        ED_corr = 0
        ED_mag = []
        for i in range(N):
            ED_mag.append(ED.ED_correl(ED.GSl[:16,k],op,np.eye(2),i,N-1))
            for j in range(N):
                ED_corr += ED.ED_correl(ED.GSl[:16,k],op,op,i,j)
                
        ED_corr = ED_corr - sum(ED_mag)**2
        corr_list.append(ED_corr)
        maglist.append(sum(ED_mag))
    
    print("\nVariance of observable")
    print(corr_list)
    print("Expectation value of observable")
    print(maglist)

### MAIN PROGRAM #######################################
def main():
    model = input("Enter model (Heisen, HCboson): ")
    
    directory = "C:/Users/Gunnar/Documents/Ph.D/Learning/DMRG/Tryout code/models/" + model + ".txt"
    
    f = open(directory,"r")
    inputlist = []
    keylist = {}
    for line in f:
        if line.startswith('#'):
            continue
        d = []
        inputlist.append(line.rstrip('\n').split(" = "))
        x = inputlist[-1][0]
        yl = inputlist[-1][1].split(',')
        for y in yl:
            try:
                d.append(int(y))
            except ValueError:
                try:
                    d.append(float(y))
                except ValueError:
                    d.append(y)
        if len(d) == 1:
            d = d[0]
        keylist[x] = d
    f.close()
    
    
    J = keylist['g']
    h = keylist['U']
    dt = keylist['dt']       #Time step
    d = keylist['d']          #One-particle Hilbert space dim
    chi = keylist['chi']        #Maximum MPS dim      
    N = keylist['N']          #Site number
    T = keylist['T']         #Total time evolved
    step_number = int(T/dt) #Number of time steps
    order = keylist['order']    #Trotter order
    algo = keylist['algo']     #Which algorithm to use
    
    #H belongs to class Hamiltonian which features time evolution and model construction
    H = st.Hamiltonian(J,h,N,dt,d,chi,model,TO = order)
    
    # Build the initial state #
    Psi = st.StateChain(N, chi, d, algo)
        
    #Time evolve the state (imaginary time)
    start = t.process_time()
    Psi = H.time_evolve(Psi, step_number, algo)
    end = t.process_time()
    print("Time taken for algorithm is:", end-start)
    
    #Output data to terminal
    algo_output(Psi,H)

    #Exact diagonalization trotter
    if N <= 20:
        start = t.process_time()
        ED = ed.ExactD(J,h,N,dt,d, model, TO = order)
        end1 = t.process_time()
        ED.exact_GS()
        end2 = t.process_time()
        exact_output(ED)
        print("\nTime taken for full Hamiltonian:",end1-start)
        print("\nTime taken for ED is:", end2-end1)
    
    # if model == "HCboson":
    #     Free = st.Freeferm(J,h,N)
    #     print("\nFree fermion result:",Free.E_GS)
    
    #Other output of possible interest
    print("\nOther output")
    print("Truncation error: ",Psi.err)
    print("Final maximum bond dimension: ",H.chi)
    
    # for l in range(N+1):
    #     print(Psi.L[l])
    #     print(ED.get_lam(ED.trot_state,l))
    
    #Magnetization testing for Heisenberg model
    # sz = np.array([[1.,0.],[0.,-1.]])/2
    # obser_test(ED,sz,N)
    
    #Number operator testing for HCboson model
    # num = np.array([[1,0],[0,0]])
    # obser_test(ED,num,N)
    
    
    
    return 0
        
main()
