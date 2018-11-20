"""tMPS algorithm solving an antiferromagnetic transverse Ising chain
Written by: Gunnar Bollmark
"""

import numpy as np
import storage
import imp
imp.reload(storage)

#Function for getting the energy of a state in gamma/lambda notation        
def get_ener(Lam, B, Hchain, N):
    E = []
    for b in range(0,N-1):
        LB = np.tensordot(np.diag(Lam[b]),B[b],(1,1))
        LBL = np.tensordot(LB,np.diag(Lam[b+1]), (2,0))
        LBLB = np.tensordot(LBL,B[b+1], (2,1))
        LBLBL = np.tensordot(LBLB, np.diag(Lam[b+2]), (3,0))
        C = np.tensordot(LBLBL, np.reshape(Hchain[b], (2,2,2,2)), ([1,2],[2,3]))
        E.append(np.tensordot(np.conj(LBLBL), C, ([0,3,1,2],[0,1,2,3])))
    return E

#Assuming completely right-normalized form, changes the state to Lambda/Gamma notation
def LamGam_form(B, N, chi, d):
    Lam = []
    Gam = []
    for n in range(N+1):
        Lam.append(np.ones([1]))
    
    for i in range(N):
        chia = B[i].shape[1]
        chib = B[i+1].shape[2]
        theta = np.tensordot(np.diag(Lam[i]), B[i], (1,1))
        theta = np.tensordot(theta, B[i+1], (2,1))
        theta = np.reshape(np.transpose(theta, (1,0,2,3)), (d*chia, d*chib))
        
        U, S, V = np.linalg.svd(theta)
        V = V.T
        
        chic = min([np.sum(S>10**-12), chi])
        
        U = np.reshape(U[:d*chia,:chic], (d,chia, chic))
        V = np.transpose(np.reshape(V[:d*chib,:chic], (d,chib, chic)), (0,2,1))
        
        Gam.append(np.transpose(np.tensordot(np.diag(Lam[i]**-1), U, (1,1)), (1,0,2)))
        Lam[i+1] = S[:chic]/np.sqrt(np.sum(S[:chic]**2))
        B[i] = U
        B[i+1] = V
        
    Gam.append(np.array([1,0]))
    return Lam, Gam
    
def LamGam_test(A, N, chi, d):
    Lam = []
    Gam = []
    for n in range(N+1):
        Gam.append(np.ones([1]))
        Lam.append(np.ones([1]))
    
    for i in range(N-2,-1,-1):
        chia = A[i].shape[1]
        chib = A[i+1].shape[2]
        theta = np.tensordot(A[i], A[i+1], (2,1))
        theta = np.tensordot(theta, np.diag(Lam[i+2]), (3,0))
        theta = np.reshape(theta, (d*chia, d*chib))
        
        U, S, V = np.linalg.svd(theta)
        V = V.T
        
        chic = min([np.sum(S>10**-12), chi])
        
        U = np.reshape(U[:d*chia,:chic], (d,chia, chic))
        V = np.transpose(np.reshape(V[:d*chib,:chic], (d,chib, chic)), (0,2,1))
        
        Gam[i+1] = np.tensordot(V, np.diag(Lam[i+2]**-1), (2,0))
        Lam[i+1] = S[:chic]/np.sqrt(np.sum(S[:chic]**2))
        A[i] = U
        A[i+1] = V
    Gam[0] = A[0]
    
    return Lam, Gam

### MAIN PROGRAM #######################################
def main():
    #Model parameters
    J = 1.0        #Spin-spin coupling
    h = 1.0        #External field coupling
    d = 2          #One-particle Hilbert space dim
    chi = 8        #Maximum MPS dim
    dt = 0.01      #Time step
    N = 4          #Site number
    step_number = int(6.5/dt) #Number of time steps
    order = "second"    #Trotter order
    
    #H belongs to class Hamiltonian which features time evolution and model construction
    H = storage.Hamiltonian(J,h,N,dt,d,chi,which= "AFTIC",TO = order)
    
    #Energy operators on each bond
    Hchain = H.Hchain
    
    # Build the initial state #
    B = []
    Lam = []
    for ind in range(N+1):
        B.append(np.zeros([d, 1, 1]))
        
        #All spins pointing up
        B[-1][0,0,0] = 1
        
        #Lambda for tebd
        # Lam.append(np.zeros(1))
        # Lam[-1][0] = 1
        
    Lam = np.array([1])
        
    #Time evolve the state
    B, Lam = H.time_evolve(B,Lam,step_number, "tDMRG")
    
    #Attempts at energy calculation
    
    #Using a two-site Hamiltonian at each bond and find total energy per site
    
    #Brings state from lambda/gamma to B/A notation
    # for i in range(N):
    #     B[i] = np.tensordot(B[i], np.diag(Lam[i+1]), (2,0))
    # for i in range(N):
    #     B[i] = np.transpose(np.tensordot(np.diag(Lam[i]), B[i], (1,1)), (1,0,2))
    
    #Brings state from A/B notation to lambda/gamma notation
    Lam, B = LamGam_test(B, N, chi, d)
    
    E = get_ener(Lam, B, Hchain, N)
        
    print(sum(E))
    

    return 0
    
main()
