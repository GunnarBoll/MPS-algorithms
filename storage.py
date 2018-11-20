

import numpy as np

#Class for Hamiltonian
class Hamiltonian:
    def __init__(self,J,h,N,dt,d,chi,which = "AFTIC", TO = "first"):
        self.d = d
        self.chi = chi
        self.N = N
        self.TO = TO
        
        if which == "AFTIC":
            Hlist = get_AFTIC(J,h)
            self.Hchain = ener_chain(Hlist, N)
            self.I = np.transpose(np.reshape(np.kron(np.eye(2),np.eye(2)), (2,2,2,2)), (0,2,1,3))
            
            #If first order do not put a half factor on odd operators
            if self.TO == "first":
                plist = [1]
                FO = True 
            
            elif self.TO == "second":
                plist = [1]
                FO = False
            
            #Construct even and odd bond time evolution operators
            for p in plist:
                self.Uodd, self.Ueven = model_constructor(Hlist, p, self.I, N, dt, d, FO)
            
            return
        else:
            return
            
    def time_evolve(self, B, Lam, step_number, algo):
        #First order algorithm
        if self.TO == "first":
            
            sweep_order = [self.Uodd, self.Ueven]
            direc = True #True for forward sweep
            
            for t in range(step_number):
                for oper in sweep_order:
                    B, Lam = sweep(B, Lam, self.I, oper, self.chi, self.d, algo, forward = direc)
                    direc = not direc
                
        
        #Second order algorithm
        elif self.TO == "second":
            #Time evolution operators for full time steps
            Utful = []
            for n in range(0,self.N-1,2):
                Utful.append(np.transpose( np.tensordot(self.Uodd[n], self.Uodd[n], ([1,3],[0,2]) ), (0,2,1,3)))
                if n != self.N-2:
                    Utful.append(self.Uodd[n+1])
            
            # One half time step on odd bonds
            B, Lam = sweep(B, Lam, self.I, self.Uodd, self.chi, self.d, algo, forward = True)
            
            # Full loop for time evolution
            for t in range(step_number-1): 
                B, Lam = sweep(B, Lam, self.I, self.Ueven, self.chi, self.d, algo, forward = False)
                B, Lam = sweep(B, Lam, self.I, Utful, self.chi, self.d, algo, forward = True)
                    
            #Final even bond step
            B, Lam = sweep(B, Lam, self.I, self.Ueven, self.chi, self.d, algo, forward = False)
            
            #One final half step on odd bonds
            B, Lam = sweep(B, Lam, self.I, self.Uodd, self.chi, self.d, algo, forward = True)
        
        return B, Lam




#Creates two-site Hamiltonians for an Anti-Ferromagnetic Transverse Ising Chain (AFTIC)
def get_AFTIC(J,h):
    Hlist = []
    parlist = [[1.0,1.0,0.5],[1.0,0.5,0.5],[1.0,0.5,1.0]]
    for params in parlist:
        Hlist.append(kron_AFTIC(params[0],params[1],params[2]))
    
    return Hlist

#Construction of a single two-site Hamiltonian
def kron_AFTIC(J,h1,h2):
    sz = np.array([[1.,0.],[0.,-1.]]) 
    sx = np.array([[0.,1.],[1.,0.]])
    
    H = J*np.kron(sz,sz) + h1*np.kron(sx,np.eye(2)) + h2*np.kron(np.eye(2),sx)
    return H

#Constructs the bond Hamiltonians for energy calculation
def ener_chain(Hlist,N):
    Hchain = []
    Hchain.append(Hlist[0])
    for i in range(N-3):
        Hchain.append(Hlist[1])
    Hchain.append(Hlist[-1])
    return Hchain

#Creates the time evolution for even and odd sites
def model_constructor(Hlist, p,I,N,dt,d,FO = False):
    podd = peven = p
    if not FO:
        podd = podd/2
    
    time_ops = []
    for ps, odd in [[podd,peven],[True, False]]:
        time_ops.append(timeop_chain(Hlist,I,N,dt,d,ps,odd))
    return time_ops
    
#Constructs a chain of time evolution operators (even or odd spaces)
def timeop_chain(Hlist, I, N, dt, d, p, odd = True):
    Ulist = []
    expHlist = []
    for H in Hlist:
        expHlist.append(get_timeop(H,dt,d,p))
    
    iter_list = [I, expHlist[1]]
    i = 0
    if not odd:
        expHlist[0] = expHlist[-1] = I
        iter_list[0], iter_list[1] = iter_list[1], iter_list[0]
    
    Ulist.append(expHlist[0])
    while len(Ulist)<N-2:
        Ulist.append(iter_list[np.mod(i,2)])
        i+=1
    Ulist.append(expHlist[-1])
    
    return Ulist

#Constructs a time evolution operator given a Hamiltonian
def get_timeop(H, dt, d, p = 1):
    e, v = np.linalg.eigh(H)
    Ut = np.matmul(v,np.matmul(np.diag(np.exp(-dt*p*e)), np.linalg.inv(v)))
    Ut = np.transpose(np.reshape(Ut, (d,d,d,d)), (0,2,1,3))
    return Ut

#Applies time evolution to the chain (sweeping forward or backward)
def sweep(B, Lam, I, time_ops, chi, d, algo, forward = True):
    sites = range(len(time_ops))
    if not forward:
        sites = range(len(time_ops)-1, -1, -1)
    
    for i in sites:
        #Change tebd function for other algorithm
        if algo == "TEBD":
            B, Lam = tebd(B, Lam, time_ops, I, i, chi, d)
        elif algo == "tDMRG":
            B, Lam = tdmrg(B, Lam, time_ops, I, i, chi, d, forward)
    
    if algo == "tDMRG":
        if forward:
            B[len(time_ops)] = np.transpose(np.tensordot(np.diag(Lam), B[len(time_ops)], (1,1)), (1,0,2))
            Lam = np.array([1])
        else:
            B[0] = np.tensordot( B[0],np.diag(Lam), (2,0))
            Lam = np.array([1])
    return B, Lam

#Main TEBD algorithm (utilizing inverse lambda matrices)
def tebd(B,Lam,Ut, I, i, chi, d):
    chia = B[i].shape[1]
    chib = B[i+1].shape[2]
    
    theta = get_theta(Lam, B, i)
    
    if np.all(Ut[i] != I):
        phi = np.transpose(np.tensordot(Ut[i], theta, ([1,3], [0,2])), (0,2,1,3))
    elif np.all(Ut[i] == I):
        return B, Lam
    
    phi = np.reshape(phi, (d*chia,d*chib))
    
    #Singular value decomposition    
    U, S, V, chic = svd_truncator(phi, chia, chib, chi, d)
    
    #State update
    B[i+1] = np.tensordot(V, np.diag(Lam[i+2]**(-1)), (2,0))
    B[i] = np.transpose(np.tensordot(np.diag(Lam[i]**(-1)),U, (1,1)), (1,0,2))
    Lam[i+1] = S/np.sqrt(sum(S**2))
    
    
    return B, Lam

#Constructs the theta tensor (Psi tensor) for TEBD algorithm
def get_theta(Lam, B, i):
    theta = np.tensordot(np.diag(Lam[i]),B[i], (1,1))
    theta = np.tensordot(theta, np.diag(Lam[i+1]), (2,0))
    theta = np.tensordot(theta, B[i+1], (2,1))
    theta = np.tensordot(theta, np.diag(Lam[i+2]), (3,0))
    theta = np.transpose(theta, (1,0,2,3))
    return theta

#tDMRG algorithm for time evolution
def tdmrg(B, L, Ut, I, i, chi, d, forward):
    chia = B[i].shape[1]
    chib = B[i+1].shape[2]
    
    
    if forward:
        theta = np.tensordot(np.diag(L), B[i], (1,1) )
        theta = np.tensordot(theta, B[i+1], (2,1))
        theta = np.transpose(theta, (1,0,2,3))
    else:
        theta = np.tensordot(B[i],B[i+1], (2,1))
        theta = np.tensordot(theta, np.diag(L), (3,0))
    
    
    phi = np.transpose(np.tensordot(Ut[i], theta, ([1,3],[0,2])), (0,2,1,3))
    
    phi = np.reshape(phi,(d*chia,d*chib))
    
    U, S, V, chic = svd_truncator(phi, chia, chib, chi ,d)
    
    B[i] = U
    B[i+1] = V
    L = S/np.sqrt(sum(S**2))
    
    return B, L

#Performs an SVD and truncates the singular values to specified bond dimension
def svd_truncator(phi, chia, chib, chi, d):
    
    U, S, V = np.linalg.svd(phi)
    V = V.T
    
    chic = min([np.sum(S>10**-12), chi])
    
    U = np.reshape(U[:d*chia,:chic], (d,chia,chic))
    V = np.transpose(np.reshape(V[:d*chib, :chic], (d,chib,chic)), (0,2,1))
    S = S[:chic]
    
    return [U, S ,V, chic]



    




