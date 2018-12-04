""" ED for short quantum Heisenberg chain
Written by: Gunnar Bollmark"""

import numpy as np
import itertools
import storage
import imp
import scipy.sparse.linalg.eigen.arpack as arp
imp.reload(storage)

class ExactD(storage.Hamiltonian):
    def __init__(self, J, h, N ,dt, d,which = "AFTIC", TO = "second"):
        sz = np.array([[1.,0.],[0.,-1.]])
        sy = np.array([[0,-complex(0,1)],[complex(0,1),0]])
        sx = np.array([[0.,1.],[1.,0.]])
        
        super(ExactD,self).__init__(J,h,N,dt,d, 8, which, TO, ED = True) #Inheritance seems to function strangely??
                
        initspinor = np.array([1,1])/np.sqrt(2)
        self.i_state = Nkron(list(itertools.repeat(initspinor,self.N)))
        
        #Create full Hilbert space Hamiltonian
        self.Id = np.eye(d)
        self.op = list(itertools.repeat(self.Id,N)) #One operator on each site
        self.H = 0
        
        S = np.array([sx, sy, sz])
        
        for i in range(N-1):
            for k in range(len(S)):
                self.op[i] = self.op[i+1] = S[k]
                self.H += J*Nkron(self.op).real
            self.op[i] = sz
            self.op[i+1] = self.Id
            self.H += h*Nkron(self.op)
            self.op[i] = self.op[i+1] = self.Id
        self.op[-1] = sz
        self.H += h*Nkron(self.op)
        self.op[-1] = self.Id
        
        #List of time operators for each bond
        self.Ulist = []
        for i in range(N-1):
            if np.mod(i,2) == 0:
                self.Ulist.append(self.Uodd[i])
            elif np.mod(i,2) == 1:
                self.Ulist.append(self.Ueven[i])
        
        return
        
    def trotter_time(self,step_number):
        state = self.i_state
        Iden = self.op[:self.N-1]
        
        if self.TO == "first":
            for t in range(step_number):
                state = self.ED_sweep(self.Ulist, Iden, state, self.Id, forw = True)
                state = self.ED_sweep(self.Ulist, Iden, state, self.Id, forw = False)
        
        elif self.TO == "second":
            Utfull = []
            for i in range(0,self.N-1,2):
                Utfull.append( np.matmul(self.Ulist[i],self.Ulist[i]))
                if i != self.N-2:
                    Utfull.append(self.Ulist[i+1])
            
            state = self.ED_sweep(self.Ulist, Iden, state, self.Id, forw = True)
            
            for t in range(step_number-1):
                state = self.ED_sweep(Utfull, Iden, state, self.Id, forw = False)
                state = self.ED_sweep(Utfull, Iden, state, self.Id, forw = True)
                
            state = self.ED_sweep(Utfull, Iden, state, self.Id, forw = False)
            
            state = self.ED_sweep(self.Ulist, Iden, state, self.Id, forw = True)
            
        self.trot_state = state
        
        return
    
    def exact_time(self, step_number):
        state = self.i_state
        e, v = np.linalg.eig(self.H)
        
        expH = np.matmul(v, np.matmul(np.diag(np.exp(-self.dt*e)),np.linalg.inv(v)))
        
        for t in range(step_number):
            state = np.dot(expH,state)
            state = state/np.sqrt(sum(state**2))
        
        self.E_evol = np.dot(state,np.dot(self.H,state))
        self.evol_state = state
        
        return
        
    def exact_GS(self):
        self.elist, self.GS = arp.eigsh(self.H, k=6, v0 = self.i_state, which = "SA")
        self.E_GS = min(self.elist)
        return self.E_GS, self.GS

    def get_ener(self,state):
        E = []
        Iden = self.op[:self.N-1]
        for k in range(self.N-1):
            Iden[k] = self.Hchain[k]
            Enop = Nkron(Iden)
            E.append(np.dot(state,np.dot(Enop,state)))
            Iden[k] = self.Id
        
        return E


    def ED_sweep(self, Ut, Iden, Istate, I, forw = True):
        start = 0
        end = self.N-1
        inc = 2
        if not forw:
            start = self.N-3
            end = 0
            inc = -2
        for i in range(start,end,inc):
            Iden[i] = Ut[i]
            U = Nkron(Iden)
            Istate = np.dot(U,Istate)
            Istate = Istate/np.sqrt(sum(Istate**2))
            Iden[i] = I
        return Istate
        
    def ED_correl(self,psi,op1,op2,i,j):
        self.op[i] = op1
        op1 = Nkron(self.op)
        self.op[i] = self.Id
        self.op[j] = op2
        op2 = Nkron(self.op)
        self.op[j] = self.Id
        
        phi = np.dot(op2,np.dot(op1,psi))
        
        corr = np.dot(psi,phi)
        
        
        return corr

#Kronecker product of N matrices
def Nkron(matlist):
    if matlist == []:
        return 1
    else:
        return np.kron(matlist[0],Nkron(matlist[1:]))
