""" ED for short antiferromagnetic Ising chain 
with transversal external magnetic field
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
        self.I2 = np.eye(d)
        self.op = list(itertools.repeat(self.I2,N)) #One operator on each site
        self.H = 0
        
        S = np.array([sx, sy, sz])
        
        for i in range(N-1):
            for k in range(len(S)):
                self.op[i] = self.op[i+1] = S[k]
                self.H += J*Nkron(self.op).real
            self.op[i] = sz
            self.op[i+1] = self.I2
            self.H += h*Nkron(self.op)
            self.op[i] = self.op[i+1] = self.I2
        self.op[-1] = sz
        self.H += h*Nkron(self.op)
        self.op[-1] = self.I2
        
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
                state = self.ED_sweep(self.Ulist, Iden, state, self.I2, forw = True)
                state = self.ED_sweep(self.Ulist, Iden, state, self.I2, forw = False)
        
        elif self.TO == "second":
            Utfull = []
            for i in range(0,self.N-1,2):
                Utfull.append( np.matmul(self.Ulist[i],self.Ulist[i]))
                if i != self.N-2:
                    Utfull.append(self.Ulist[i+1])
            
            state = self.ED_sweep(self.Ulist, Iden, state, self.I2, forw = True)
            
            for t in range(step_number-1):
                state = self.ED_sweep(Utfull, Iden, state, self.I2, forw = False)
                state = self.ED_sweep(Utfull, Iden, state, self.I2, forw = True)
                
            state = self.ED_sweep(Utfull, Iden, state, self.I2, forw = False)
            
            state = self.ED_sweep(self.Ulist, Iden, state, self.I2, forw = True)
            
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
            Iden[k] = self.I2
        
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

#Kronecker product of N matrices
def Nkron(matlist):
    if matlist == []:
        return 1
    else:
        return np.kron(matlist[0],Nkron(matlist[1:]))

# def get_4OTS():
#     Hfirst = J*np.kron(sz,sz) + h*np.kron(sx,np.eye(2)) + h/2*np.kron(np.eye(2),sx)
#     Hbond = J*np.kron(sz,sz) + h/2*np.kron(sx,np.eye(2)) + h/2*np.kron(np.eye(2),sx)
#     Hlast = J*np.kron(sz,sz) + h/2*np.kron(sx,np.eye(2)) + h*np.kron(np.eye(2),sx)
#     
#     vec = np.array([1,1])/np.sqrt(2)
#     
#     Istate = Nkron(list(itertools.repeat(vec,N)))
#     
#     
#     Iden = list(itertools.repeat(I,N-1))
#     
#     ef, vf = np.linalg.eig(Hfirst)
#     eb, vb = np.linalg.eig(Hbond)
#     el, vl = np.linalg.eig(Hlast)
#     
#     p = 1/(4-4**(1/3))
#     expHf1 = np.matmul(vf,np.matmul(np.diag(np.exp(-dt*p*ef/2)), np.linalg.inv(vf)))
#     expHf2 = np.matmul(vf,np.matmul(np.diag(np.exp(-dt*(1-4*p)*ef/2)), np.linalg.inv(vf)))
#     expHodd1 = np.matmul(vb,np.matmul(np.diag(np.exp(-dt*p*eb/2)), np.linalg.inv(vb)))
#     expHodd2 = np.matmul(vb,np.matmul(np.diag(np.exp(-dt*(1-4*p)*eb/2)), np.linalg.inv(vb)))
#     expHeven1 = np.matmul(vb,np.matmul(np.diag(np.exp(-dt*p*eb)), np.linalg.inv(vb)))
#     expHeven2 = np.matmul(vb,np.matmul(np.diag(np.exp(-dt*(1-4*p)*eb)), np.linalg.inv(vb)))
#     expHl1 = np.matmul(vl,np.matmul(np.diag(np.exp(-dt*p*el)), np.linalg.inv(vl)))
#     expHl2 = np.matmul(vl,np.matmul(np.diag(np.exp(-dt*(1-4*p)*el)), np.linalg.inv(vl)))
#     
#     Ulist1 = []
#     Ulist2 = []
#     Hchain = []
#     
#     for n in range(N-1):
#         if n == 0:
#             Ulist1.append(expHf1)
#             Ulist2.append(expHf2)
#             Hchain.append(Hfirst)
#         elif n == N-2:
#             Ulist1.append(expHl1)
#             Ulist2.append(expHl2)
#             Hchain.append(Hlast)
#         elif np.mod(n,2) == 0:
#             Ulist1.append(expHodd1)
#             Ulist2.append(expHodd2)
#             Hchain.append(Hbond)
#         else:
#             Ulist1.append(expHeven1)
#             Ulist2.append(expHeven2)
#             Hchain.append(Hbond)
#             
#             
#     Utfull1 = []
#     Utfull2 = []
#     
#     for i in range(0,N-1,2):
#         Utfull1.append( np.matmul(Ulist1[i],Ulist1[i]))
#         Utfull2.append( np.matmul(Ulist2[i],Ulist2[i]))
#         if i != N-2:
#             Utfull1.append(Ulist1[i+1])
#             Utfull2.append(Ulist2[i+1])
#     
#     for t in range(step_number):
#         for ind in range(2):
#             for i in range(0,N-1,2):
#                 Iden[i] = Ulist1[i]
#                 U = Nkron(Iden)
#                 Istate = np.dot(U,Istate)
#                 Istate = Istate/np.sqrt(sum(Istate**2))
#                 Iden[i] = I
#             for j in range(N-3,0,2):
#                 Iden[i] = Ulist1[i]
#                 U = Nkron(Iden)
#                 Istate = np.dot(U,Istate)
#                 Istate = Istate/np.sqrt(sum(Istate**2))
#                 Iden[i] = I
#             for i in range(0,N-1,2):
#                 Iden[i] = Ulist1[i]
#                 U = Nkron(Iden)
#                 Istate = np.dot(U,Istate)
#                 Istate = Istate/np.sqrt(sum(Istate**2))
#                 Iden[i] = I
#                 
#         for i in range(0,N-1,2):
#             Iden[i] = Ulist2[i]
#             U = Nkron(Iden)
#             Istate = np.dot(U,Istate)
#             Istate = Istate/np.sqrt(sum(Istate**2))
#             Iden[i] = I
#         for j in range(N-3,0,2):
#             Iden[i] = Ulist2[i]
#             U = Nkron(Iden)
#             Istate = np.dot(U,Istate)
#             Istate = Istate/np.sqrt(sum(Istate**2))
#             Iden[i] = I
#         for i in range(0,N-1,2):
#             Iden[i] = Ulist2[i]
#             U = Nkron(Iden)
#             Istate = np.dot(U,Istate)
#             Istate = Istate/np.sqrt(sum(Istate**2))
#             Iden[i] = I
#             
#         for ind in range(2):
#             for i in range(0,N-1,2):
#                 Iden[i] = Ulist1[i]
#                 U = Nkron(Iden)
#                 Istate = np.dot(U,Istate)
#                 Istate = Istate/np.sqrt(sum(Istate**2))
#                 Iden[i] = I
#             for j in range(N-3,0,2):
#                 Iden[i] = Ulist1[i]
#                 U = Nkron(Iden)
#                 Istate = np.dot(U,Istate)
#                 Istate = Istate/np.sqrt(sum(Istate**2))
#                 Iden[i] = I
#             for i in range(0,N-1,2):
#                 Iden[i] = Ulist1[i]
#                 U = Nkron(Iden)
#                 Istate = np.dot(U,Istate)
#                 Istate = Istate/np.sqrt(sum(Istate**2))
#                 Iden[i] = I
#     
#     
#     return Istate, Hchain

"""
J = 1.0 #Spin-spin coupling
h = 10.0 #External field coupling
N = 4  #Number of sites
d = 2
dt = 0.1
step_number = int(10/dt)





ED = ExactD(J,h,N,dt,d,which = "AFTIC", TO = "second")

ED.trotter_time(step_number)
E = ED.get_ener(ED.trot_state)
print("Energy with Trotter:", sum(E))

if N <10:
    ED.exact_time(step_number)
    print("Exact time evolution:", ED.E_evol)

ED.exact_GS()
print("GS energy:", ED.E_GS)


print(E)
"""