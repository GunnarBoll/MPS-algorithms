""" ED for short antiferromagnetic Ising chain 
with transversal external magnetic field
Written by: Gunnar Bollmark"""

import numpy as np
import itertools

J = 1.0 #Spin-spin coupling
h = 1.0 #External field coupling
N = 4  #Number of sites
dt = 0.01
step_number = int(6.5/dt)

#Kronecker product of N matrices
def Nkron(matlist):
    if matlist == []:
        return 1
    else:
        return np.kron(matlist[0],Nkron(matlist[1:]))



# Construct Hamiltonian
sz = np.array([[1.,0.],[0.,-1.]])
sx = np.array([[0.,1.],[1.,0.]])
I = np.eye(2)

def get_FOTS():
    Hfirst = J*np.kron(sz,sz) + h*np.kron(sx,np.eye(2)) + h/2*np.kron(np.eye(2),sx)
    Hbond = J*np.kron(sz,sz) + h/2*np.kron(sx,np.eye(2)) + h/2*np.kron(np.eye(2),sx)
    Hlast = J*np.kron(sz,sz) + h/2*np.kron(sx,np.eye(2)) + h*np.kron(np.eye(2),sx)
    
    vec = np.array([1,0])
    
    Istate = Nkron(list(itertools.repeat(vec,N)))
    
    
    Iden = list(itertools.repeat(I,N-1))
    
    ef, vf = np.linalg.eig(Hfirst)
    eb, vb = np.linalg.eig(Hbond)
    el, vl = np.linalg.eig(Hlast)
    
    expHf = np.matmul(vf,np.matmul(np.diag(np.exp(-dt*ef)), np.linalg.inv(vf)))
    expHb = np.matmul(vb,np.matmul(np.diag(np.exp(-dt*eb)), np.linalg.inv(vb)))
    expHl = np.matmul(vl,np.matmul(np.diag(np.exp(-dt*el)), np.linalg.inv(vl)))
    
    Ulist = []
    Hchain = []
    
    for n in range(0,N-1):
        if n == 0:
            Ulist.append(expHf)
            Hchain.append(Hfirst)
        elif n == N-2:
            Ulist.append(expHl)
            Hchain.append(Hlast)
        else:
            Ulist.append(expHb)
            Hchain.append(Hbond)
    
    for t in range(step_number):
        for i in range(0,N-1,2):
            Iden[i] = Ulist[i]
            U = Nkron(Iden)
            Istate = np.dot(U,Istate)
            Istate = Istate/np.sqrt(sum(Istate**2))
            Iden[i] = I
            
        for i in range(N-3,0,-2):
            Iden[i] = Ulist[i]
            U = Nkron(Iden)
            Istate = np.dot(U,Istate)
            Istate = Istate/np.sqrt(sum(Istate**2))
            Iden[i] = I
    return Istate, Hchain
    
    
def get_SOTS():
    Hfirst = J*np.kron(sz,sz) + h*np.kron(sx,np.eye(2)) + h/2*np.kron(np.eye(2),sx)
    Hbond = J*np.kron(sz,sz) + h/2*np.kron(sx,np.eye(2)) + h/2*np.kron(np.eye(2),sx)
    Hlast = J*np.kron(sz,sz) + h/2*np.kron(sx,np.eye(2)) + h*np.kron(np.eye(2),sx)
    
    vec = np.array([1,0])
    
    Istate = Nkron(list(itertools.repeat(vec,N)))
    
    
    Iden = list(itertools.repeat(I,N-1))
    
    ef, vf = np.linalg.eig(Hfirst)
    eb, vb = np.linalg.eig(Hbond)
    el, vl = np.linalg.eig(Hlast)
    
    expHf = np.matmul(vf,np.matmul(np.diag(np.exp(-dt*ef/2)), np.linalg.inv(vf)))
    expHodd = np.matmul(vb,np.matmul(np.diag(np.exp(-dt*eb/2)), np.linalg.inv(vb)))
    expHeven = np.matmul(vb,np.matmul(np.diag(np.exp(-dt*eb)), np.linalg.inv(vb)))
    expHl = np.matmul(vl,np.matmul(np.diag(np.exp(-dt*el/2)), np.linalg.inv(vl)))
    
    Ulist = []
    Hchain = []
    
    for n in range(N-1):
        if n == 0:
            Ulist.append(expHf)
            Hchain.append(Hfirst)
        elif n == N-2:
            Ulist.append(expHl)
            Hchain.append(Hlast)
        elif np.mod(n,2) == 0:
            Ulist.append(expHodd)
            Hchain.append(Hbond)
        else:
            Ulist.append(expHeven)
            Hchain.append(Hbond)
    
    Utfull = []
    for i in range(0,N-1,2):
        Utfull.append( np.matmul(Ulist[i],Ulist[i]))
        if i != N-2:
            Utfull.append(Ulist[i+1])
    
    for i in range(0,N-1,2):
        Iden[i] = Ulist[i]
        U = Nkron(Iden)
        Istate = np.dot(U,Istate)
        Istate = Istate/np.sqrt(sum(Istate**2))
        Iden[i] = I
    
    for t in range(step_number-1):
        for i in range(N-3,0,-2):
            Iden[i] = Utfull[i]
            U = Nkron(Iden)
            Istate = np.dot(U,Istate)
            Istate = Istate/np.sqrt(sum(Istate**2))
            Iden[i] = I
        for i in range(0,N-1,2):
            Iden[i] = Utfull[i]
            U = Nkron(Iden)
            Istate = np.dot(U,Istate)
            Istate = Istate/np.sqrt(sum(Istate**2))
            Iden[i] = I
    
    for i in range(N-3,0,-2):
        Iden[i] = Utfull[i]
        U = Nkron(Iden)
        Istate = np.dot(U,Istate)
        Istate = Istate/np.sqrt(sum(Istate**2))
        Iden[i] = I
    
    for i in range(0,N-1,2):
        Iden[i] = Ulist[i]
        U = Nkron(Iden)
        Istate = np.dot(U,Istate)
        Istate = Istate/np.sqrt(sum(Istate**2))
        Iden[i] = I
        
        
    
    return Istate, Hchain
    
def get_4OTS():
    Hfirst = J*np.kron(sz,sz) + h*np.kron(sx,np.eye(2)) + h/2*np.kron(np.eye(2),sx)
    Hbond = J*np.kron(sz,sz) + h/2*np.kron(sx,np.eye(2)) + h/2*np.kron(np.eye(2),sx)
    Hlast = J*np.kron(sz,sz) + h/2*np.kron(sx,np.eye(2)) + h*np.kron(np.eye(2),sx)
    
    vec = np.array([1,0])
    
    Istate = Nkron(list(itertools.repeat(vec,N)))
    
    
    Iden = list(itertools.repeat(I,N-1))
    
    ef, vf = np.linalg.eig(Hfirst)
    eb, vb = np.linalg.eig(Hbond)
    el, vl = np.linalg.eig(Hlast)
    
    p = 1/(4-4**(1/3))
    expHf1 = np.matmul(vf,np.matmul(np.diag(np.exp(-dt*p*ef/2)), np.linalg.inv(vf)))
    expHf2 = np.matmul(vf,np.matmul(np.diag(np.exp(-dt*(1-4*p)*ef/2)), np.linalg.inv(vf)))
    expHodd1 = np.matmul(vb,np.matmul(np.diag(np.exp(-dt*p*eb/2)), np.linalg.inv(vb)))
    expHodd2 = np.matmul(vb,np.matmul(np.diag(np.exp(-dt*(1-4*p)*eb/2)), np.linalg.inv(vb)))
    expHeven1 = np.matmul(vb,np.matmul(np.diag(np.exp(-dt*p*eb)), np.linalg.inv(vb)))
    expHeven2 = np.matmul(vb,np.matmul(np.diag(np.exp(-dt*(1-4*p)*eb)), np.linalg.inv(vb)))
    expHl1 = np.matmul(vl,np.matmul(np.diag(np.exp(-dt*p*el)), np.linalg.inv(vl)))
    expHl2 = np.matmul(vl,np.matmul(np.diag(np.exp(-dt*(1-4*p)*el)), np.linalg.inv(vl)))
    
    Ulist1 = []
    Ulist2 = []
    Hchain = []
    
    for n in range(N-1):
        if n == 0:
            Ulist1.append(expHf1)
            Ulist2.append(expHf2)
            Hchain.append(Hfirst)
        elif n == N-2:
            Ulist1.append(expHl1)
            Ulist2.append(expHl2)
            Hchain.append(Hlast)
        elif np.mod(n,2) == 0:
            Ulist1.append(expHodd1)
            Ulist2.append(expHodd2)
            Hchain.append(Hbond)
        else:
            Ulist1.append(expHeven1)
            Ulist2.append(expHeven2)
            Hchain.append(Hbond)
            
            
    Utfull1 = []
    Utfull2 = []
    
    for i in range(0,N-1,2):
        Utfull1.append( np.matmul(Ulist1[i],Ulist1[i]))
        Utfull2.append( np.matmul(Ulist2[i],Ulist2[i]))
        if i != N-2:
            Utfull1.append(Ulist1[i+1])
            Utfull2.append(Ulist2[i+1])
    
    for t in range(step_number):
        for ind in range(2):
            for i in range(0,N-1,2):
                Iden[i] = Ulist1[i]
                U = Nkron(Iden)
                Istate = np.dot(U,Istate)
                Istate = Istate/np.sqrt(sum(Istate**2))
                Iden[i] = I
            for j in range(N-3,0,2):
                Iden[i] = Ulist1[i]
                U = Nkron(Iden)
                Istate = np.dot(U,Istate)
                Istate = Istate/np.sqrt(sum(Istate**2))
                Iden[i] = I
            for i in range(0,N-1,2):
                Iden[i] = Ulist1[i]
                U = Nkron(Iden)
                Istate = np.dot(U,Istate)
                Istate = Istate/np.sqrt(sum(Istate**2))
                Iden[i] = I
                
        for i in range(0,N-1,2):
            Iden[i] = Ulist2[i]
            U = Nkron(Iden)
            Istate = np.dot(U,Istate)
            Istate = Istate/np.sqrt(sum(Istate**2))
            Iden[i] = I
        for j in range(N-3,0,2):
            Iden[i] = Ulist2[i]
            U = Nkron(Iden)
            Istate = np.dot(U,Istate)
            Istate = Istate/np.sqrt(sum(Istate**2))
            Iden[i] = I
        for i in range(0,N-1,2):
            Iden[i] = Ulist2[i]
            U = Nkron(Iden)
            Istate = np.dot(U,Istate)
            Istate = Istate/np.sqrt(sum(Istate**2))
            Iden[i] = I
            
        for ind in range(2):
            for i in range(0,N-1,2):
                Iden[i] = Ulist1[i]
                U = Nkron(Iden)
                Istate = np.dot(U,Istate)
                Istate = Istate/np.sqrt(sum(Istate**2))
                Iden[i] = I
            for j in range(N-3,0,2):
                Iden[i] = Ulist1[i]
                U = Nkron(Iden)
                Istate = np.dot(U,Istate)
                Istate = Istate/np.sqrt(sum(Istate**2))
                Iden[i] = I
            for i in range(0,N-1,2):
                Iden[i] = Ulist1[i]
                U = Nkron(Iden)
                Istate = np.dot(U,Istate)
                Istate = Istate/np.sqrt(sum(Istate**2))
                Iden[i] = I
    
    
    return Istate, Hchain

Istate, Hchain = get_SOTS()
Iden = list(itertools.repeat(I,N-1))

E = []
for k in range(N-1):
    Iden[k] = Hchain[k]
    Enop = Nkron(Iden)
    E.append(np.dot(Istate,np.dot(Enop,Istate)))
    Iden[k] = I



op = list(itertools.repeat(I,N))
vec = np.array([1,0])
H = 0
Istate2 = Nkron(list(itertools.repeat(vec,N)))

for i in range(N-1):
    op[i] = op[i+1] = sz
    H += J*Nkron(op)
    op[i] = sx
    op[i+1] = I
    H += h*Nkron(op)
    op[i] = op[i+1] = I

op[-1] = sx
H += h*Nkron(op)

e, v = np.linalg.eig(H)

expH = np.matmul(v, np.matmul(np.diag(np.exp(-dt*e)),np.linalg.inv(v)))

for t in range(step_number):
    Istate2 = np.dot(expH,Istate2)
    Istate2 = Istate2/np.sqrt(sum(Istate2**2))

ener = np.dot(Istate2,np.dot(H,Istate2))
print("Energy with Trotter:", sum(E))
print("Energy of time evolution: ", ener)
print("GS energy = ", min(e))
