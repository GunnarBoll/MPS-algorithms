import numpy as np
import scipy.linalg as scl
import itertools
import ED_AFTIC as ed

def tebd(G,L,Ut,i, forw):
    chia = G[i].shape[1]
    chib = G[i+1].shape[2]
    
    theta = np.tensordot(np.diag(L[i]), G[i], (1,1))
    theta = np.tensordot(theta, np.diag(L[i+1]), (2,0))
    theta= np.tensordot(theta, G[i+1], (2,1))
    theta = np.transpose(np.tensordot(theta, np.diag(L[i+2]), (3,0)), (1,0,2,3))
    
    # if not forw:
    #     theta = np.tensordot(G[i],G[i+1],(2,1))
    #     theta = np.tensordot(theta,np.diag(L[i+2]), (3,0))
    # elif forw:
    #     theta = np.tensordot(np.diag(L[i]),G[i], (1,1))
    #     theta = np.transpose(np.tensordot(theta, G[i+1], (2,1)), (1,0,2,3))
    
    phi = np.transpose(np.tensordot(Ut,theta,([1,3],[0,2])), (0,2,1,3))
    phi = np.reshape(phi, (d*chia,d*chib))
    
    U, S, V = np.linalg.svd(phi)
    V = V.T
    
    chic = min([np.sum(S>10**-16), chi])
    
    U = np.reshape(U[:chia*d,:chic], (d,chia,chic))
    V = np.transpose(np.reshape(V[:chib*d,:chic], (d,chib,chic)), (0,2,1))
    S = S[:chic]
    
    G[i] = np.transpose(np.tensordot(np.diag(L[i]**(-1)),U, (1,1)), (1,0,2))
    G[i+1] = np.tensordot(V, np.diag(L[i+2]**(-1)), (2,0))
    L[i+1] = S/np.sqrt(np.sum(S**2))
    
    # G[i] = U
    # G[i+1] = V
    # L[i+1] = S/np.sqrt(np.sum(S**2))
    
    return G, L
    
def get_timeop(H, dt, d, p = 1):
    Ut = scl.expm(-dt*p*H)
    
    Ut = np.reshape(np.transpose(np.reshape(Ut, (d,d,d,d)), (0,2,1,3)), (4,4))
    U, S, V = np.linalg.svd(Ut)
    U1 = np.reshape(np.tensordot(U,np.diag(np.sqrt(S)), (1,0)), (d,d,len(S)))
    U2 = np.reshape(np.tensordot(np.diag(np.sqrt(S)), V, (1,0)), (len(S),d,d))
    Ut = np.tensordot(U1,U2, (2,0))
    return Ut
    
def get_ener(G,L,Hlist):
    E = []
    for i in range(N-1):
        theta = np.tensordot(np.diag(L[i]), G[i], (1,1))
        theta = np.tensordot(theta, np.diag(L[i+1]), (2,0))
        theta= np.tensordot(theta, G[i+1], (2,1))
        theta = np.tensordot(theta, np.diag(L[i+2]), (3,0))
        
        # theta = np.tensordot(G[i],G[i+1],(2,1))
        # theta = np.transpose(np.tensordot(theta,np.diag(L[i+2]), (3,0)), (1,0,2,3))
        
        # theta = np.tensordot(np.diag(L[i]),G[i], (1,1))
        # theta = np.tensordot(theta, G[i+1], (2,1))
        
        C = np.tensordot(theta,np.reshape(Hlist[i], (2,2,2,2)), ([1,2],[2,3]))
        E.append(np.tensordot(np.conj(theta),C, ([0,3,1,2],[0,1,2,3])))
    return E
    
J = 1.
h = 0.1
N = 4
d = 2
chi = 8
dt = 0.05
t = 30
step_number = int(t/dt)

G = []
L = []
for i in range(N+1):
    if i < N:
        G.append(np.zeros([d,1,1]))
        
        G[-1][0,0,0] = 1/np.sqrt(2)
        G[-1][1,0,0] = 1/np.sqrt(2)
    
    L.append(np.ones([1]))
    
sz = np.array([[1.,0.],[0.,-1.]])
sy = np.array([[0,-complex(0,1)],[complex(0,1),0]])
sx = np.array([[0.,1.],[1.,0.]])

S2 = np.array([np.kron(sx,sx),np.kron(sy,sy),np.kron(sz,sz)]).real

H_int = J*(S2[0]+S2[1]+S2[2])
Hf = H_int + h*np.kron(sz,np.eye(2)) + h/2*np.kron(np.eye(2), sz)
Hb = H_int + h/2*np.kron(sz,np.eye(2))+h/2*np.kron(np.eye(2), sz)
Hl = H_int + h/2*np.kron(sz,np.eye(2)) + h*np.kron(np.eye(2), sz)

p = 1/2
Uf = get_timeop(Hf,dt,d,p)
Ub = get_timeop(Hb,dt,d)
Uo = get_timeop(Hb,dt,d,p)
Ul = get_timeop(Hl,dt,d,p)

I = np.eye(2)
I4 = np.transpose(np.reshape(np.kron(I,I), (2,2,2,2)), (0,2,1,3))
Ulist = [Uf]
Uodd = [Uf]
Ueven = [I4]
Hlist = [Hf]
for i in range(1,N-2):
    if np.mod(i,2) == 0:
        Ulist.append(Uo)
        Uodd.append(Uo)
        Ueven.append(I4)
    else:
        Ulist.append(Ub)
        Uodd.append(I4)
        Ueven.append(Ub)
    Hlist.append(Hb)
Ulist.append(Ul)
Uodd.append(Ul)
Ueven.append(I4)
Hlist.append(Hl)

forw = True

for l in range(step_number):
    if forw:
        for i in range(0,N-1,2):
            G,L = tebd(G,L,Ulist[i],i,forw)
        for j in range(N-3,-1,-2):
            G,L = tebd(G,L,Ulist[j],j, not forw)
        for k in range(0,N-1,2):
            G,L =tebd(G,L,Ulist[k],k, forw)
    else:
        for i in range(N-2,-1,-2):
            G,L = tebd(G,L,Ulist[i],i, forw)
        for j in range(1,N-1,2):
            G,L = tebd(G,L,Ulist[j],j, not forw)
        for k in range(N-2,-1,-2):
            G,L =tebd(G,L,Ulist[k],k, forw)
    forw = not forw

# end = N-1
# st = 0
# for l in range(step_number):
#     if forw:
#         for i in range(0,N-1):
#             G,L = tebd(G,L,Uodd[i],i,forw)
#         G[end] = np.transpose(np.tensordot(np.diag(L[end]),G[end], (1,1)), (1,0,2))
#         for j in range(N-2,-1,-1):
#             G,L = tebd(G,L,Ueven[j],j, not forw)
#         G[st] = np.tensordot(G[st],np.diag(L[st+1]), (2,0))
#         for k in range(0,N-1):
#             G,L =tebd(G,L,Uodd[k],k, forw)
#         G[end] = np.transpose(np.tensordot(np.diag(L[end]),G[end], (1,1)), (1,0,2))
#     else:
#         for i in range(N-2,-1,-1):
#             G,L = tebd(G,L,Uodd[i],i, forw)
#         G[st] = np.tensordot(G[st],np.diag(L[st+1]), (2,0))
#         for j in range(0,N-1):
#             G,L = tebd(G,L,Ueven[j],j, not forw)
#         G[end] = np.transpose(np.tensordot(np.diag(L[end]),G[end], (1,1)), (1,0,2))
#         for k in range(N-2,-1,-1):
#             G,L =tebd(G,L,Uodd[k],k, forw)
#         G[st] = np.tensordot(G[st],np.diag(L[st+1]), (2,0))
#     forw = not forw

# for i in range(0,N-1):
#     G, L = tebd(G,L,I4,i,forw)

# if forw:
#     for i in range(0,N-1):
#         G, L = tebd(G,L,I4,i,forw)
#     G[end] = np.transpose(np.tensordot(np.diag(L[end]),G[end], (1,1)), (1,0,2))
# elif not forw:
#     for i in range(N-2,-1,-1):
#         G, L = tebd(G,L,I4,i,forw)
#     G[st] = np.tensordot(G[st],np.diag(L[st+1]), (2,0))

psi = G[0]
for i in range(len(G)-1):
    psi = np.tensordot(psi,np.diag(L[i+1]), (i+2,0))
    psi = np.tensordot(psi,G[i+1], (i+2,1))

# psi = G[0]
# for i in range(len(G)-1):
#     psi = np.tensordot(psi,G[i+1], (i+2,1))

psi = np.reshape(psi, (d**N))
E = get_ener(G,L,Hlist)
print(E)
print(sum(E))

ED = ed.ExactD(J,h,N,dt,d)

print(ED.get_ener(psi))

