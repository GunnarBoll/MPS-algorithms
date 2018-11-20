
import numpy as np
from tMPS import tebd
from tMPS import get_ener

J = 1.0
h = 1.0
dt = 0.01
step_number = int(6.5/dt)
chi = 8
N = 4
d = 2

sz = np.array([[1.,0.],[0.,-1.]]) 
sx = np.array([[0.,1.],[1.,0.]])

Hfirst = J*np.kron(sz,sz) + h*np.kron(sx,np.eye(2)) + h/2*np.kron(np.eye(2),sx)
Hbond = J*np.kron(sz,sz) + h/2*np.kron(sx,np.eye(2)) + h/2*np.kron(np.eye(2),sx)
Hlast = J*np.kron(sz,sz) + h/2*np.kron(sx,np.eye(2)) + h*np.kron(np.eye(2),sx)

ef, vf = np.linalg.eig(Hfirst)
eb, vb = np.linalg.eig(Hbond)
el, vl = np.linalg.eig(Hlast)

expHf = np.matmul(vf, np.matmul(np.diag(np.exp(-dt*ef/2)), np.linalg.inv(vf)))
expHo = np.matmul(vb, np.matmul(np.diag(np.exp(-dt*eb/2)),np.linalg.inv(vb)))
expHe = np.matmul(vb, np.matmul(np.diag(np.exp(-dt*eb)),np.linalg.inv(vb)))
expHl = np.matmul(vl, np.matmul(np.diag(np.exp(-dt*el/2)), np.linalg.inv(vl)))

U = []
Hc = []

for n in range(N-1):
    if n == 0:
        U.append(expHf)
        Hc.append(Hfirst)
    elif n == N-2:
        U.append(expHl)
        Hc.append(Hlast)
    elif np.mod(n,2) == 0:
        U.append(expHo)
        Hc.append(Hbond)
    else:
        U.append(expHe)
        Hc.append(Hbond)

Utfull = []
for i in range(0,N-1,2):
    Utfull.append(np.matmul(U[i],U[i]))
    Utfull[i] = np.transpose(np.reshape(Utfull[i],(2,2,2,2)), (0,2,1,3))
    U[i] = np.transpose(np.reshape(U[i],(2,2,2,2)), (0,2,1,3))
    if i != N-2:
        Utfull.append(U[i+1])
        Utfull[i+1] = np.transpose(np.reshape(Utfull[i+1],(2,2,2,2)), (0,2,1,3))
        U[i+1] = np.transpose(np.reshape(U[i+1],(2,2,2,2)), (0,2,1,3))


B = []
Lam = []
for i in range(N+1):
    B.append(np.zeros([2,1,1]))
    
    B[-1][0,0,0] = 1
    
    Lam.append(np.ones([1]))



for k in range(0,N-1,2):
    B, Lam = tebd(B, Lam, U, k, chi, d)

for t in range(step_number-1):
    for i in range(N-3,0,-2):
        B, Lam = tebd(B, Lam, Utfull, i, chi, d)
    for j in range(0,N-1,2):
        B, Lam = tebd(B, Lam, Utfull, j, chi, d)
        
for l in range(N-3,0,-2):
    B, Lam = tebd(B, Lam, Utfull, l, chi, d)

for m in range(0,N-1,2):
    B, Lam = tebd(B, Lam, U, m, chi, d)

E = get_ener(Lam, B, Hc, N)
print(sum(E))