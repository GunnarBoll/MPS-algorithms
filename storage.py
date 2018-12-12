"""Storage for classes used in the tMPS algorithm
The class state_chain stores the current state in whatever notation is used and contains
measurement of energy (to be moved). It also contains the update procedure given the state notation

The class Hamiltonian constructs the two-site Hamiltonians of a given model and the time evolution
opertors which they generate. Further contains the time evolution algorithms.
"""

import numpy as np
import itertools

#Class for the state
class StateChain:
    def __init__(self,N, chi, d, algo):
        self.N = N
        self.d = d
        self.chi = chi
        self.L = []
        self.B = []
        self.err = 0
        for i in range(N+1):
            if i < N:
                self.B.append(np.ones([d,1,1]))
                for j in range(d):
                    self.B[-1][j,0,0] = self.B[-1][j,0,0]/np.sqrt(d)
            self.L.append(np.ones([1]))
        
        
        if algo == "TEBD":
            self.notation = "LG"
        elif algo == "tDMRG":
            self.notation = "B"
        
        return
        
        
    
    def update(self, U, S, V, i, forward):
        if self.notation == "LG":
            self.B[i+1] = np.tensordot(V, np.diag(self.L[i+2]**(-1)), (2,0))
            self.B[i] = np.transpose(np.tensordot(np.diag(self.L[i]**(-1)),U, (1,1)), (1,0,2))
            self.L[i+1] = S/np.sqrt(np.sum(S**2))
        else:
            self.B[i] = U
            self.B[i+1] = V
            self.L[i+1] = S/np.sqrt(np.sum(S**2))
            
            if i==0 and not forward:
                phi = np.reshape(np.transpose(np.tensordot(self.B[i], np.diag(self.L[i+1]), (2,0)), (1,0,2)), (1,self.d*len(self.L[i+1])))
                U1, S1, V1 = np.linalg.svd(phi,full_matrices = False)
                self.B[i] = np.transpose(np.reshape(V1, (1,self.d,len(self.L[i+1]))), (1,0,2))
                self.notation = "B"

            elif i == self.N-2 and forward:
                phi = np.reshape(np.transpose(np.tensordot(np.diag(self.L[i+1]), self.B[i+1], (1,1)), (1,0,2)), (self.d*len(self.L[i+1]), 1))
                U1, S1, V1 = np.linalg.svd(phi, full_matrices = False)
                self.B[i+1] = np.reshape(U1, (self.d, len(self.L[i+1]), 1))
                self.notation = "A"
        return 
    
    def get_ener(self, Hchain):        
        E = []
        for b in range(0,self.N-1):
            LBLBL = np.transpose(self.get_theta(b), (1,0,2,3))
            C = np.tensordot(LBLBL, np.reshape(Hchain[b], (2,2,2,2)), ([1,2],[2,3]))
            E.append(np.tensordot(np.conj(LBLBL), C, ([0,3,1,2],[0,1,2,3])))
        return E
        
    def get_theta(self, b):
        if self.notation == "LG":
            LB = np.tensordot(np.diag(self.L[b]),self.B[b],(1,1))
            LBL = np.tensordot(LB,np.diag(self.L[b+1]), (2,0))
            LBLB = np.tensordot(LBL,self.B[b+1], (2,1))
            LBLBL = np.tensordot(LBLB, np.diag(self.L[b+2]), (3,0))    
            theta = np.transpose(LBLBL, (1,0,2,3))
        elif self.notation == "A":
            AA = np.tensordot(self.B[b],self.B[b+1], (2,1))
            AAL = np.tensordot(AA, np.diag(self.L[b+2]), (3,0))
            theta = AAL
        elif self.notation == "B":
            LB = np.tensordot(np.diag(self.L[b]), self.B[b], (1,1))
            LBB = np.tensordot(LB, self.B[b+1], (2,1))
            theta = np.transpose(LBB, (1,0,2,3))
        return theta
        

#Class for Hamiltonian
class Hamiltonian:
    def __init__(self,J,h,N,dt,d, chi, model, TO, ED = False):
        self.J = J
        self.h = h
        self.d = d
        self.chi = 8
        self.N = N
        self.dt = dt
        self.TO = TO
        self.model = model
        self.chi_max = chi
        self.Uall = []
        
        if model == "Heisen":
            self.Hlist = self.get_heisen(J,h)
        elif model == "HCboson":
            self.Hlist = self.get_HC_boson(J,h)
        else:
            return
        
        for n in range(len(self.Hlist)):
            if np.all(self.Hlist[n].imag == np.zeros(self.Hlist[n].shape)):
                self.Hlist[n] = self.Hlist[n].real            
        
        self.Hchain = self.ener_chain(self.Hlist, N)
        self.I = np.transpose(np.reshape(np.kron(np.eye(d),np.eye(d)), (2,2,2,2)), (0,2,1,3))
        if ED:
            self.I = np.kron(np.eye(self.d),np.eye(self.d))
        
        #If first/fourth order do not put a half factor on odd operators
        if self.TO == "first":
            plist = [1]
            evenodd = True 
        
        elif self.TO == "second":
            plist = [1]
            evenodd = False
            
        elif self.TO == "fourth":
            plist = [1/12,-1/6]
            evenodd = True
        
        #Construct even and odd bond time evolution operators
        for p in plist:
            Uodd, Ueven = self.model_constructor(self.Hlist, p, self.I, evenodd, ED)
            self.Uall.append([Uodd,Ueven])
        return
    
    #Creates two-site Hamiltonians for an Anti-Ferromagnetic Transverse Ising Chain (AFTIC)
    def get_heisen(self,J,h):
        Hlist = []
        parlist = [[J,h,h/2],[J,h/2,h/2],[J,h/2,h]]
        for params in parlist:
            Hlist.append(self.kron_heisen(params[0],params[1],params[2]))
        return Hlist
        
    def get_HC_boson(self, t, mu):
        Hlist = []
        parlist = [[t,mu,mu/2],[t,mu/2,mu/2],[t,mu/2,mu]]
        for params in parlist:
            Hlist.append(self.kron_hcb(params[0],params[1],params[2]))
        return Hlist
    
    def kron_hcb(self,t,mu1,mu2):
        adag = np.array([[0,1],[0,0]])
        a = np.array([[0,0],[1,0]])
        n = np.matmul(adag,a)

        H = (-t[0]*(np.kron(adag,a) + np.kron(a,adag)) - mu1*np.kron(n-np.eye(self.d)/2,np.eye(self.d)) - 
        mu2*np.kron(np.eye(self.d),n-np.eye(self.d)/2) + t[1]*np.kron(n-np.eye(self.d)/2, n-np.eye(self.d)/2) )
        return H
    
    #Construction of a single two-site Hamiltonian (AFTIC)
    def kron_heisen(self,J,h1,h2):
        sz = np.array([[1.,0.],[0.,-1.]])/2
        sy = np.array([[0,-complex(0,1)],[complex(0,1),0]])/2
        sx = np.array([[0.,1.],[1.,0.]])/2
        S = np.array([sx, sy, sz])
        H_int = 0
        for i in range(3):
            H_int += J[i]*np.kron(S[i],S[i])
        H = H_int + h1*np.kron(sz,np.eye(2)) + h2*np.kron(np.eye(2),sz)
        return H
        
    #Constructs the bond Hamiltonians for energy calculation
    def ener_chain(self,Hlist,N):
        Hchain = []
        Hchain.append(Hlist[0])
        for i in range(N-3):
            Hchain.append(Hlist[1])
        Hchain.append(Hlist[-1])
        return Hchain
    
    #Creates the time evolution for even and odd sites
    def model_constructor(self, Hlist, p,I,evenodd = False, ED = False):
        podd = peven = p
        if not evenodd:
            podd = podd/2
        time_ops = []
        for ps, odd in [[podd,True],[peven, False]]:
            time_ops.append(self.timeop_chain(Hlist,I,ps,odd, ED))
        return time_ops
        
    #Constructs a chain of time evolution operators (even or odd spaces)
    def timeop_chain(self,Hlist, I, p, odd = True, ED = False):
        Ulist = []
        expHlist = []
        for H in Hlist:
            expHlist.append(self.get_timeop(H,p, ED))
        
        iter_list = [I, expHlist[1]]
        i = 0
        if not odd:
            expHlist[0] = expHlist[-1] = I
            iter_list[0], iter_list[1] = iter_list[1], iter_list[0]
        
        Ulist.append(expHlist[0])
        while len(Ulist)<self.N-2:
            Ulist.append(iter_list[np.mod(i,2)])
            i+=1
        Ulist.append(expHlist[-1])
        
        return Ulist
    
    #Constructs a time evolution operator given a Hamiltonian
    def get_timeop(self,H, p = 1, ED = False):
        e, v = np.linalg.eig(H)
        Ut = np.tensordot(v,np.tensordot(np.diag(np.exp(-self.dt*p*e)),np.linalg.inv(v),(1,0)),(1,0))
        
        Ut = np.reshape(np.transpose(np.reshape(Ut, (self.d,self.d,self.d,self.d)), (0,2,1,3)), (4,4))
        U, S, V = np.linalg.svd(Ut)
        U1 = np.reshape(np.tensordot(U,np.diag(np.sqrt(S)), (1,0)), (self.d,self.d,len(S)))
        U2 = np.reshape(np.tensordot(np.diag(np.sqrt(S)), V, (1,0)), (len(S),self.d,self.d))
        Ut = np.tensordot(U1,U2, (2,0))
                
        if ED:
            Ut = np.reshape(np.transpose(Ut, (0,2,1,3)), (self.d**2,self.d**2))
        
        return Ut
    
    
    
    def time_evolve(self,Psi, step_number, algo):
        #First order algorithm
        if self.TO == "first":
            
            sweep_order = [self.Uall[0][0], self.Uall[0][1]]
            direc = True #True for forward sweep
            
            for t in range(step_number):
                for oper in sweep_order:
                    Psi = self.sweep(Psi, oper, algo, forward = direc)
                    direc = not direc
                
        
        #Second order algorithm
        elif self.TO == "second":
            #Time evolution operators for full time steps
            Utful = []
            for n in range(0,self.N-1,2):
                Utful.append(np.transpose( np.tensordot(self.Uall[0][0][n], self.Uall[0][0][n], ([1,3],[0,2]) ), (0,2,1,3)))
                if n != self.N-2:
                    Utful.append(self.Uall[0][0][n+1])
            
            #Second order trotter sweeping order
            order = [self.Uall[0][1], Utful]
            Psi = self.sweep(Psi,self.Uall[0][0],algo, forward = True)
            
            Psi = self.sweeping_order(Psi, step_number-1, algo, order, forward = False)
            
            Psi = self.sweep(Psi,self.Uall[0][1],algo, forward = False)
            Psi = self.sweep(Psi,self.Uall[0][0],algo, forward = True)
            
            #Fixes the error??
            Psi= self.sweep(Psi, list(itertools.repeat(self.I,self.N-1)), algo, forward = False)
                    
        elif self.TO == "fourth":
            fourU = [[],[]]
            for i in range(len(self.Uall)):
                for k in range(self.N-1):
                    if np.mod(k,2) == 0:
                        fourU[i].append(self.Uall[i][0][k])
                    else:
                        fourU[i].append(self.Uall[i][1][k])
            
            Uid = list(itertools.repeat(self.I,self.N-1))
            
            order = [fourU[0], fourU[0], fourU[0], fourU[1], fourU[0], Uid, fourU[0], Uid,fourU[0], Uid, fourU[0], fourU[0], fourU[0],
            fourU[0], Uid, fourU[0], Uid, fourU[0], Uid, fourU[0], fourU[1], fourU[0], fourU[0], fourU[0]]
            
            Psi = self.sweeping_order(Psi, step_number, algo, order, forward = True)
            Psi= self.sweep(Psi, list(itertools.repeat(self.I,self.N-1)), algo, forward = True)
            Psi= self.sweep(Psi, list(itertools.repeat(self.I,self.N-1)), algo, forward = False)
            
        return Psi

    def sweeping_order(self, Psi, step_number, algo, order, forward = True):
        t = 0
        operlist = order
        while t < step_number:
            for oper in operlist:
                Psi = self.sweep(Psi, oper, algo, forward)
                forward = not forward
                if Psi.err > 10**-6:
                    break
            
            t += 1
        return Psi

    #Applies time evolution to the chain (sweeping forward or backward)
    def sweep(self,Psi, time_ops, algo, forward = True):
        sites = range(self.N-1)
        if not forward:
            sites = range(self.N-2, -1, -1)
        
        for i in sites:
            #Change tebd function for other algorithm
            if algo == "TEBD":
                Psi = self.tebd(Psi, time_ops, i, forward)
            elif algo == "tDMRG":
                Psi = self.tdmrg(Psi, time_ops, i, forward)
            
            if Psi.err > 10**-6:
                print("Truncation error too big! Evolution stopped.")
                break
        return Psi
    
    #Main TEBD algorithm (utilizing inverse lambda matrices)
    def tebd(self,Psi,Ut, i, forward):
        chia = Psi.B[i].shape[1]
        chib = Psi.B[i+1].shape[2]
        
        theta = Psi.get_theta(i)
        # if np.all(Ut[i] == self.I):
        #     return Psi
        # else:
        phi = np.transpose(np.tensordot(Ut[i], theta, ([1,3], [0,2])), (0,2,1,3))
        
        phi = np.reshape(phi, (self.d*chia,self.d*chib))
        
        #Singular value decomposition    
        U, S, V, err, chic = self.svd_truncator(phi, chia, chib)
        
        Psi.err += err
        #State update
        Psi.update(U,S,V,i,forward)
        
        
        return Psi
    
    #tDMRG algorithm for time evolution
    def tdmrg(self,Psi, Ut, i, forward):
        chia = Psi.B[i].shape[1]
        chib = Psi.B[i+1].shape[2]
        
        theta = Psi.get_theta(i)
        
        phi = np.transpose(np.tensordot(Ut[i], theta, ([1,3],[0,2])), (0,2,1,3))
        
        phi = np.reshape(phi,(self.d*chia,self.d*chib))
        
        U, S, V, err, chic = self.svd_truncator(phi, chia, chib)
        
        Psi.err += err
        Psi.update(U,S,V, i, forward)
    
        return Psi
    
    #Performs an SVD and truncates the singular values to specified bond dimension
    def svd_truncator(self,phi, chia, chib):
        
        U, S, V = np.linalg.svd(phi,full_matrices=False)
        V = V.T
        
        chic = min([np.sum(S>10**-14), self.chi])
        
        err = 2*np.sum(S[chic:]**2)
        
        while err > 10**-10:
            if self.chi < self.chi_max:
                self.chi += 1
            chic = min([np.sum(S>10**-16), self.chi])
            err = 2*np.sum(S[chic:]**2)
        
        err = np.sum(S[chic:]**2)
        S = S[:chic]
        U = np.reshape(U[:self.d*chia,:chic], (self.d,chia,chic))
        V = np.transpose(np.reshape(V[:self.d*chib, :chic], (self.d,chib,chic)), (0,2,1))
        
        return [U, S ,V, err, chic]


class Measure:
    def __init__(self):
        
        return
    
    #Calculates correlation of two single site operators (op1, op2) at site i and j resp.
    def correl(self, Psi, op1, op2, i, j):
        if j < i:
            temp = i
            i = j
            j = temp
        elif j == i:
            oper = np.matmul(op1,op2)
            return self.expec(Psi,oper,i)
        
        B_bar1 = np.tensordot(op1,Psi.B[i], (1,0))
        B_bar2 = np.tensordot(op2,Psi.B[j], (1,0))
        
        Phi = [B_bar1]
        psi_bra = [np.conj(Psi.B[i])]
        for l in range(i+1,j):
            Phi.append(Psi.B[l])
            psi_bra.append(np.conj(Psi.B[l]))
        psi_bra.append(np.conj(Psi.B[j]))
        Phi.append(B_bar2)
        
        if Psi.notation == "B":
            Phi[0] = np.transpose(np.tensordot(np.diag(Psi.L[i]),Phi[0], (1,1)), (1,0,2))
            psi_bra[0] = np.transpose(np.tensordot(np.conj(np.diag(Psi.L[i])),psi_bra[0], (1,1)), (1,0,2))
        elif Psi.notation == "A":
            Phi[-1] = np.tensordot(Phi[-1],np.diag(Psi.L[j+1]), (2,0))
            psi_bra[-1] = np.tensordot(psi_bra[-1],np.conj(np.diag(Psi.L[j+1])), (2,0))
        
        
        corr = np.tensordot(psi_bra[0],Phi[0], ([0,1],[0,1]))
        for k in range(1,j-i+1):
            corr = np.tensordot(corr, psi_bra[k], (0,1))
            corr = np.tensordot(corr, Phi[k],([1,0],[0,1]))
        
        corr = np.trace(corr)
        
        return corr
    
    def expec(self, Psi, op, i):
        B_bar = np.tensordot(np.diag(Psi.L[i]),np.tensordot(op,Psi.B[i], (1,0)), (1,1))
        psi = np.tensordot(np.conj(np.diag(Psi.L[i])),np.conj(Psi.B[i]),(1,1))
        exval = np.tensordot(psi,B_bar, ([1,0],[1,0]))
        exval = np.trace(exval)
        
        return exval
    
    def struc_func(self,Psi, i, op):
        corrs = []
        
        for l in range(Psi.N-1,i,-1):
            corrs.append(self.correl(Psi,op,op,i,l))
        
        
        return corrs
    
class Freeferm:
    def __init__(self,t,mu,N):
        self.t = t
        self.mu = mu
        self.N = N
        
        self.E_GS = 0
        self.H = np.zeros([N,N])
        
        for i in range(N-1):
            self.H[i,i] = -mu
            self.H[i,i+1] = -t
            self.H[i+1,i] = -t
        self.H[N-1,N-1] = -mu
        
        e, v = np.linalg.eigh(self.H)
        
        for ener in e:
            if not ener < 0:
                break
            self.E_GS += ener
        return