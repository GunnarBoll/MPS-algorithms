"""Storage for classes used in the tMPS algorithm
The class StateChain stores the current state in whatever notation is used 
and contains measurement of energy. It also contains the update procedure
given the state notation.

The class Hamiltonian constructs the two-site Hamiltonians of a given model
and the time evolution opertors which they generate. Further contains the time
evolution algorithms.

The class Measure contains measurement routines for (mostly) correlators.

The FreeFerm class contains algorithms for finding a free fermion solution
given a set of parameters (features its own measurement methods).
"""

import numpy as np
import scipy as sp
import itertools

# Class storing the state of a 1-D chain
class StateChain:
    def __init__(self, N, d, algo, bis_err=10**-10):
        self.N = N
        self.d = d
        self.L = []
        self.B = []
        self.err = 0
        self.bis_err = bis_err
        for i in range(N+1):
            if i < N:
                self.B.append(np.ones([d, 1, 1]))
                for j in range(d):
                    self.B[-1][j, 0, 0] = self.B[-1][j, 0, 0] / np.sqrt(d)
            self.L.append(np.ones([1]))
        if algo == "TEBD":
            self.notation = "LG"
        elif algo == "tDMRG":
            self.notation = "B"
        return
        
    def update(self, U, S, V, i, forward):
        if self.notation == "LG":
            self.B[i+1] = np.tensordot(V, np.diag(self.L[i+2] ** (-1)), (2, 0))
            self.B[i] = np.tensordot(np.diag(self.L[i] ** (-1)), U, (1, 1))
            self.B[i] = np.transpose(self.B[i], (1, 0, 2))
            self.L[i+1] = S / np.sqrt(np.sum(S ** 2))
        else:
            self.B[i] = U
            self.B[i + 1] = V
            self.L[i + 1] = S / np.sqrt(np.sum(S ** 2))
            if i == 0 and not forward:
                phi = np.tensordot(self.B[i], np.diag(self.L[i + 1]), (2, 0))
                phi = np.transpose(phi, (1, 0, 2))
                phi = np.reshape(phi, (1, self.d*len(self.L[i + 1])))
                U1, S1, V1 = np.linalg.svd(phi, full_matrices=False)
                self.B[i] = np.reshape(V1, (1, self.d, len(self.L[i + 1])))
                self.B[i] = np.transpose(self.B[i], (1, 0, 2))
                self.notation = "B"
            elif i == self.N-2 and forward:
                phi = np.tensordot(np.diag(self.L[i + 1]),
                                   self.B[i + 1], (1, 1))
                phi = np.reshape(np.transpose(phi, (1, 0, 2)),
                                 (self.d * len(self.L[i + 1]), 1))
                U1, S1, V1 = np.linalg.svd(phi, full_matrices=False)
                self.B[i + 1] = np.reshape(U1, (self.d, len(self.L[i + 1]), 1))
                self.notation = "A"
        return 
    
    def get_ener(self, Hchain):        
        E = []
        for b in range(0, self.N - 1):
            LBLBL = np.transpose(self.get_theta(b), (1, 0, 2, 3))
            C = np.tensordot(LBLBL, np.reshape(Hchain[b], (2, 2, 2, 2)),
                             ([1, 2], [2, 3]))
            C = np.tensordot(np.conj(LBLBL), C, ([0, 3, 1, 2], [0, 1, 2, 3]))
            E.append(C)
        return E
        
    def get_theta(self, b):
        if self.notation == "LG":
            LB = np.tensordot(np.diag(self.L[b]), self.B[b], (1, 1))
            LBL = np.tensordot(LB, np.diag(self.L[b + 1]), (2, 0))
            LBLB = np.tensordot(LBL, self.B[b + 1], (2, 1))
            LBLBL = np.tensordot(LBLB, np.diag(self.L[b + 2]), (3, 0))    
            theta = np.transpose(LBLBL, (1, 0, 2, 3))
        elif self.notation == "A":
            AA = np.tensordot(self.B[b], self.B[b + 1], (2, 1))
            AAL = np.tensordot(AA, np.diag(self.L[b + 2]), (3, 0))
            theta = AAL
        elif self.notation == "B":
            LB = np.tensordot(np.diag(self.L[b]), self.B[b], (1, 1))
            LBB = np.tensordot(LB, self.B[b + 1], (2, 1))
            theta = np.transpose(LBB, (1, 0, 2, 3))
        return theta



# Class containing MPS Hamiltonian and time evolution
class Hamiltonian:
    def __init__(self, g1, g2, N, dt, d, chi, model, TO, ED=False, grow_chi = True):
        self.g1 = g1
        self.g2 = g2
        self.d = d
        if grow_chi:    
            self.chi = 8
        else:
            self.chi = chi
        self.N = N
        self.dt = dt
        self.TO = TO
        self.model = model
        self.chi_max = chi
        self.Uall = []
        if model == "Heisen":
            self.Hlist = self.get_heisen(g1, g2)
        elif model == "HCboson":
            self.Hlist = self.get_HC_boson(g1, g2)
        else:
            return
        for n in range(len(self.Hlist)):
            if np.all(self.Hlist[n].imag == np.zeros(self.Hlist[n].shape)):
                self.Hlist[n] = self.Hlist[n].real
        self.Hchain = self.ener_chain(self.Hlist, N)
        I4 = np.reshape(np.kron(np.eye(d), np.eye(d)), (2, 2, 2, 2))
        self.I = np.transpose(I4, (0, 2, 1, 3))
        if ED:
            self.I = np.kron(np.eye(self.d), np.eye(self.d))
        # If first/fourth order do not put a half factor on odd operators
        if self.TO == "first":
            plist = [1]
            evenodd = True 
        elif self.TO == "second":
            plist = [1]
            evenodd = False
        elif self.TO == "fourth":
            plist = [1/12, -1/6]
            evenodd = True
        # Construct even and odd bond time evolution operators
        for p in plist:
            Uodd, Ueven = self.model_constructor(self.Hlist, p, self.I,
                                                 evenodd, ED)
            self.Uall.append([Uodd, Ueven])
    
    # Creates two-site Hamiltonians for a Heisenberg model in 1-D
    def get_heisen(self, g1, g2):
        Hlist = []
        parlist = [[g1, g2, g2/2], [g1, g2/2, g2/2], [g1, g2/2, g2]]
        for params in parlist:
            Hlist.append(self.kron_heisen(params[0], params[1], params[2]))
        return Hlist
        
    def get_HC_boson(self, g1, g2):
        Hlist = []
        g2 = np.array(g2)
        parlist = [[g1, g2, g2/2], [g1, g2/2, g2/2], [g1, g2/2, g2]]
        for params in parlist:
            Hlist.append(self.kron_hcb(params[0], params[1], params[2]))
        return Hlist
    
    def kron_hcb(self, t, mu1, mu2):
        adag = np.array([[0, 1], [0, 0]])
        a = np.array([[0, 0], [1, 0]])
        num_op = np.matmul(adag, a)
        #num_op = num_op - np.eye(self.d / 2)
        
        H = (- t[0] * (np.kron(adag, a)+np.kron(a, adag))
             - mu1[0] * np.kron(num_op, (np.eye(self.d)))
             - mu2[0] * np.kron(np.eye(self.d), num_op)
             + t[1] * np.kron(num_op, num_op)
             + mu1[1] * np.kron(adag+a, np.eye(self.d))
             + mu2[1] * np.kron(np.eye(self.d), adag+a)
             )
        return H
    
    # Construction of a single two-site Hamiltonian (Heisenberg model)
    def kron_heisen(self, J, h1, h2):
        sz = np.array([[1., 0.],[0., -1.]]) / 2
        sy = np.array([[0, -complex(0., 1.)],[complex(0., 1.), 0.]]) / 2
        sx = np.array([[0., 1.], [1., 0.]]) / 2
        S = np.array([sx, sy, sz])
        H_int = 0
        for i in range(3):
            H_int += J[i] * np.kron(S[i], S[i])
        H = H_int + h1 * np.kron(sz, np.eye(2)) + h2 * np.kron(np.eye(2), sz)
        return H
        
    # Constructs the bond Hamiltonians for energy calculation
    def ener_chain(self, Hlist, N):
        Hchain = []
        Hchain.append(Hlist[0])
        for i in range(N - 3):
            Hchain.append(Hlist[1])
        Hchain.append(Hlist[-1])
        return Hchain
    
    # Creates the time evolution for even and odd sites
    def model_constructor(self, Hlist, p, I, evenodd=False, ED=False):
        podd = peven = p
        if not evenodd:
            podd = podd / 2
        time_ops = []
        for ps, odd in [[podd, True], [peven, False]]:
            time_ops.append(self.timeop_chain(Hlist, I, ps, odd, ED))
        return time_ops
        
    # Constructs a chain of time evolution operators (even or odd spaces)
    def timeop_chain(self, Hlist, I, p, odd=True, ED=False):
        Ulist = []
        expHlist = []
        for H in Hlist:
            expHlist.append(self.get_timeop(H, p, ED))
        iter_list = [I, expHlist[1]]
        i = 0
        if not odd:
            expHlist[0] = expHlist[-1] = I
            iter_list[0], iter_list[1] = iter_list[1], iter_list[0]
        Ulist.append(expHlist[0])
        while len(Ulist) < self.N-2:
            Ulist.append(iter_list[np.mod(i, 2)])
            i += 1
        Ulist.append(expHlist[-1])
        return Ulist
    
    # Constructs a time evolution operator given a Hamiltonian
    def get_timeop(self, H, p = 1, ED=False):
        e, v = np.linalg.eig(H)
        DO = np.tensordot(np.diag(np.exp(-self.dt * p * e)),
                          np.linalg.inv(v), (1, 0))
        Ut = np.tensordot(v, DO, (1, 0))
        
        Ut = np.reshape(Ut, (self.d, self.d, self.d, self.d))
        Ut = np.reshape(np.transpose(Ut, (0, 2, 1, 3)), (4, 4))
        U, S, V = np.linalg.svd(Ut)
        U1 = np.reshape(np.tensordot(U, np.diag(np.sqrt(S)), (1, 0)),
                        (self.d, self.d, len(S)))
        U2 = np.reshape(np.tensordot(np.diag(np.sqrt(S)), V, (1, 0)),
                        (len(S), self.d, self.d))
        Ut = np.tensordot(U1, U2, (2, 0))
        
        assert Ut.shape == (self.d, self.d, self.d, self.d)
        
        if ED:
            Ut = np.reshape(np.transpose(Ut, (0, 2, 1, 3)),
                            (self.d**2, self.d**2))
        return Ut
    
    
    
    def time_evolve(self,Psi, step_number, algo):
        # First order algorithm
        if self.TO == "first":
            
            sweep_order = [self.Uall[0][0], self.Uall[0][1]]
            direc = True #True for forward sweep
            
            for t in range(step_number):
                for oper in sweep_order:
                    Psi = self.sweep(Psi, oper, algo, forward=direc)
                    direc = not direc
                
        
        # Second order algorithm
        elif self.TO == "second":
            # Time evolution operators for full time steps
            Utful = []
            for n in range(0, self.N-1, 2):
                fullstep_op = np.tensordot(self.Uall[0][0][n],
                                           self.Uall[0][0][n],
                                           ([1, 3], [0, 2]))
                Utful.append(np.transpose(fullstep_op, (0, 2, 1, 3)))
                if n != self.N-2:
                    Utful.append(self.Uall[0][0][n + 1])
            
            # Second order trotter sweeping order
            order = [self.Uall[0][1], Utful]
            Psi = self.sweep(Psi, self.Uall[0][0], algo, forward=True)
            
            Psi = self.sweeping_order(Psi, step_number-1, algo, order,
                                      forward=False)
            
            Psi = self.sweep(Psi, self.Uall[0][1], algo, forward=False)
            Psi = self.sweep(Psi, self.Uall[0][0], algo, forward=True)
            
            # Curative sweep
            Psi= self.sweep(Psi, list(itertools.repeat(self.I, self.N-1)),
                            algo, forward=False)
                    
        elif self.TO == "fourth":
            fourU = [[], []]
            for i in range(len(self.Uall)):
                for k in range(self.N - 1):
                    if np.mod(k, 2) == 0:
                        fourU[i].append(self.Uall[i][0][k])
                    else:
                        fourU[i].append(self.Uall[i][1][k])
            
            Uid = list(itertools.repeat(self.I, self.N-1))
            
            order = [fourU[0], fourU[0], fourU[0], fourU[1], fourU[0], Uid,
                     fourU[0], Uid,fourU[0], Uid, fourU[0], fourU[0], fourU[0],
                     fourU[0], Uid, fourU[0], Uid, fourU[0], Uid, fourU[0],
                     fourU[1], fourU[0], fourU[0], fourU[0]
                     ]
            
            Psi = self.sweeping_order(Psi, step_number, algo, order,
                                      forward=True)
            Psi= self.sweep(Psi, list(itertools.repeat(self.I, self.N-1)),
                            algo, forward=True)
            Psi= self.sweep(Psi, list(itertools.repeat(self.I, self.N-1)),
                            algo, forward=False)
        return Psi

    def sweeping_order(self, Psi, step_number, algo, order, forward=True):
        t = 0
        operlist = order
        while t < step_number:
            for oper in operlist:
                Psi = self.sweep(Psi, oper, algo, forward)
                forward = not forward
            t += 1
        if Psi.err > 10**-3:
            print("Warning: Accumulated truncation error is:", Psi.err)
        return Psi

    # Applies time evolution to the chain (sweeping forward or backward)
    def sweep(self, Psi, time_ops, algo, forward=True):
        sites = range(self.N - 1)
        if not forward:
            sites = range(self.N-2, -1, -1)
        
        for i in sites:
            # Change tebd function for other algorithm
            if algo == "TEBD":
                Psi = self.tebd(Psi, time_ops, i, forward)
            elif algo == "tDMRG":
                Psi = self.tdmrg(Psi, time_ops, i, forward)
            
        return Psi
    
    # TEBD algorithm (utilizing inverse lambda matrices)
    def tebd(self, Psi, Ut, i, forward):
        chia = Psi.B[i].shape[1]
        chib = Psi.B[i + 1].shape[2]
        
        theta = Psi.get_theta(i)
        # if np.all(Ut[i] == self.I):
        #     return Psi
        # else:
        phi = np.transpose(np.tensordot(Ut[i], theta, ([1, 3], [0, 2])),
                           (0, 2, 1, 3))
        phi = np.reshape(phi, (self.d * chia, self.d * chib))
        
        # Singular value decomposition    
        U, S, V, err, chic = self.svd_truncator(phi, chia, chib, Psi.bis_err)
        # Truncation error accumulation
        Psi.err += err
        # State update
        Psi.update(U, S, V, i, forward)
        return Psi
    
    # tDMRG algorithm for time evolution
    def tdmrg(self, Psi, Ut, i, forward):
        chia = Psi.B[i].shape[1]
        chib = Psi.B[i + 1].shape[2]
        
        theta = Psi.get_theta(i)
        # Time evolve tensor
        phi = np.transpose(np.tensordot(Ut[i], theta, ([1, 3], [0, 2])),
                           (0, 2, 1, 3))
        phi = np.reshape(phi, (self.d * chia, self.d * chib))
        
        U, S, V, err, chic = self.svd_truncator(phi, chia, chib, Psi.bis_err)
        # Truncation error accumulation
        Psi.err += err
        # State update
        Psi.update(U, S, V, i, forward)
        return Psi
    
    # Performs an SVD and truncates the singular values to specified bond
    # dimension.
    def svd_truncator(self, phi, chia, chib, max_err):
        U, S, V = sp.linalg.svd(phi, full_matrices=False)
        V = V.T
        
        chic = min([np.sum(S > 10**-14), self.chi])
        
        err = np.sum(S[chic:] ** 2)
        
        while err > max_err and self.chi < self.chi_max:
            self.chi += 1
            chic = min([np.sum(S > 10**-14), self.chi])
            err = np.sum(S[chic:] ** 2)
        
        S = S[: chic]
        U = np.reshape(U[:self.d*chia, :chic], (self.d, chia, chic))
        VT = np.reshape(V[:self.d*chib, :chic], (self.d, chib, chic))
        V = np.transpose(VT, (0, 2, 1))
        return [U, S ,V, err, chic]


class Measure:
    def __init__(self):
        return
    
    # Calculates correlation of two single site operators (op1, op2) at
    # site i and j resp.
    def correl(self, Psi, op1, op2, i, j):
        if j < i:
            i, j = j, i
        elif j == i:
            oper = np.matmul(op1, op2)
            return self.expec(Psi, oper, i)
        
        B_bar1 = np.tensordot(op1, Psi.B[i], (1, 0))
        B_bar2 = np.tensordot(op2, Psi.B[j], (1, 0))
        Phi = [B_bar1]
        psi_bra = [np.conj(Psi.B[i])]
        for l in range(i+1, j):
            Phi.append(Psi.B[l])
            psi_bra.append(np.conj(Psi.B[l]))
        psi_bra.append(np.conj(Psi.B[j]))
        Phi.append(B_bar2)
        
        if Psi.notation == "B":
            Phi[0] = np.tensordot(np.diag(Psi.L[i]), Phi[0], (1, 1))
            Phi[0] = np.transpose(Phi[0], (1, 0, 2))
            psi_bra[0] = np.tensordot(np.conj(np.diag(Psi.L[i])),
                                      psi_bra[0], (1, 1))
            psi_bra[0] = np.transpose(psi_bra[0], (1, 0, 2))
        elif Psi.notation == "A":
            Phi[-1] = np.tensordot(Phi[-1], np.diag(Psi.L[j + 1]), (2, 0))
            psi_bra[-1] = np.tensordot(psi_bra[-1],
                                       np.conj(np.diag(Psi.L[j + 1])), (2, 0))
        
        corr = np.tensordot(psi_bra[0], Phi[0], ([0, 1], [0, 1]))
        for k in range(1, j-i+1):
            corr = np.tensordot(corr, psi_bra[k], (0, 1))
            corr = np.tensordot(corr, Phi[k], ([1, 0], [0, 1]))
        corr = np.trace(corr)
        return corr
    
    def expec(self, Psi, op, i):
        B_bar = np.tensordot(np.diag(Psi.L[i]),
                             np.tensordot(op, Psi.B[i], (1, 0)), (1, 1))
        B_psi = np.tensordot(np.conj(np.diag(Psi.L[i])), np.conj(Psi.B[i]),
                           (1, 1))
        exval = np.tensordot(B_psi, B_bar, ([1, 0], [1, 0]))
        exval = np.trace(exval)
        return exval
    
    def corr_mat(self, Psi, op1, op2):
        mat = np.zeros([Psi.N, Psi.N])
        expec_op = np.matmul(op1, op2)
        oplist = [op1, op2]
        revoplist = [op2, op1]
            
        
        for ind1 in range(Psi.N):
            left = []
            for op in oplist:
                B_bar = np.tensordot(op, Psi.B[ind1], (1, 0))
                if Psi.notation == "B":
                    B_bar = np.tensordot(np.diag(Psi.L[ind1]), B_bar, (1, 1))
                    B_star = np.tensordot(np.diag(Psi.L[ind1]), Psi.B[ind1],
                                          (1, 1))
                    B_bar = np.transpose(B_bar, (1, 0, 2))
                    B_star = np.transpose(B_star, (1, 0, 2))
                else:
                    B_star = Psi.B[ind1]
                B_pair = np.tensordot(B_star.conj(), B_bar, ([0, 1], [0, 1]))
                left.append(B_pair)
            mat[ind1, ind1] = self.expec(Psi, expec_op, ind1)
            for ind2 in range(ind1+1, Psi.N):
                for ind3 in range(2):
                    Rp = np.tensordot(revoplist[ind3], Psi.B[ind2], (1, 0))
                    if Psi.notation == "A":
                        B_right = np.tensordot(Psi.B[ind2],
                                               np.diag(Psi.L[ind2 + 1]),
                                               (2, 0))
                        Rp = np.tensordot(Rp, Psi.L[ind2 + 1], (2, 0))
                    else:
                        B_right = Psi.B[ind2]
                    Lp = np.tensordot(B_right.conj(), left[ind3], (1, 0))
                    mat[ind1, ind2] = np.trace(np.tensordot(Lp, Rp,
                                               ([0, 2], [0, 1])))
                    left[ind3] = np.tensordot(left[ind3], Psi.B[ind2], (1, 1))
                    left[ind3] = np.tensordot(Psi.B[ind2].conj(), left[ind3],
                                              ([0, 1], [1, 0]))
                    mat = mat.T
        return mat
    
class FreeFerm:
    def __init__(self, t, mu, N):
        self.t = t[0]
        self.mu = mu[0]
        self.N = N
        
        self.E_GS = 0
        self.H = np.zeros([N, N])
        
        for i in range(N - 1):
            self.H[i, i] = -self.mu
            self.H[i, i+1] = -self.t
            self.H[i+1, i] = -self.t
        self.H[N-1, N-1] = -self.mu
        
        e, v = np.linalg.eigh(self.H)
        e.sort()
        
        N_f = sum(e < 0)
        
        self.P = v[:, :N_f]

        for ind in range(N_f):
            self.E_GS += e[ind]
        # self.E_GS += N * self.mu/2
        return
    
    def measure(self, opl, indl):
        chain = OpChain(opl, indl, sign=1)
        chains = self.order(chain)
        res = 0
        for opchain in chains:
            res += opchain.sign * self.correl(opchain.ind[:opchain.acount],
                                              opchain.ind[opchain.acount:])
        return res
    
    def order(self, opers):
        rest = []        
        for i in range(len(opers.op)-1, -1, -1):
            if opers.op[i] == 'adag':
                for j in range(i+1, len(opers.op)):
                    if opers.op[j] == 'adag':
                        break
                    opers.swap(j-1, j)
                    if opers.ind[j - 1] == opers.ind[j]:
                        rest.append(OpChain(opers.op[:j-1] + opers.op[j+1:],
                                             opers.ind[:j-1] + opers.ind[j+1:],
                                             opers.sign))
                        opers.sign = -opers.sign
        new_stuff = [opers]
        for k in range(len(rest)):
            new_stuff = new_stuff + self.order(rest[k])
        return new_stuff
        
    def correl(self,a_list,dag_list):
        P_A = self.P
        P_B = self.P
        dag_list.reverse()
        for a_ind in a_list:
            P_A = self.get_Pmat(P_A, a_ind)
        for dag_ind in dag_list:
            P_B = self.get_Pmat(P_B, dag_ind)
        
        corr = np.linalg.det(np.matmul(np.conj(P_A).T, P_B))
        return corr
    
    def get_Pmat(self, P, row):
        part = -P[:row]
        Zcol = np.zeros([self.N, 1])
        Zcol[row] = 1
        
        P = np.hstack([np.concatenate([part, P[row:]]), Zcol])
        return P

class OpChain:
    def __init__(self, oplist, indlist, sign):
        self.op = oplist
        self.ind = indlist
        self.sign = sign
        self.acount = 0
        for op in oplist:
            if op == "a":
                self.acount += 1
    
    def swap(self, i, j):
        self.op[i], self.op[j] = self.op[j], self.op[i]
        self.ind[i], self.ind[j] = self.ind[j], self.ind[i]
        