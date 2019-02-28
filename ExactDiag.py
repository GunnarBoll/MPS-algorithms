""" 
File contains a class with methods that feature the diagonalization of a
Hamiltonian describing a 1D chain. The class has methods for simulating a
Trotter step to compare with a DMRG algorithm if needed.
Written by: Gunnar Bollmark
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sp_linalg
import itertools
import imp

import storage

imp.reload(storage)

# The class ExactD has Hamiltonian from storage.py as parent
class ExactD(storage.Hamiltonian):
    def __init__(self, g1, g2, N, d, model, TO="second", dt=0.1):
        
        super(ExactD,self).__init__(g1, g2, N, dt, d, 8, model, TO, ED=True)
                
        initspinor = np.array([1, 1]) / np.sqrt(2)
        self.i_state = Nkron(list(itertools.repeat(initspinor, self.N)))
        
        self.Id = sp.eye(d) # Identity matrix in sparse matrix format
        self.op = list(itertools.repeat(self.Id, N))
        
        # Create full Hilbert space Hamiltonian based on what model is used
        self.H = self.full_ham(model)
        
        # List of the time operators that will be applied each bond when using
        # Trotter stepping.
        self.Ulist = []
        self.Ulist2 = []
        for i in range(N - 1):
            if np.mod(i, 2) == 0:
                self.Ulist.append(sp.coo_matrix(self.Uall[0][0][i]))
                try:
                    self.Ulist2.append(sp.coo_matrix(self.Uall[1][0][i]))
                except IndexError:
                    pass
            elif np.mod(i,2) == 1:
                self.Ulist.append(sp.coo_matrix(self.Uall[0][1][i]))
                try:
                    self.Ulist2.append(sp.coo_matrix(self.Uall[1][1][i]))
                except IndexError:
                    pass
        return
    
    # Function which creates a Hamiltonian matrix in coordinate sparse
    # formulation dependent on what model is passed.
    def full_ham(self, model):
        H = 0
        if model == "Heisen":
            ssop, tsop, glist, hlist = self.heisen_opers()
        elif model == "HCboson":
            ssop, tsop, glist, hlist = self.hcb_opers()
        for i in range(self.N - 1):
            for k in range(len(tsop)):
                H += glist[k] * spkron([sp.eye(self.d ** i), spkron(tsop[k]),
                                        sp.eye(self.d ** (self.N-i-2))
                                        ]
                                      ).real
            for j in range(len(ssop)):
                H += hlist[j] * spkron([sp.eye(self.d ** i), ssop[j],
                                      sp.eye(self.d ** (self.N-i-1))
                                      ]
                                    )
        for l in range(len(ssop)):
            H += hlist[l] * spkron([sp.eye(self.d ** (self.N-1)), ssop[l]])
        return H
    
    # Function which returns single particle and two-particle operators and
    # a list of coupling constants for the Heisenberg model.
    def heisen_opers(self):
        sz = sp.coo_matrix([[1., 0.], [0., -1.]]) / 2
        sy = sp.coo_matrix([[0, -complex(0, 1)], [complex(0, 1), 0]]) / 2
        sx = sp.coo_matrix([[0., 1.], [1., 0.]]) / 2
        
        tsop = [[sx, sx], [sy, sy], [sz, sz]]
        ssop = [sz]
        return ssop, tsop, self.g1, [self.g2]
    
    # Function which returns single particle and two-particle operators and
    # a list of coupling constants for the hardcore boson model.
    def hcb_opers(self):
        adag = sp.coo_matrix([[0, 1], [0, 0]])
        a = sp.coo_matrix([[0, 0], [1, 0]])
        n_op = adag*a # - sp.eye(self.d)/2
        alp1 = 4 * self.g2[1] / np.sqrt(self.N/2)
        alp2 = 4 * self.g2[1] / np.sqrt(self.N/2)
        
        tsop = [[adag, a], [a, adag], [n_op, n_op]]
        ssop = [n_op, adag, a]
        glist = [-self.g1[0], -self.g1[0], self.g1[1]]
        hlist = [-self.g2[0], self.g2[1], self.g2[1]]
        return ssop, tsop, glist, hlist
    
    # Function which time evolves an initial state in time steps by splitting 
    # up the time evolution operator. This can be done to different orders in 
    # the error. Possible orders: first, second and fourth
    def trotter_time(self, step_number):
        state = self.i_state
        Iden = self.op[:self.N-1]
        if self.TO == "first":
            for t in range(step_number):
                state = self.ED_sweep(self.Uall[0][0], Iden, state, self.Id,
                                      forw=True)
                state = self.ED_sweep(self.Uall[0][1], Iden, state, self.Id,
                                      forw=False)
        elif self.TO == "second":
            Utfull = []
            for i in range(0, self.N-1, 2):
                full_mat = np.matmul(self.Uall[0][0][i], self.Uall[0][0][i])
                Utfull.append(full_mat)
                if i != self.N-2:
                    Utfull.append(self.Uall[0][0][i + 1])
            
            state = self.ED_sweep(self.Uall[0][0], Iden, state, self.Id,
                                  forw=True)
            order = [self.Uall[0][1], Utfull]
            
            state = self.ED_sweeping_order(state, step_number-1, order, Iden,
                                           forward=False)
            
            state = self.ED_sweep(self.Uall[0][1], Iden, state, self.Id,
                                  forw=False)
            state = self.ED_sweep(self.Uall[0][0], Iden, state, self.Id,
                                  forw=True)
        
        elif self.TO == "fourth":
            Uid = list(itertools.repeat(np.kron(self.Id,self.Id),self.N-1))
            order = [self.Ulist, self.Ulist, self.Ulist, self.Ulist2,
                     self.Ulist, Uid, self.Ulist, Uid, self.Ulist, Uid,
                     self.Ulist, self.Ulist, self.Ulist, self.Ulist, Uid,
                     self.Ulist, Uid, self.Ulist, Uid, self.Ulist, self.Ulist2,
                     self.Ulist, self.Ulist, self.Ulist]
            
            state = self.ED_sweeping_order(state, step_number, order, Iden,
                                           forward=False)
            
        self.trot_state = state
        return
    
    # Function which sweeps through the state with time evolution given some
    # order of time evolution operators.
    def ED_sweeping_order(self, state, step_number, order, Iden,
                          forward = True):
        t = 0
        while t < step_number:
            for oper in order:
                state = self.ED_sweep(oper, Iden, state, self.Id, forward)
                forward = not forward
            t+=1
        return state
    
    # Exact time evolution by matrix multiplication
    def exact_time(self, step_number):
        state = self.i_state
        e, v = np.linalg.eigh(self.H)
        
        DO = np.matmul(np.diag(np.exp(-self.dt * e)), np.linalg.inv(v))
        expH = np.matmul(v, DO)
        for t in range(step_number):
            state = np.dot(expH,state)
            state = state / np.sqrt(sum(state ** 2))
        self.E_evol = np.dot(state, np.dot(self.H, state))
        self.evol_state = state
        return
    
    # Exact ground state by diagonalization with a Lanczos algorithm. Stores
    # the three lowest states as class properties and returns the ground state
    # energy and ground state
    def exact_GS(self):
        self.elist, self.GSl = sp_linalg.eigsh(self.H, k=3, v0 = self.i_state,
                                         which="SA")
        self.E_GS = min(self.elist)
        self.GS = self.GSl[: 2**self.N, 0]
        return self.E_GS, self.GS
    
    # Retrieves the energy of a state in the format of states given by method
    # exact_GS (vector).
    def get_ener(self,state):
        E = []
        Id = np.eye(self.d)
        Iden = list(itertools.repeat(Id, self.N-1))
        for k in range(self.N - 1):
            Iden[k] = self.Hchain[k]
            Enop = Nkron(Iden)
            E.append(np.dot(state, np.dot(Enop, state)))
            Iden[k] = Id
        return E
    
    # Performs a Trotter sweep and is called by function ED_sweeping_order().
    # Returns the time evolved state.
    def ED_sweep(self, Ut, Iden, Istate, I, forw=True):
        start = 0
        end = self.N - 1
        inc = 1
        if not forw:
            start = self.N - 2
            end = -1
            inc = -1
        for i in range(start, end, inc):
            U = spkron([sp.eye(self.d ** i), Ut[i],
                        sp.eye(self.d ** (self.N-i-2))
                        ]
                      )
            Istate = U.dot(Istate)
            Istate = Istate / np.sqrt(sum(Istate ** 2))
        return Istate
    
    # Calculates the correlator of single site operator op1 at i and op2 at j
    # in a state psi and returns the result.
    def ED_correl(self, psi, op1, op2, i, j):
        op1 = spkron([sp.eye(self.d ** i), op1,
                      sp.eye(self.d ** (self.N-i-1))
                      ]
                    )
        op2 = spkron([sp.eye(self.d ** j), op2,
                      sp.eye(self.d ** (self.N-j-1))
                      ]
                    )
        phi = op1.dot(op2.dot(psi))
        corr = np.dot(psi, phi)
        return corr
    
    # Using a state in vector format calculates the lambda tensor in MPS
    # formalism at bond i.
    def get_lam(self, state, i):
        phi = np.reshape(state, (self.d ** i, self.d ** (self.N-i)))
        lam = np.linalg.svd(phi, full_matrices=False, compute_uv=False)
        return lam

# Kronecker product of N matrices
def Nkron(matlist):
    if matlist == []:
        return 1
    else:
        return np.kron(matlist[0], Nkron(matlist[1 :]))

# Kronecker product of N matrices in sparse matrix format (coordinate)
def spkron(smatlist):
    if smatlist == []:
        return 1
    else:
        return sp.kron(smatlist[0], spkron(smatlist[1 :]))
