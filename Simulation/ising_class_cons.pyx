#!/usr/bin/python
#cython: boundscheck=False,wraparound=False,cdivision=True
import numpy as np
from libc.stdlib cimport rand, srand, RAND_MAX
from libc.time cimport time
cimport numpy
srand(time(NULL))

cdef int[:,:] adjacent = np.zeros((4,2),dtype = np.intc)
adjacent[0,0] = 1
adjacent[1,0] = -1
adjacent[2,1] = 1
adjacent[3,1] = -1  
cdef class ising:
    def __init__(self,int dim,int[:,:] stuck,float temp,float comp):
        self.D = dim 
        self.L = np.zeros((self.D,self.D),dtype = np.intc)
        self.B = np.zeros((dim,dim),dtype = np.intc)
        self.update_stuck(stuck)
        self.m = comp
        self.init_lattice()
        self.T = temp
        self.probs = np.zeros(17,dtype = np.double)
        self.make_probs()
        self.flips = 0
    cdef float crandom(ising self):
        return <float>rand()/<float>RAND_MAX
    cdef int crandint(ising self, int lower, int upper):
        return (rand() % (upper - lower + 1)) + lower

    cpdef init_lattice(ising self):
        cdef int total_up = int(self.m*(self.D**2))
        cdef int total_down = int((1-self.m)*(self.D**2))
        cdef int i = 0
        cdef int j,s
        for i in range(self.D):
            j = 0
            for j in range(self.D):
                if self.B[i,j] != 0:
                    self.L[i,j] = -1
                    total_down -= 1
        for i in range(self.D):
            j = 0
            for j in range(self.D):
                if self.L[i,j] == 0:
                    if total_up != 0 and total_down != 0:
                        s = np.random.choice([-1,1])
                        self.L[i,j] = s
                        if s == -1:
                            total_down -= 1
                        else: 
                            total_up -= 1
                    else:
                        if total_down == 0:
                            self.L[i,j] = 1
                            total_up -= 1
                        if total_up == 0:
                            self.L[i,j] = -1
                            total_down -= 1
        print total_up, total_down, np.sum(self.B)
        return 

    #makes probability array 
    cpdef make_probs(ising self):
        Jc = np.log(np.sqrt(2) + 1) /2
        cdef int i=0;
        for i in range(17):
            self.probs[i] = np.exp((Jc*2*(-i))/self.T)

    #updates stuck, storing points of the bound polymers
    cpdef update_stuck(ising self,int[:,:] sticky):
        cdef int i = 0;
        cdef int j;
        self.B = np.zeros((self.D,self.D),dtype=np.intc)
        for i in range(self.D):
            for j in range(self.D):
                self.B[i,j] = sticky[i,j]
    cpdef set_temp(self,float t):
        self.T = t
        self.make_probs()
    # returns spin lattice 
    cpdef int[:,:] get_spins(ising self):
        return self.L
    cpdef set_spins(ising self, spin):
        self.L = spin
    #sweeps through the lattice
    cpdef sweep_lattice(ising self):
        cdef int idx1,idx2;
        cdef int i = 0;
        cdef int N = self.D**2   
        cdef int[:] prop_flip;
        sites = np.arange(N,dtype = np.intc)
        prop_flip = np.random.permutation(sites)
        for i in range(N):
            idx1 = prop_flip[i]/self.D
            idx2 = prop_flip[i]%self.D
            if self.B[idx1,idx2] == 0:
                self.kawasaki(idx1,idx2)
    cpdef int get_flips(ising self):
        return self.flips
    cpdef kawasaki(ising self,int idx1, int idx2):
        cdef float energy = 0
        cdef int ridx = self.crandint(0,3)
        cdef int jdx1 = (idx1 + adjacent[ridx,0])%self.D
        cdef int jdx2 = (idx2 + adjacent[ridx,1])%self.D
        if jdx1 == -1:
            jdx1 = self.D - 1
        if jdx2 == -1:
            jdx2 = self.D - 1
        cdef int spin_i = self.L[idx1,idx2]
        cdef int spin_j = self.L[jdx1,jdx2]
        cdef int nnx,nny,si,sj,k,Eswitch;
        if spin_i != spin_j:
            if self.B[jdx1,jdx2] != 0: 
                return 0
            si,sj,k = 0,0,0
            for k in range(4):
                nnx = adjacent[k,0]
                nny = adjacent[k,1]
                nidx1 = (idx1 + nnx)%self.D
                nidx2 = (idx2 + nny)%self.D
                if nidx1 == -1:
                    nidx1 = self.D - 1
                if nidx2 == -1:
                    nidx2 = self.D - 1
                njdx1 = (jdx1 + nnx)%self.D
                njdx2 = (jdx2 + nny)%self.D
                if njdx1 == -1:
                    njdx1 = self.D - 1
                if njdx2 == -1:
                    njdx2 = self.D - 1
                nspin_i = self.L[nidx1,nidx2]
                nspin_j = self.L[njdx1,njdx2]
                si += (spin_i*nspin_i)
                sj += (spin_j*nspin_j)
            Eswitch = (-(si+1)) + (-(sj+1))
            if Eswitch >= 0:
                self.L[idx1,idx2] = spin_j
                self.L[jdx1,jdx2] = spin_i
                return 1
            elif self.probs[-Eswitch] > self.crandom():
                self.L[idx1,idx2] = spin_j
                self.L[jdx1,jdx2] = spin_i
                return 1
            else:
                return 0
        else:
            return 0
                
    def print_spins(self,fn,iter):
        to_print = ""
        for i in range(self.D):
            for j in range(self.D):
                to_print += str(self.L[i,j]) + "\t"
        fn.write("%d\t%s\n" % (iter,to_print))
        return
