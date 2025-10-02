#!/usr/bin/python
cdef class ising:
     cdef int D
     cdef float m
     cdef int[:,:] stuck
     cdef double[:] probs
     cdef float T
     cdef int[:,:] L
     cdef int[:,:] B
     cdef int flips;
     cdef int crandint(ising self,int lower, int upper)
     cdef float crandom(ising self)

     cpdef int get_flips(ising self)
     cpdef init_lattice(ising self)
     cpdef make_probs(ising self)
     cpdef set_spins(ising self,spin)
     cpdef update_stuck(self,int[:,:] sticky)
     cpdef set_temp(ising self, float t)
     cpdef int[:,:] get_spins(ising self)  
     cpdef sweep_lattice(ising self)
     cpdef kawasaki(ising self, int idx1, int idx2)
     
