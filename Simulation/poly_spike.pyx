#!/usr/bin/python
import numpy as np
from libc.stdlib cimport rand, RAND_MAX
cimport numpy as np

cdef int[:,:] moves = np.zeros((4,2),dtype =np.intc)
moves[0,0] =1
moves[1,0] = -1
moves[2,1] = 1
moves[3,1] = -1
cdef class poly_spike:
    def __init__(self,int dim,int N):
        self.length = N
        self.D = dim
        # Iniitialize and set positions
        self.pos = np.zeros((self.length,3),dtype=np.intc)
        self.set_pos(rand()%dim,rand()%dim)
    cpdef get_pos(poly_spike self):
        return self.pos[0,1],self.pos[0,2]
    cpdef set_pos(poly_spike self,int y, int z):
        cdef int i = 0
        for i in range(self.length):
            self.pos[i,0] = i
            self.pos[i,1] = y
            self.pos[i,2] = z
        return
    cdef int crandint(poly_spike self,int lower, int upper):
        return (rand() % (upper - lower + 1)) + lower
    cdef (int,int) move_translate(poly_spike self):
        cdef int direction = self.crandint(0,3)
        cdef int yPos = (self.pos[0,1]+moves[direction,0])%self.D
        cdef int zPos = (self.pos[0,2]+moves[direction,1])%self.D
        return yPos,zPos
    cpdef make_move(poly_spike self):
        return self.move_translate()


