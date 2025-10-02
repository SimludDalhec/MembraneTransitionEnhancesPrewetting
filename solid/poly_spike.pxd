#!/usr/bin/python
cdef class poly_spike:
    cdef int D
    cdef int length
    cdef int[:,:] pos
    cpdef get_pos(poly_spike self)
    cpdef set_pos(poly_spike self, int y, int z)
    cdef (int,int) move_translate(poly_spike self)
    cpdef make_move(poly_spike self)
    cdef int crandint(poly_spike self,int lower,int upper)
