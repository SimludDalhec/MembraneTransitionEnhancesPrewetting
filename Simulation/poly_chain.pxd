#!/usr/bin/python
import cython
cdef class poly_chain:
     cdef int init_x
     cdef int init_y
     cdef int init_z
     cdef int sys
     cdef int D1
     cdef int D2
     cdef int length
     cdef object pos
     cdef _set_positions_fast(poly_chain self)
     cdef int crandint(poly_chain self,int lower, int upper)
     cdef float crandom(poly_chain self)
     cpdef get_pos(poly_chain self)
     cpdef set_pos_local(poly_chain self,pF,pI,int idx_F,int idx_I)
     cpdef int get_sys(poly_chain self)
     cpdef set_sys(poly_chain self, int i)
     cpdef set_pos(poly_chain self, p)
     cdef move_reptamer(poly_chain self)
     cdef move_kink(poly_chain self)
     cpdef make_move(poly_chain self)
     cpdef make_move_non_local(poly_chain self)
