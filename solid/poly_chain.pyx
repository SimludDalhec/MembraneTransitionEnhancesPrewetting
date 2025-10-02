#!/usr/bin/python
#cython boundscheck=False,wraparound=False
import numpy as np
from libc.stdlib cimport rand, srand, RAND_MAX
from libc.time cimport time
cimport numpy as np
srand(time(NULL))
## Polymer Class for a chain
## Can move
## 3 Dimensional
## Observes "Sticky" Boundary on X = 0 plane. Periodic boundaries for Z 
cdef int[:,:] moves = np.zeros((6,3),dtype =np.intc)
moves[0,0] =1
moves[1,0] = -1
moves[2,1] = 1
moves[3,1] = -1
moves[4,2] = 1
moves[5,2] = -1 

cdef int[:] unallowed = np.zeros(6,dtype =np.intc)
unallowed[0] = 1
unallowed[1] = 0
unallowed[2] = 3
unallowed[3] = 2
unallowed[4] = 5
unallowed[5] = 4 


cdef class poly_chain:
    #initializes polymer
    def __init__(self,dim1,dim2,N,s):
        self.init_x = np.random.randint(dim2)
        self.init_y = np.random.randint(dim1)
        self.init_z = np.random.randint(dim1)
        self.length = N
        self.pos =[[self.init_x,self.init_y,self.init_z] for i in range(N)]# np.zeros((N,3),dtype=np.intc)
        self.D1 = dim1
        self.D2 = dim2
        self.sys = s
        # Set positions_
        self._set_positions_fast()
    #prints positions
    cpdef get_pos(poly_chain self):
        return self.pos[:]
    cpdef int get_sys(poly_chain self):
        return self.sys
    cpdef set_sys(poly_chain self,int i):
        self.sys = i
    cpdef set_pos(poly_chain self,p):
        cdef int i = 0
        cdef int N = len(p)
        for i in range(N):
            self.pos[i][0] = p[i][0]
            self.pos[i][1] = p[i][1]
            self.pos[i][2] = p[i][2]
    cpdef set_pos_local(poly_chain self,pF,pI,int idx_F,int idx_I):
        #pos_init = self.pos
        if idx_I == 0:
            #print "front to end",self.pos
            for k in range(idx_I,self.length-1):
                self.pos[k][0] = pI[k+1][0]
                self.pos[k][1] = pI[k+1][1]
                self.pos[k][2] = pI[k+1][2]
            self.pos[idx_F][0] = pF[0]
            self.pos[idx_F][1] = pF[1]
            self.pos[idx_F][2] = pF[2]
            #print "front to end", self.pos
            return
        elif idx_I == self.length-1:
            #print "End to front",self.pos
            for k in range(idx_I,-1,-1):
                self.pos[k][0] = pI[k-1][0]
                self.pos[k][1] = pI[k-1][1]
                self.pos[k][2] = pI[k-1][2]
            self.pos[idx_F][0] = pF[0]
            self.pos[idx_F][1] = pF[1]
            self.pos[idx_F][2] = pF[2]
            #print "end to front",self.pos
            return
        elif idx_I == idx_F:
            self.pos[idx_F][0] = pF[0]
            self.pos[idx_F][1] = pF[1]
            self.pos[idx_F][2] = pF[2]
            return
        else:
            return
        
    cdef float crandom(poly_chain self):
        return <float>rand()/<float>RAND_MAX
    cdef int crandint(poly_chain self,int lower, int upper):
        return (rand() % (upper - lower + 1)) + lower    

    #sets random positions to start, starts as straight line     
    cdef _set_positions_fast(self):
        cdef int i = 0
        dir = np.random.randint(3)
        if(dir == 1):
            if np.random.randint(2) == 0:
                self.init_x = np.random.randint(0,self.D2-self.length-1)                
                for i in range(self.length):
                    self.pos[i][0] = self.init_x+i
            else:
                self.init_x = np.random.randint(self.length+1,self.D2-1)                
                for i in range(self.length):
                    self.pos[i][0] = self.init_x-i
        elif dir == 2:
            if(np.random.randint(2)) == 0:
                for i in range(self.length):
                    self.pos[i][1] = (self.init_y+i) % self.D1
            else:
                for i in range(self.length):
                    self.pos[i][1] = (self.init_y-i) % self.D1

        else:
            if(np.random.randint(2)) == 0:
                for i in range(self.length):
                    self.pos[i][2] = (self.init_z+i) % self.D1
            else:
                for i in range(self.length):
                    self.pos[i][2] = (self.init_z-i) % self.D1
        
        #moves head-to-tail motion
    # next step for improvement would be making these all memoryviews? 
    # or reworking the output to give pos_from, pos replace, and use that in simulate
    cdef move_reptamer(poly_chain self):
        cdef int direction,ht,j,d,mx,my,mz;
        d = self.D1
        ht = self.crandint(0,1)        
        direction = self.crandint(0,5)
        #print ht,direction
        prev = self.pos
        cdef int N = self.length-1
        mx = moves[direction,0]
        my = moves[direction,1]
        mz = moves[direction,2]
        if ht == 1:
            move_final = [(prev[N][0] + mx),(prev[N][1] + my) % d,(prev[N][2] + mz) % d]
            move_initial = [self.pos[0][0],self.pos[0][1],self.pos[0][2]]
            #print move_final,move_initial,N,0,"rept th"
            return move_final,move_initial,N,0
        else:
            move_final = [(prev[0][0]+mx),(prev[0][1]+my)%d,(prev[0][2]+mz)%d]
            move_initial = [self.pos[N][0],self.pos[N][1],self.pos[N][2]]
            #print move_final,move_initial,0,N,"rept ht"
            return move_final,move_initial,0,N

    #Pops a corner: from any orientation where b(i-1) != bi+1, b(i)
    cdef move_kink(poly_chain self):
        cdef int loc = self.crandint(1,self.length-2)
        cdef int pX,pY,pZ,pXn,pYn,pZn,pXp,pYp,pZp
        pos_init = self.pos[loc]
        pX = pos_init[0]
        pY = pos_init[1]
        pZ = pos_init[2]
        nextP=self.pos[loc+1]
        pXn = nextP[0]
        pYn = nextP[1]
        pZn = nextP[2]
        prev=self.pos[loc-1]
        pXp = prev[0]
        pYp = prev[1]
        pZp = prev[2]        
        bVec1 = ((pX - pXn),(pY-pYn)%self.D1,(pZ-pZn)%self.D1)
        bVec2 = ((pX - pXp),(pY-pYp)%self.D1,(pZ-pZp)%self.D1)
        move = [(pX - (bVec1[0] + bVec2[0])), (pY - (bVec1[1] +bVec2[1]))%self.D1,(pZ - (bVec1[2]+bVec2[2]))%self.D1]
        #print "move kink"
        #print move,pos_init,loc,loc,"kink"
        return move,pos_init,loc,loc

    cpdef make_move(poly_chain self):
        #if self.crandom() < 0.5:
        #    return self.move_kink()        
        #else:
        return self.move_reptamer()
    
    cpdef make_move_non_local(poly_chain self):
        ### choose a random position
        cdef int pX,pY,pZ,pXn,pYn,pZn,pXp,pYp,pZp
        cdef int idx_init = self.crandint(0,self.length-1)        
        pos_init = self.pos[idx_init]
        pX = pos_init[0]
        pY = pos_init[1]
        pZ = pos_init[2]
        ### translate D spaces in x, D spaces in Y dimennsions
        cdef int transY = self.crandint(0,self.D1-1)
        cdef int transZ = self.crandint(0,self.D1-1)
        cdef int i = 0
        pos_new = []
        #### set new position, Z position is the same
        for i in range(self.length):
            pos_newY = (self.pos[i][1] + transY)%self.D1
            pos_newZ = (self.pos[i][2] + transZ)%self.D1
            pos_new.append((self.pos[i][0],pos_newY,pos_newZ))
            #print(pos_new)
        return pos_new
