#!/usr/bin/python
#cython: profile=True
import numpy as np
from networkx.algorithms.components.connected import connected_components
import random
import networkx as nx
import os,sys
from libc.stdlib cimport rand, RAND_MAX,srand
from libc.time cimport time
cimport cython
cimport numpy as np 
cimport poly_spike as ps 
cimport ising_class_cons as I
cimport poly_chain as pc

##################################################
############## CONSTANTS #########################
##################################################
cdef int print_interval = 50000
#cdef int iters = 50000000
cdef int iters =  10000000
cdef int D1 = 30;                                   # maximum length of box from surfae
cdef int L = 64;                                  # Dimensions of surface
cdef int N_res = 150;                               # Polymers in reservoir
cdef int N_restot = 2*N_res;                               # Polymers in reservoir
cdef int N_sys = 25;                                # initial polys in system
cdef double[:] p_cluster_arr = make_cluster_probs() # prob to propose cluster move
cdef double[:,:,:,:] prob;#array of switching probabilities
out_path = "/home/mnr29/project/PrewettingExperiments/JobsShuffleFigS4_2/"
# these will be much different than in the other sims
cdef float J_nn = 0.1
cdef int Lbulk,Lt1,Lt2;
Lbulk,Lt1 = 10,5
## For checking nearest neighbors 
cdef int[:,:] nn = np.array([[1,0,0],[-1,0,0],[0,-1,0],[0,1,0],[0,0,1],[0,0,-1]],dtype=np.intc) #nearest neighbor coordinates
cdef int[:,:] nn_tether = np.array([[0,1],[0,-1],[1,0],[-1,0]],dtype =np.intc)
cdef int[:,:] nn_lo = np.array([[1,0,0],[0,-1,0],[0,1,0],[0,0,1],[0,0,-1]],dtype=np.intc) #nearest neighbor coordinates
cdef int[:,:] nn_hi = np.array([[-1,0,0],[0,-1,0],[0,1,0],[0,0,1],[0,0,-1]],dtype=np.intc) #nearest neighbor coordinates 
##################################################?

#########################################################

#### Main loop
# runs through, first decreasing bulk temperature, then surface bulk coupling, then ising temp
def main(float cT1, float Jbulk, float chem_potent1, float chem_potent2, float comp, float Ti, int fidx, double res_speedup):
    cdef int Nt1 = int((L**2.0)*cT1)
    srand(int(time(NULL)*(-chem_potent1+fidx)))
    dir_name = out_path + "d{0}-{1}_n2t{2}_l{3}_jnn{4}".format(D1,L,Nt1,Lbulk,J_nn)
    if os.path.isdir(dir_name) == False:
        os.mkdir(dir_name)
    # out directory containst all of the thermodynamic paramteres
    fn = dir_name+"/Ti{0}_Jb{1}_u{2}_u2{3}_M{4}_idx{5}_speedup{6:.3f}".format(Ti,Jbulk,chem_potent1,chem_potent2,comp,fidx,res_speedup)
    fn_i = dir_name+"/Ti{0}_Jb{1}_u{2}_u2{3}_M{4}_idx{5}_speedup{6:.3f}_ising.txt".format(Ti,Jbulk,chem_potent1,chem_potent2,comp,fidx,res_speedup)
    cdef int[:,:] tetherLattice  = np.zeros((L,L),dtype = np.intc)    
    switch_probs = make_switch_probs(chem_potent1,chem_potent2, J_nn)
    probs = make_probs(Jbulk,J_nn)
    bulk_state_sys,bulk_state_res,tether_state = init_state(Nt1) 
    for k in range(Nt1):
        pos = tether_state[k].get_pos()
        g,f = pos
        tetherLattice[g,f] = 1
    ising = I.ising(L,tetherLattice,Ti,comp)    
    cdef int[:,:,:,:] lattice,lattice_res;
    lattice= make_lattice(bulk_state_sys,0)
    lattice_res = make_lattice(bulk_state_res,1)
    sim(bulk_state_sys,bulk_state_res,tether_state,fn,fn_i,ising,lattice,lattice_res,Nt1,probs,switch_probs,tetherLattice,res_speedup)

    return

#########################
####### FUNCTIONS #######
#########################

#Initialized each polymer to straight lines
# puts even amount of red-blue polymers in reservoir and system
def init_state(Nt1):
    polys1,polys2 = [[pc.poly_chain(L,D1,Lbulk,0)],[pc.poly_chain(L,D1,Lbulk,0)]]
    polys1_res,polys2_res = [[pc.poly_chain(L,D1,Lbulk,0)],[pc.poly_chain(L,D1,Lbulk,0)]]
    # Placing bulk polymers in the system
    a = True    
    while len(polys1) < N_sys:
        poly = pc.poly_chain(L,D1,Lbulk,0)
        a = check_intersect_polys(poly,polys1)
        if a == False:
            polys1.append(poly)
    a = True
    while len(polys2) < N_sys:
        poly = pc.poly_chain(L,D1,Lbulk,0)   
        a = check_intersect_polys(poly,polys2)
        if a == False:
            polys2.append(poly)

    # Placing bulk polymers in the particle reservoir
    a = True 
    while len(polys1_res) < N_res:
        poly = pc.poly_chain(L,D1,Lbulk,1)
        a = check_intersect_polys(poly,polys1_res)
        if a == False:
            polys1_res.append(poly)
    a = True 
    while len(polys2_res) < N_res:
        poly = pc.poly_chain(L,D1,Lbulk,1)   
        a = check_intersect_polys(poly,polys2_res)
        if a == False:
            polys2_res.append(poly)

    # Placing tethers 
    tethers1 = [ps.poly_spike(L,Lt1)]
    a = True 
    while len(tethers1) < Nt1:
        tether= ps.poly_spike(L,Lt1) 
        a = check_intersect_tethers(tether,[tethers1])
        if a == False:
            a =  check_intersect_poly_tethers(tether,polys1)
            if a == False:
                tethers1.append(tether)

    bulk_state_sys = [polys1,polys2]
    bulk_state_res = [polys1_res,polys2_res]
    tether_state = tethers1
    return bulk_state_sys,bulk_state_res,tether_state

def check_intersect_polys(poly,state):
    pos2 = tuple(map(tuple,poly.get_pos()))
    for p in state:
        pos_occ = tuple(map(tuple,p.get_pos()))
        inter = len(list(set(pos2) & set(pos_occ)))
        if inter > 0:
            return True
    return False
def check_intersect_tethers(tether,state):
    tethers1 = state[0]
    tx,ty = tether.get_pos()
    for t1 in tethers1:
        tx1,ty1 = t1.get_pos()
        if tx == tx1 and ty == ty1:
            return True
    return False

def check_intersect_poly_tethers(tether,polys):
    tx,ty = tether.get_pos()
    for p in polys:
        posi = p.get_pos()
        for i in posi:
            px,py,pz = i
            if px < Lt1:
                if py == tx and pz == ty:
                    return True
    return False


# simulates a system with the particle reservoir
# Moves for each polymer are proposed then acc/rej 
# Ising model is sweeped over once per sweep through polymers
# A poisson number of polymers are selected to be switched from sys -> res/res -> switch. 
cdef sim(bulk_state_sys,bulk_state_res,tether_state,fn,fn_i,ising, int[:,:,:,:] lattice,int[:,:,:,:] lattice_res,int Nt1, double[:,:] probs, double[:,:,:] switch_probs, int[:,:] tetherLattice,double res_speedup):
    cdef int i,j,a,p_total,b,n,s,Nsys,k,ridx,finalY,finalZ,initY,initZ,px,py
    cdef int[:] pseq;
    cdef int[:] pseq_res = np.arange(N_restot,dtype=np.intc)
    cdef int pnum,ptype
    cdef double p_swap;
    b,n,s,i,j= 0,0,0,0,0
    cdef int res_iterations = int(max(res_speedup,1))
    cdef int poly_iterations = min(int(N_restot*res_speedup),N_restot)
    print(res_iterations,poly_iterations)
    f_i = open(fn_i,'w')
    f_state = open(fn+"_pos.txt",'w')
    for i in range(iters):
        ### update the system state, tethers and bulk ####
        N1 = len(bulk_state_sys[0])
        N2 = len(bulk_state_sys[1])
        Nsys = N1 + N2
        p_total = Nsys + Nt1
        pseq = np.arange(p_total,dtype=np.intc)
        fisher_yates_shuffle(pseq,p_total)
        a = 0;
        for a in range(p_total):
            bulk_tether,pnum,ptype = idx_to_idx_sys(pseq[a],N1,N2,Nt1)
            if bulk_tether == 1: # Moving a bulk polymer                
                b,n,s = 0,0,0
                final,init,fIdx,iIdx = bulk_state_sys[ptype][pnum].make_move()
                poly_pos = bulk_state_sys[ptype][pnum].get_pos()
                b,n = get_move_energy(final,init,lattice,tetherLattice,ptype)
                if b != -100:
                    if b > 0 and n > 0:
                        lattice = update_lattice(final,init,lattice,ptype)
                        bulk_state_sys[ptype][pnum].set_pos_local(final,poly_pos,fIdx,iIdx)
                    elif probs[b+15,n+35] > crandom():
                        lattice = update_lattice(final,init,lattice,ptype)
                        bulk_state_sys[ptype][pnum].set_pos_local(final,poly_pos,fIdx,iIdx)                       
            else:
                b,s,n = 0,0,0
                finalY,finalZ = tether_state[pnum].make_move()
                initY,initZ = tether_state[pnum].get_pos()
                s,n = get_move_energy_tether(finalY,finalZ,initY,initZ,lattice,tetherLattice,ising)
                if s != -100:
                    if s > 0 and n > 0:
                        tetherLattice[initY,initZ] = 0
                        tetherLattice[finalY,finalZ] = 1
                        tether_state[pnum].set_pos(finalY,finalZ)
                    elif probs[s+15,n+35] > crandom():
                        tetherLattice[initY,initZ] = 0
                        tetherLattice[finalY,finalZ] = 1
                        tether_state[pnum].set_pos(finalY,finalZ)
        ####  reservoir ###
        j=0
        for j in range(res_iterations):
            fisher_yates_shuffle(pseq_res,N_restot)
            a = 0 
            for a in range(poly_iterations):
                pnum,ptype = idx_to_idx_res(pseq_res[a])
                final,init,fIdx,iIdx = bulk_state_res[ptype][pnum].make_move()
                poly_pos = bulk_state_res[ptype][pnum].get_pos()
                b = get_move_energy_res(final,init,lattice_res,ptype)
                if b != 0:                        
                    lattice_res = update_lattice(final,init,lattice_res,ptype)
                    bulk_state_res[ptype][pnum].set_pos_local(final,poly_pos,fIdx,iIdx)
           
        ising.update_stuck(tetherLattice)
        ising.sweep_lattice()

        ## PROPOSE SWITCH per sweep###
        N1 = len(bulk_state_sys[0])+len(bulk_state_res[0])
        N2 = len(bulk_state_sys[1])+len(bulk_state_res[1])
        p_swap = (float(N1))/(2*(N_sys+N_res))
        if crandom() < p_swap:
            b = 0 
            num = crandint(0,N1-1)
            idx,sys = idx_to_idx_swap(num)
            if sys == 0:
                pos = bulk_state_sys[0][idx].get_pos()
                add_rem = prop_switch(0,sys,pos,lattice,lattice_res,tetherLattice,switch_probs)
                if add_rem == -1:
                    lattice = update_lattice_exchange(pos,0,add_rem,lattice)
                    bulk_state_sys[0].remove(bulk_state_sys[0][idx])
            else:
                pos = bulk_state_res[0][idx].get_pos()
                add_rem = prop_switch(0,sys,pos,lattice,lattice_res,tetherLattice,switch_probs)
                if add_rem == 1:
                    polycopy = pc.poly_chain(L,D1,Lbulk,0)
                    polycopy.set_pos(pos)
                    lattice = update_lattice_exchange(pos,0,add_rem,lattice)
                    bulk_state_sys[0].append(polycopy)
                ### Shuffle the polymer positions
                pos_nonlocal = []
                while b == 0:
                    pos_nonlocal = bulk_state_res[0][idx].make_move_non_local()
                    b = check_pos_non_local(pos_nonlocal,lattice_res,0)
                lattice_res=update_lattice_non_local(pos_nonlocal,pos,lattice_res,0)
                bulk_state_res[0][idx].set_pos(pos_nonlocal)


        ## PROPOSE SWITCH per sweep###
        N1 = len(bulk_state_sys[0])+len(bulk_state_res[0])
        N2 = len(bulk_state_sys[1])+len(bulk_state_res[1])
        p_swap = (float(N2))/(2*(N_sys+N_res))
        if crandom() < p_swap:
            num = crandint(0,N2-1)
            idx,sys = idx_to_idx_swap(num)
            b=0
            if sys == 0:
                pos = bulk_state_sys[1][idx].get_pos()
                add_rem = prop_switch(1,sys,pos,lattice,lattice_res,tetherLattice,switch_probs)
                if add_rem == -1:
                    lattice = update_lattice_exchange(pos,1,add_rem,lattice)
                    bulk_state_sys[1].remove(bulk_state_sys[1][idx])
            else:
                pos = bulk_state_res[1][idx].get_pos()
                add_rem = prop_switch(1,sys,pos,lattice,lattice_res,tetherLattice,switch_probs)
                if add_rem == 1:
                    polycopy = pc.poly_chain(L,D1,Lbulk,0)
                    polycopy.set_pos(pos)
                    lattice = update_lattice_exchange(pos,1,add_rem,lattice)
                    bulk_state_sys[1].append(polycopy)
                ### Shuffle the polymer positions
                pos_nonlocal = []
                while b == 0:
                    pos_nonlocal = bulk_state_res[1][idx].make_move_non_local()
                    b = check_pos_non_local(pos_nonlocal,lattice_res,1)
                lattice_res=update_lattice_non_local(pos_nonlocal,pos,lattice_res,1)
                bulk_state_res[1][idx].set_pos(pos_nonlocal)

        # PRINTING
        if i % print_interval == 0:
            #print('Reservoir size,',len(bulk_state_res[0]),len(bulk_state_res[1]),np.sum(lattice_res[0])/Lbulk,np.sum(lattice_res[1])/Lbulk)
            #print('Reservoir state,',np.mean(bulk_state_res[0][0].get_pos(),axis=0),np.mean(bulk_state_res[1][0].get_pos(),axis=0))
            #print('System size,',len(bulk_state_sys[0]),len(bulk_state_sys[1]),np.sum(lattice[0])/Lbulk,np.sum(lattice[1])/Lbulk)
            ising.print_spins(f_i,i)
            print_state(bulk_state_sys,tether_state,f_state,i)
    f_i.close()
    f_state.close()
    return

## TRANSLATES POLYMER COUNT IDX TO PE-TYPE IDX
cdef (int,int,int) idx_to_idx_sys(int a,int N1,int N2,int Nt1):
    cdef int pnum,ptype,bulk;
    pnum,ptype,bulk = 0,0,0
    if a >= N2 + N1:
        ptype = 0
        bulk = 0
        pnum = a - (N2 + N1)
        return bulk,pnum,ptype
    elif a >= N1 and a < N2 + N1:
        ptype = 1
        bulk = 1
        pnum = a - N1
        return bulk,pnum,ptype
    else:
        ptype = 0
        bulk = 1
        pnum = a
        return bulk,pnum,ptype

## TRANSLATES POLYMER COUNT IDX TO PE-TYPE IDX
cdef (int,int) idx_to_idx_res(int a):
    if a >= N_res:
        return a-N_res,0
    else:
        return a, 1


## TRANSLATES POLYMER COUNT IDX TO PE-TYPE IDX
cdef (int,int) idx_to_idx_swap(int a):
    if a >= N_res:
        return a - N_res,0
    else:
        return a,1

## energy of each move. Rejecting (-100) if self/type overlaps or violates boundary conditions
# otherwise return change in bonds form initial to final position
@cython.boundscheck(False)
@cython.wraparound(False)
cdef (int,int) get_move_energy_tether(int finalY, int finalZ,int initialY,int initialZ, int[:,:,:,:] lattice, int[:,:] tether_lattice, ising):
    cdef int i,x,y,z,xi,yi,zi,bonds_tether,nn;
    cdef int[:,:] isingLat;
    i,nn,bonds_tether = 0,0,0;
    isingLat = ising.get_spins()
    if isingLat[finalY,finalZ] != isingLat[initialY,initialZ]:
        return -100,-100
    if tether_lattice[finalY,finalZ] != 0:
        return -100,-100
    for i in range(Lt1):
        if lattice[0,i,finalY,finalZ] != 0:
            return -100,-100
        bonds_tether += (lattice[1,i,finalY,finalZ] -  lattice[1,i,initialY,initialZ])
    nn = check_nn_tether(finalY, initialY , finalZ,initialZ, lattice,Lt1)
    return bonds_tether,nn

## energy of each move. Rejecting (-100) if self/type overlaps or violates boundary conditions
# otherwise return change in bonds form initial to final position
@cython.boundscheck(False)
@cython.wraparound(False)
cdef (int,int) get_move_energy(final,initial,int[:,:,:,:] lattice, int[:,:] tetherLattice, int polytype):
    cdef int x,y,z,xi,yi,zi
    #check that doens't overlap with type
    #print final,initial
    x = final[0]
    y = final[1]
    z = final[2]
    xi = initial[0]
    yi = initial[1]
    zi = initial[2]
    #print final,initial
    if lattice[polytype,x,y,z] != 0:
        return -100,0
    if x == xi and y == yi and z == zi:
        return -100,0
    if x > D1-1 or x < 0:
        return -100,0
    # In this scheme polytype 0 does not interact with tethers
    if polytype == 0:
        if x < Lt1:
            if tetherLattice[y,z] != 0:
                return -100,0
    #spikes cant move off of ising spins
    cdef int bonds_bulk = 0
    cdef int bonds_tether = 0
    cdef int nn =0
    #all polys are self avoiding but blue interacts with y/r and red only with blue 
    # This is where the interaction scheme can be changed if desired
    bonds_bulk = (lattice[1 - polytype,x,y,z] ) -  (lattice[1 - polytype,xi,yi,zi])
    nn = check_nn(x,y,z,xi,yi,zi,lattice,tetherLattice)
    if polytype == 1:
        if x < Lt1:
            bonds_bulk += tetherLattice[y,z]
        if xi < Lt1:
            bonds_bulk += -tetherLattice[yi,zi]
    return bonds_bulk,nn

## energy of each move. Rejecting (-100) if self/type overlaps or violates boundary conditions

## energy of each move. Rejecting (-100) if self/type overlaps or violates boundary conditions
# otherwise return change in bonds form initial to final position
@cython.boundscheck(False)
@cython.wraparound(False)
cdef int get_move_energy_res(final,initial,int[:,:,:,:] lattice,int polytype):
    cdef int x,y,z;
    x = final[0]
    y = final[1]
    z = final[2]
    if x == -1 or x == D1:
        return 0
    if lattice[polytype,x,y,z] != 0:
        return 0
    return 1


@cython.boundscheck(False)
@cython.wraparound(False)
cdef int check_pos_non_local(pos_nonlocal,int[:,:,:,:] lattice_res,int ptype):
    cdef int x,y,z,k;
    k=0
    for k  in range(Lbulk):
        x = pos_nonlocal[k][0]
        y = pos_nonlocal[k][1]
        z = pos_nonlocal[k][2]
        if lattice_res[ptype,x,y,z] != 0 or lattice_res[1-ptype,x,y,z] != 0:
            return 0
    return 1


## nn interactions satisfies metropolis criterion,no sort of overlaps allowed on either side.  
cdef int prop_switch(int ptype,int sys,pos,int[:,:,:,:] lattice,int[:,:,:,:] lattice_res, int[:,:] tetherLattice, double[:,:,:] switch_probs):
    cdef int i =0;
    cdef double p_change;
    cdef int resCollision,sysCollision
    #check that the spaces aren't occupied 
    for i in range(Lbulk):
        x,y,z = pos[i]
        resCollision = lattice_res[1,x,y,z]+ lattice_res[0,x,y,z]
        sysCollision = lattice[1,x,y,z] + lattice[0,x,y,z]
        if x < Lt1:
            if tetherLattice[y,z] != 0:
                return 0
        # no collisions on either side
        if sys == 0:
            if resCollision != 0:
                return 0
            if lattice[1 - ptype,x,y,z] == 1:
                return 0 
        else:
            if sysCollision != 0:
                #print('rejected sys collison' )
                return 0
            if lattice_res[1 - ptype,x,y,z] == 1:
                #print('rejected res collison' )
                return 0
    #get neighbor counts so can compute p(switch)
    n_count = get_neighbor_counts(pos,lattice,tetherLattice)
    p_change = switch_probs[1-sys,ptype,n_count]
    if sys == 1:
        if crandom() < p_change:
            return 1
        else:
            return 0
    #"Moving" from system into reservoir
    else:
        if crandom() < p_change:
            #print('rejected energy')
            return -1
        return 0

## Probs for chemcial potential
# idx1 = res/to sys (0_) or sys to res (1)
# idx2 = ammount of nearest neighbors
cdef double[:,:,:] make_switch_probs(float chem_potent1, float chem_potent2, float J_nn):
    cdef double[:,:,:] probs = np.zeros((2,2,Lbulk*6*2),dtype = np.double)
    n_cutoff = 10
    for i in range(Lbulk*6*2):
        #moving from reservoir to system
        if i > n_cutoff:
            probs[0,1,i] = 0
            probs[0,0,i] = 0
        else:
            energy = -chem_potent1+(-J_nn*i)
            probs[0,0,i] = np.exp(-energy,dtype =np.double)
            energy = -chem_potent2+(-J_nn*i)
            probs[0,1,i] = np.exp(-energy,dtype =np.double)
            #moving from system to reservoir Energy = u*(-1) + (J_nn*(-i))
            energy = chem_potent1+(J_nn*i)
            probs[1,0,i] = np.exp(-energy,dtype= np.double)
            energy = chem_potent2+(J_nn*i)
            probs[1,1,i] = np.exp(-energy,dtype= np.double)
    return probs 


#Getting full neighbor counts of a polymer 
#BCs make this somewhat long to get through
@cython.boundscheck(False)
@cython.wraparound(False)
cdef int get_neighbor_counts(poly,int[:,:,:,:] l, int[:,:] lt):
    cdef int neighbor_count,i,k,j
    neighbor_count,i,k,j = 0,0,0,0;
    cdef int px,py,pz,pxn,pyn,pzn;
    cdef int overlap = 0;
    #cdef int N = len(poly)
    for i in range(Lbulk):
        px,py,pz = poly[i] 
        if (px < D1 - 1 and px > 0):
            for j in range(6):
                p = nn[j]
                overlap,k=0,0;
                pxn = px + p[0]
                pyn = (py + p[1])% L
                pzn = (pz + p[2])% L
                for k in range(Lbulk):
                    if poly[k][0] == pxn and poly[k][1] == pyn and poly[k][2] == pzn:
                        overlap =1
                        break
                if overlap != 1:
                    neighbor_count += l[0,pxn,pyn,pzn] + l[1,pxn,pyn,pzn]
                    if pxn < Lt1 :
                        neighbor_count += lt[pyn,pzn]
        elif px == D1 - 1:
            for j in range(5):
                p = nn_hi[j]
                k,overlap =0,0
                pxn = px + p[0]
                pyn = (py + p[1])% L
                pzn = (pz + p[2])% L
                for k in range(Lbulk):
                    if poly[k][0] == pxn and poly[k][1] == pyn and poly[k][2] == pzn:
                        overlap =1
                        break
                if overlap != 1:
                    neighbor_count += l[0,pxn,pyn,pzn] + l[1,pxn,pyn,pzn]
        else:
            for j in range(5):
                p = nn_lo[j]
                overlap,k = 0,0
                pxn = px + p[0]
                pyn = (py + p[1])% L
                pzn = (pz + p[2])% L
                for k in range(Lbulk):
                    if poly[k][0] == pxn and poly[k][1] == pyn and poly[k][2] == pzn:
                        overlap =1
                        break
                if overlap != 1:
                    neighbor_count += l[0,pxn,pyn,pzn] + l[1,pxn,pyn,pzn]
                    if pxn < Lt1 :
                        neighbor_count += lt[pyn,pzn]
    return neighbor_count

### Check NN and add energies for a snake-like move
@cython.boundscheck(False)
@cython.wraparound(False)
cdef int check_nn(int x_f, int y_f, int z_f,int x_i, int y_i, int z_i,int[:,:,:,:] lattice, int[:,:] tetherLattice):
    cdef int i,n_f,n_i,xfn,yfn,zfn,xin,yin,zin;
    i,n_f,n_i = 0,0,0
    # initial neigboring stats
    if x_i == D1 - 1: 
        for i in range(5):
            n = nn_hi[i]
            xin = (x_i+n[0])
            yin = (y_i+n[1])%L
            zin = (z_i+n[2])%L
            n_i += lattice[0,xin,yin,zin] + lattice[1,xin,yin,zin] 
    elif x_i == 0:
        for i in range(5):
            n = nn_lo[i]
            xin = (x_i+n[0])
            yin = (y_i+n[1])%L
            zin = (z_i+n[2])%L
            n_i += lattice[0,xin,yin,zin] + lattice[1,xin,yin,zin]
            if xin < Lt1:
                n_i += tetherLattice[yin,zin]
    else:
        for i in range(6):
            n = nn[i]
            xin = (x_i+n[0])
            yin = (y_i+n[1])%L
            zin = (z_i+n[2])%L
            n_i += lattice[0,xin,yin,zin] + lattice[1,xin,yin,zin]
            if xin < Lt1:
                n_i += tetherLattice[yin,zin]
    # Final Neighboring states
    if x_f == D1 - 1:
        for i in range(5):
            n = nn_hi[i]
            xfn = (x_f+n[0])
            yfn = (y_f+n[1])%L
            zfn = (z_f+n[2])%L
            n_f += lattice[0,xfn,yfn,zfn] + lattice[1,xfn,yfn,zfn]
    elif x_f == 0:
        for i in range(5):
            n = nn_lo[i]
            xfn = (x_f+n[0])
            yfn = (y_f+n[1])%L
            zfn = (z_f+n[2])%L
            n_f += lattice[0,xfn,yfn,zfn] + lattice[1,xfn,yfn,zfn]
            if xfn < Lt1:
                n_f += tetherLattice[yfn,zfn]
    else:
        for i in range(6):
            n = nn[i]
            xfn = (x_f+n[0])
            yfn = (y_f+n[1])%L
            zfn = (z_f+n[2])%L
            n_f += lattice[0,xfn,yfn,zfn] + lattice[1,xfn,yfn,zfn]
            if xfn < Lt1:
                n_f += tetherLattice[yfn,zfn]
    return n_f - n_i

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int check_nn_tether(int y_f, int y_i, int z_f, int z_i, int[:,:,:,:] lattice, int N):
    cdef int i,j,n_f,n_i,yf,yi,zf,zi,ny,nz;
    i = 0 
    j = 0
    n_f = 0
    n_i = 0
    # initial neigboring stats
    for j in range(4):
        ny = nn_tether[j,0]
        nz = nn_tether[j,1]
        yi = (y_i + ny)%L
        yf = (y_f + ny)%L
        zi = (z_i + nz)%L
        zf = (z_f + nz)%L
        for i in range(N):
            n_i += lattice[0,i,yi,zi] + lattice[1,i,yi,zi]
            n_f += lattice[0,i,yf,zf] + lattice[1,i,yf,zf]
    return n_f - n_i


#makes a graph and returns the connected component that "a" is in 
# a slow but sure way to get connected groups of polymers
cdef make_graph(bulk_state,int N1,int N2,int a):
    g = nx.Graph()
    polys1,polys2 = bulk_state
    cdef int i = 0;
    cdef int j = 0;
    cdef int idx,jdx;
    ## N1 and N2 interact -- get these connections 
    for i in range(N1):
        if polys1[i].get_sys() == 0:
            g.add_node(i)
            g.add_edge(i,i)
            j = 0;
            for j in range(N2):
                idx = j+N1
                if polys2[j].get_sys() == 0:
                    g.add_node(idx)
                    g.add_edge(idx,idx)
                    if len(list(set(polys2[j].get_pos()) & set( polys1[i].get_pos())) ) > 0:
                        g.add_edge(i,idx)
                        g.add_edge(idx,i)
    c = connected_components(g)
    sizes = sorted(c, key = len,reverse = True)
    for p in sizes:
        if a in p:
            return p
    return [] 


#Moves a cluster of polymers
# returns a rejection if any bonds are formed when moving 
# otherwise 
cdef move_cluster(state, r_c, int[:,:,:,:] lattice, int[:,:] tetherLattice):
    # Translate cluster by 1 lattices square
    cdef int i,n_c,j,posY,posZ;
    polys1,polys2= state
    n_c,i = len(r_c),0
    cdef int n1 = len(state[0])
    cdef int n2 = len(state[1])
    l_clust_i,l_clust_f = np.zeros((D1,L,L)),np.zeros((D1,L,L))
    m = nn[crandint(0,5)]
    init_pos,move = [[],[]],[[],[]]
    move[0] = [p1.get_pos() for p1 in polys1]
    move[1] = [p2.get_pos() for p2 in polys2]

    for i in range(n_c):
        c = list(r_c)[i]
        if c >=  n1  and c < n1 + n2:
            pos = state[1][c-n1].get_pos()
            x,y,z = zip(*pos)
            l_clust_i[x,y,z] = 1            
            move[1][c - n1] = [(p[0] + m[0],(p[1]+m[1])%L,(p[2]+m[2])%L) for p in pos]
            if(check_x_bounds(move[1][c-len(polys1)])):
                return state
            x,y,z = zip(*move[1][c-len(polys1)])
            l_clust_f[x,y,z] = 1
        elif c < n1:
            pos = move[0][c]
            x,y,z = zip(*pos)
            l_clust_i[x,y,z] = 1
            move[0][c] = [(p[0] + m[0],(p[1]+m[1])%L,(p[2]+m[2])%L) for p in pos]   
            if(check_x_bounds(move[0][c])):
                return state
            x,y,z = zip(*move[0][c])
            l_clust_f[x,y,z] = 1
        else: 
            return state
            #if m[0] != 0:
            #    return state
            #posY = tetherActiveArr[c - n1 - n2,0]
            #posZ = tetherActiveArr[c - n1 - n2,1]            
            #for k in range(Lt1):
            #    l_clust_i[k,posY,posZ] = 1
            #move_tether[c - n1 - n2] = [(posY + m[1])%L,(posZ + m[2])%L]
            #for k in range(Lt1):
            #    l_clust_f[k,(posY+m[1])%L,(posZ+m[2])%L] = 1
            #return state
    #diff = points in moved cluster that were not in original cluster
    diff = zip(*np.where((l_clust_f - l_clust_i) == 1))
    i = 0
    for i in range(len(diff)):
        x,y,z = diff[i]
        if lattice[0,x,y,z] - l_clust_i[x,y,z] == 1:
            return state
        if lattice[1,x,y,z] - l_clust_i[x,y,z] == 1:
            return state
        if x < Lt1:
            if tetherLattice[y,z] != 0:
                return state
    for i in range(2):
        j = 0;
        n_c = len(move[i])
        for j in range(n_c):
            state[i][j].set_pos(move[i][j])
    return state

cdef check_x_bounds(pos):
    cdef int i = 0;
    cdef int N = len(pos)
    for i in range(N):
        x = pos[i][0]
        if x > D1-1:
            return True
        if x < 0 :
            return True
    return False


#updates lattice positions by moving from final to initial
@cython.boundscheck(False)
@cython.wraparound(False)
cdef int[:,:,:,:] update_lattice(final, initial,int[:,:,:,:] lattice,int ptype):
    cdef int xf,yf,zf,xi,yi,zi;
    xi = initial[0]
    yi = initial[1]
    zi = initial[2]
    lattice[ptype,xi,yi,zi] = 0
    xf = final[0]
    yf = final[1]
    zf = final[2]
    lattice[ptype,xf,yf,zf] = 1
    return lattice


#updates lattice positions by moving from final to initial
@cython.boundscheck(False)
@cython.wraparound(False)
cdef int[:,:,:,:] update_lattice_non_local(final, initial,int[:,:,:,:] lattice,int ptype):
    cdef int xf,yf,zf,xi,yi,zi;
    cdef int k =0
    cdef int N = len(final)
    #print('Before',np.sum(lattice[ptype]))
    for k in range(N):
        xi = initial[k][0]
        yi = initial[k][1]
        zi = initial[k][2]
        lattice[ptype,xi,yi,zi] = 0
    #print('After removal',np.sum(lattice[ptype]))
    for k in range(N):
        xf = final[k][0]
        yf = final[k][1]
        zf = final[k][2]
        lattice[ptype,xf,yf,zf] = 1
    #print('After addition',np.sum(lattice[ptype]))
    return lattice 


## Makes state of polys into a 3D lattice with axis for each poly
#    type. 1/0 indicating occupide/empty
@cython.boundscheck(False)
@cython.wraparound(False)
cdef int[:,:,:,:] make_lattice(state, int sys):
    cdef int i = 0;
    cdef int j = 0;
    cdef int N = 0;
    cdef int k = 0;
    cdef int x,y,z
    cdef int[:,:,:,:] LAT = np.zeros((2,D1,L,L),dtype = np.intc)
    nList = [len(state[0]),len(state[1])]
    N = len(state[0])
    for i in range(N):
        p = state[0][i].get_pos()
        for j in range(Lbulk):
            x = p[j][0]
            y = p[j][1]
            z = p[j][2]
            LAT[0,x,y,z] = 1
    N = len(state[1])
    for i in range(N):
        p = state[1][i].get_pos()
        for j in range(Lbulk):
            x = p[j][0]
            y = p[j][1]
            z = p[j][2]
            LAT[1,x,y,z] = 1
    return LAT

## Remove/add poly to lattice
@cython.boundscheck(False)
@cython.wraparound(False)
cdef int[:,:,:,:] update_lattice_exchange(pos, int pt,int add_rem,int[:,:,:,:] LAT):
    cdef int j,x,y,z;
    for j in range(Lbulk):
        x = pos[j][0]
        y = pos[j][1]
        z = pos[j][2]
        if add_rem == -1:
            LAT[pt,x,y,z] = 0
        elif add_rem == 1:
            LAT[pt,x,y,z] = 1
    return LAT
            

##################################
def print_state(bulk_state,tether,fh,iter):
    print(iter,len(bulk_state[0]),len(bulk_state[1]))
    for i in range(len(bulk_state)):
        polys = bulk_state[i]
        for j in range(len(polys)):
            to_print = ""
            for x,y,z in polys[j].get_pos():
                to_print += '\t'.join([str(x),str(y),str(z)]) + "\t"
            fh.write("%d\t%d\t%d\t%s\n" % (iter,i,j,to_print))
    for j in range(len(tether)):
        y,z = tether[j].get_pos()
        to_print = "\t".join([str(y),str(z)])+"\t"
        fh.write("%d\t%d\t%s\n" % (iter,j,to_print))
    return
                
##################################
#metropolis probabilities for bulk moves
cdef double[:,:] make_probs(float J,float J_nn):
    cdef double[:,:] p = np.zeros((30,70),dtype = np.double)
    cdef int j = 0;
    cdef int k = 0;
    for j in np.arange(-15,15,1):
        for k in np.arange(-35,35,1):
            e = -(J*j + J_nn*k)
            if e <= 0:
                p[j+15,k+35] = 1
            else:
                p[j+15,k+35] = np.exp(-e,dtype = np.double)
    return p

#####################################
cdef double[:] make_cluster_probs():
    cdef double p_not = 0.01;
    cdef double [:] pclus = np.zeros(1000,dtype=np.double)
    for i in np.arange(1,1000,1):
        pclus[i-1] = (1 / i)*p_not
    return pclus

cdef double[:] generate_rands(int randCount):
    cdef double[:] rand_array = np.zeros(randCount,dtype=np.double)
    for i in range(randCount):
        rand_array[i] = crandom()
    return rand_array

############################################
@cython.cdivision(True)
cdef inline  double crandom() except -1:
    return <double>rand()/<double>RAND_MAX
@cython.cdivision(True)
cdef inline int crandint(int lower, int upper) except -1:
    return (rand() % (upper - lower + 1)) + lower

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void fisher_yates_shuffle(int[:] arr, int n):
    cdef int i, j, temp
    for i in range(n-1, 0, -1):
        j = crandint(0, i)
        temp = arr[i]
        arr[i] = arr[j]
        arr[j] = temp
