#!/usr/bin/python
import numpy as np
import os
import sys
import matplotlib
import re
matplotlib.use('Agg')
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import networkx as nx
from networkx.algorithms.components.connected import connected_components

D,L = 30,64
Res = False
def main():
    flist = sys.argv[1]
    for f in open(flist,'r'):
        fpoly,fising = f.strip().split()
        fn = re.match('(.*)\.txt',fpoly).group(1)
        p1pos,tether_arr = read_pos_file(fpoly)
        data_lattice = read_ising(fising)
        print len(tether_arr),len(tether_arr[0]),len(p1pos)
        plot_2d_slices(p1pos,tether_arr,data_lattice,fn)
        idx_max = min(len(p1pos),len(tether_arr),len(data_lattice))
        for i in range(idx_max-1,idx_max-3,-1):
            iter = i
            plot_data(p1pos[iter],tether_arr[iter],data_lattice[iter],iter,fn)
    return

def read_ising(fn):
    data = []
    for line in open(fn):
        t = line.strip().split()
        spins = t[1:]
        Lat = np.zeros((L,L))
        #if len(t) > L**2:
        for i in range(len(spins)):
            x,y = i/L,i%L
            try:
                float(spins[i])
            except:
                return data
            if x < 64 and y < 64:
                Lat[x][y] = int(float(spins[i]))
        data.append(Lat)
    return data


def read_pos_file(fn):
    p1pos,tether_arr= [],[]
    step_prev = 0
    printing_iter = 0
    for line in open(fn,'r'):
        terms = line.strip().split()
        if int(terms[0]) > printing_iter:
            printing_iter = int(terms[0])
            break
    for line in open(fn,'r'):
        terms = line.strip().split()
        if len(terms) == 4:
            iter,pnum,y,z = [int(terms[i]) for i in range(len(terms))]
            step = iter / printing_iter
            while step >= len(tether_arr):
                tether_arr.append([[]])
            tether_arr[step].append([y,z])
            if step > step_prev:
                step_prev = step
        elif len(terms) > 4:
            iter,ptype,pnum = int(terms[0]),int(terms[1]),int(terms[2])
            step = iter / printing_iter
            poly = []
            while step >= len(p1pos):
                p1pos.append([[],[]])
            for i in range(3,len(terms)-2,3):
                x,y,z = terms[i:i+3]
                poly.append((int(x),int(y),int(z)))
            p1pos[step][ptype].append(poly)
            if step > step_prev:
                step_prev = step
            #print len(p1pos[len(p1pos)-1][0]), len(p1pos[len(p1pos)-1][1]), len(tether_arr[len(p1pos)-2]),len(tether_arr[2])
    return p1pos,tether_arr
def get_spatial_density(lattice_i,lattice_ti):
    grid_size=6
    d_avP = np.zeros((D,L,L))
    d_avT = np.zeros((D,L,L))
    for k in range(0,D-grid_size):    
        for j in range(0,L-grid_size):
            for l in range(0,L-grid_size):
                d_avP[k,j,l] = np.mean(lattice_i[k:(k+grid_size),j:(j+grid_size)%L,l:(l+grid_size)%L])
                #d_avT[k,j,l] = np.mean(lattice_ti[k:(k+grid_size),j:(j+grid_size)%L,l:(l+grid_size)%L])
    return d_avP#,d_avT


def plot_data(data_poly1,data_tether,data_lattice,num,fn):
    fig = plt.figure(dpi=200,figsize=(8,6))
    ax = fig.add_subplot(111, projection='3d')
    polys1,polys2 = data_poly1
    tethers1 = data_tether
    angle = 45
    #Setting view stuff
    ax.view_init(30,angle)
    ax.set_xlim(0,L-1)
    ax.set_ylim(0,L-1)
    ax.set_zlim(0,D-1)
    fname = fn + "{:05d}.png".format(num)
    ZZ = np.zeros((L,L))   
    XX, YY = np.meshgrid(np.arange(0,L,1), np.arange(0,L,1))
    RELOC = True
    lattice_i,tethers_i = make_lattice(data_poly1,data_tether)
    spatial_av = get_spatial_density(lattice_i,tethers_i)
    loc_max_tethers = np.where(spatial_av == np.max(spatial_av))
    loc_maxX = loc_max_tethers[0][0]
    loc_maxY = loc_max_tethers[1][0]
    loc_maxZ = loc_max_tethers[2][0]
    print loc_max_tethers
    print loc_maxX,loc_maxY,loc_maxZ


    data_latticeP = np.zeros((L,L))
    if RELOC == True:
        for i in range(L):
            yp = (i - loc_maxY + L/2)%L
            for j in range(L):
                zp = (j - loc_maxZ + L/2)%L
                #data_latticeP[yp,zp] = data_lattice[i,j]
                data_latticeP[yp,zp] = data_lattice[i,j]
    ZZ = np.zeros((L,L))   
    XX, YY = np.meshgrid(np.arange(0,L,1), np.arange(0,L,1))

    #ax.plot_surface(XX,YY, ZZ, alpha = 0.5,linewidth=0,facecolors=plt.cm.Greys(np.transpose(data_latticeP)))
    ax.contourf(XX,YY,np.transpose(data_latticeP),zdir='z',offset = 0,alpha = 0.5,cmap=cm.Greys,corner_mask=False)

                
    #plotting
    for i in range(len(polys1)):
        lines_connect = []
        for x1,y1,z1 in polys1[i]:
            if len(lines_connect) == 0:
                lines_connect.append([])
            curr_line = lines_connect[len(lines_connect)-1]
            if len(curr_line) > 0:
                last_point = curr_line[len(curr_line)-1]
                if (last_point[0] == L-1 and y1 ==0) or (last_point[0] == 0 and y1 == L-1): 
                    lines_connect.append([])
                elif (last_point[1] == L-1 and z1 ==0) or (last_point[1] == 0 and z1 == L-1): 
                    lines_connect.append([])
            if RELOC:
                xP = x1
                yP = (y1 - loc_maxY + L/2)%L 
                zP = (z1 - loc_maxZ + L/2)%L 
                lines_connect[len(lines_connect)-1].append((yP,zP,xP))
            else:
                lines_connect[len(lines_connect)-1].append((y1,z1,x1))
        for lines in lines_connect:
            for j in range(1,len(lines)):
                #ax.scatter([lines[j-1][0],lines[j][0]],[lines[j-1][1],lines[j][1]],[lines[j-1][2],lines[j][2]],c='steelblue',alpha=0.7,edgecolor='k')
                ax.scatter([lines[j-1][0],lines[j][0]],[lines[j-1][1],lines[j][1]],[lines[j-1][2],lines[j][2]],c='indianred',alpha=0.7,edgecolor='k')

    for i in range(len(polys2)):
        lines_connect = []
        for x2,y2,z2 in polys2[i]:
            if len(lines_connect) == 0:
                lines_connect.append([])
            curr_line = lines_connect[len(lines_connect)-1]
            if len(curr_line) > 0:
                last_point = curr_line[len(curr_line)-1]
                if (last_point[0] == L-1 and y2 ==0) or (last_point[0] == 0 and y2 == L-1): 
                    lines_connect.append([])
                elif (last_point[1] == L-1 and z2 ==0) or (last_point[1] == 0 and z2 == L-1): 
                    lines_connect.append([])
            if RELOC:
                xP = x2
                yP = (y2 - loc_maxY +L/2)%L 
                zP = (z2 - loc_maxZ +L/2)%L
                lines_connect[len(lines_connect)-1].append((yP,zP,xP))
            else:
                lines_connect[len(lines_connect)-1].append((y2,z2,x2))
        for lines in lines_connect:
            for j in range(1,len(lines)):
                #ax.scatter([lines[j-1][0],lines[j][0]],[lines[j-1][1],lines[j][1]],[lines[j-1][2],lines[j][2]],c='indianred',alpha=0.7,edgecolor='k')
                ax.scatter([lines[j-1][0],lines[j][0]],[lines[j-1][1],lines[j][1]],[lines[j-1][2],lines[j][2]],c='lightblue',alpha=0.7,edgecolor='k')


    for i in range(1,len(tethers1)):
        y2,z2 = tethers1[i]
        for l in range(5):
            if RELOC:
                yP = (y2 - loc_maxY + L/2)%L 
                zP = (z2 - loc_maxZ + L/2)%L
                
            #ax.scatter(yP,zP,l,s=20,c="gold",alpha=0.75,edgecolors='k')
            ax.scatter(yP,zP,l,s=20,c="indianred",alpha=0.75,edgecolors='k')
    ax.plot([0,0],[0,L],zs=0,color='k')
    ax.plot([0,L],[L,L],zs=0,color='k')
    ax.plot([L,L],[L,0],zs=0,color='k')
    ax.plot([L,0],[0,0],zs=0,color='k')
    #ax.quiver(0,0,0,0,0,D-3,color='k')
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)
    plt.setp(ax.get_zticklabels(), visible=False)
    ax.grid(False)
    ax.set_yticks([])
    ax.set_zticks([])
    ax.set_xticks([])
    ax.set_axis_off()
    plt.savefig(fname,transparent=True)
    plt.close()
    return


def plot_2d_slices(polys,tethers,ising,fn):
    tether_density = np.zeros((L,L))
    surf_density = np.zeros((L,L))
    for i in range(len(polys)/2):
        lattice,tether= make_lattice(polys[i],tethers[i])
        tether_density += tether
        surf_density +=  sum(lattice[0:5,:,:]) / (len(polys)*L)

    fig,ax = plt.subplots(3,1)
    ax[0].imshow(tether_density,cmap='Purples')#,vmin=0,vmax=1)
    ax[1].imshow(np.mean(ising,axis = 0),cmap='Greys')#,vmin = -1,vmax = 1)
    ax[2].imshow(surf_density,cmap = 'RdPu')#,vmin=0,vmax=1)
    for i in range(len(ax)):
        ax[i].set_xticks([])
        ax[i].set_yticks([])
    plt.savefig(fn + '_surf_density.svg',transparent='True')
    return 


def make_lattice(d,tethers):
    lattice = np.zeros((D,L,L))
    tether = np.zeros((L,L))
    d1,d2 = d
    for i in range(len(d1)):
        p = d1[i]
        for j in range(len(p)):
            x,y,z = p[j]
            if x < D and y < L and z < L :
                lattice[x,y,z] += 1
    for i in range(len(d2)):
        p = d2[i]
        for j in range(len(p)):
            x,y,z = p[j]
            if x < D and y < L and z < L :
                lattice[x,y,z] += 1

    for i in range(1,len(tethers)):
        y,z = tethers[i]
        tether[y,z] += 1 


    return lattice,tether


#makes a graph and returns the connected component that "a" is in 
# a slow but sure way to get connected groups of polymers
def make_graph(polys1,polys2, N1, N2):
    g = nx.Graph()
    i = 0;
    j = 0;
    idx = 0;
    for i in range(N1):
        pPos1 = polys1[i]
        g.add_node(i)
        g.add_edge(i,i)
        j = 0;
        for j in range(N2):
            pPos2 = polys2[j]
            idx = j+N1
            g.add_node(idx)
            g.add_edge(idx,idx)
            if len(list(set(pPos2) & set(pPos1))) > 0:
                    g.add_edge(i,idx)
                    g.add_edge(idx,i)
    c = connected_components(g)
    sizes = sorted(c, key = len,reverse = True)
    return sizes, c

if __name__ == "__main__":
    main()

