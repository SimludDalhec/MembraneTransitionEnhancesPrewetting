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
from scipy.signal import find_peaks


##########################################
# SCRIPT TO PLOT ADSORPTION ISOTHERMS FOR SYSTEMS WITH PHASE COEXISTNCES
# APPROACH
#      AVERAGE OVER X CONFIGURRATIONS
#      SCAN A GxG GRID OVER THE AVERAGE
#      CALCULATE DISTRIBTUION OF DENSITIES FOUND IN THE GRID
#      REPEAT FOR ALL N / X SAMPLES, AND POOL DATA 
##########################################

######### PARAMS
D,L = 30,64
grid_size = 4
av_int = 5
Ltether=5
Plotting = False
Jb = 0
M = 0
Ti = 0
########################################################
### MAIN: READSLIST OF FILES, AND STORES DATA IN MASTER ARRAYS
### plots adsoprtoin vs Jbulk, mu_bulk
########################################################
def main():
    fpoly = sys.argv[1]
    fn = re.match('(.*)_pos\.txt',fpoly).group(1)

    Jb,Mu,cT,idx,speedup = re.match('Jb(.*)_u(.*)_cT(.*)_idx(.*)_speedup(.*)_pos\.txt',fpoly).groups(1)
    hPhi = float(Jb)

    p1pos,tether_arr = read_pos_file2(fpoly)

    fn_base = re.match('(.*)\.txt',fpoly).group()
    fData = fn+'_ads_full_data.txt'        
    fh_data = open(fData,'w')
    get_adsorption(p1pos,tether_arr,fn_base,fh_data)
    fh_data.close()
    return


################################
## Reads ising data and stores in nxLxL lattice
################################

#################################
# Reads polymer data and stores in Nxp1xLx3 array
#################################
def read_pos_file(fn):
    p1pos,tether_arr= [],[]
    step_prev = 0
    printing_iter = 0
    for line in open(fn,'r'):
        terms = line.strip().split()
        if int(terms[0]) > printing_iter:
            printing_iter = int(terms[0])
            break
    if printing_iter == 0:
        return -1,-1
    for line in open(fn,'r'):
        terms = line.strip().split()
        if len(terms) == 4:
            try:
                [int(terms[i]) for i in range(len(terms))]
            except:
                return p1pos,tether_arr
            iter,pnum,y,z = [int(terms[i]) for i in range(len(terms))]
            step = iter / printing_iter
            while step >= len(tether_arr):
                tether_arr.append([])
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
            if ptype < 2:
                p1pos[step][ptype].append(poly)
            if step > step_prev:
                step_prev = step
    return p1pos[0:len(p1pos)-1],tether_arr[0:len(tether_arr)-1]


#################################
# Reads polymer data and stores in Nxp1xLx3 array
#################################
def read_pos_file2(fn):
    p1pos,tether_arr= [[[],[]]],[[]]
    step_prev = 0
    printing_iter = 0
    for line in open(fn,'r'):
        terms = line.strip().split()
        if int(terms[0]) > printing_iter:
            printing_iter = int(terms[0])
            break
    if printing_iter == 0:
        return -1,-1
    count = 0
    for line in open(fn,'r'):

        terms = line.strip().split()
        if len(terms) == 4:
            try:
                [int(terms[i]) for i in range(len(terms))]
            except:
                return p1pos,tether_arr            
            iter,pnum,y,z = [int(terms[i]) for i in range(len(terms))]
            step = iter / printing_iter
            if step != step_prev:
                p1pos.append([[],[]])
                tether_arr.append([])
                step_prev = step
                count += 1
            tether_arr[count].append([y,z])
        elif len(terms) > 4:
            iter,ptype,pnum = int(terms[0]),int(terms[1]),int(terms[2])
            step = iter / printing_iter
            if step != step_prev:
                p1pos.append([[],[]])
                tether_arr.append([])
                step_prev = step
                count += 1            
            poly = []
            for i in range(3,len(terms)-2,3):
                x,y,z = terms[i:i+3]
                poly.append((int(x),int(y),int(z)))
            p1pos[count][ptype].append(poly)
    return p1pos,tether_arr


##############################
# stores polymer data onto an L^3 lattice
##############################
def make_lattice(d,tethers):
    latticeb1 = np.zeros((D,L,L))
    latticeb2 = np.zeros((D,L,L))
    tether = np.zeros((L,L))
    d1,d2 = d
    for i in range(len(d1)):
        p = d1[i]
        for j in range(len(p)):
            x,y,z = p[j]
            latticeb1[x,y,z] += 1
    for i in range(len(d2)):
        p = d2[i]
        for j in range(len(p)):
            x,y,z = p[j]
            latticeb2[x,y,z] += 1
    for i in range(len(tethers)):
        y,z = tethers[i]
        if y < L and z < L:
            tether[y,z] += 1
    return latticeb1,latticeb2,tether

#####################################
#### makes lattice with distance #####
#####################################
def DistanceLattice():
    latticeDist = np.zeros((L,L))
    for x in range(L):
        xd = min((L/2-x),abs(L-L/2+x),abs(L-x+L/2))
        for y in range(L):
            yd = min((L/2-y),abs(L-L/2+y),abs(L-y+L/2))
            dist_tot = np.sqrt(xd**2.0 + yd**2.0)
            latticeDist[x,y] = dist_tot
    return latticeDist


##################################
### Spatially average in grid
###################################
def spatial_average(lattice):
    lat_grid_av = np.zeros((L,L))
    for k in range(0,L):
        kidxs = np.arange(k-grid_size/2,k+grid_size/2,1)
        if k-grid_size/2 < 0:
            kidxs = np.concatenate((np.arange(0,k,1),np.arange(k,k+grid_size/2,1),np.arange(L-1,L+k-grid_size/2,-1)))
        elif k+grid_size/2 > L:
            kidxs = np.concatenate((np.arange(k-grid_size/2,k,1),np.arange(k,L,1),np.arange(0,L-(k+grid_size/2),1)))
        for j in range(0,L):
            
            jidxs = np.arange(j-grid_size/2,j+grid_size/2,1)
            if j-grid_size/2 < 0:
                jidxs = np.concatenate((np.arange(0,j,1),np.arange(j,j+grid_size/2,1),np.arange(L-1,L+j-grid_size/2,-1)))
            elif j+grid_size/2 > L:
                jidxs = np.concatenate((np.arange(j-grid_size/2,j,1),np.arange(j,L,1),np.arange(0,L-(j+grid_size/2),1)))

            KS,JS = np.meshgrid(kidxs,jidxs)
            lat_grid_av[k,j] = np.mean(lattice[KS,JS])
    return lat_grid_av


#################
def get_adsorption(poly_data,tethers,fn,fh_data):
    av_LAT1,av_LAT2,av_LAT_tether = calc_averaged_lattice(poly_data,tethers)
    dense_bulk_arr  = np.zeros((len(av_LAT1),2))
    dCount,gCount = 0,0
    dense_drop_arr = np.zeros((len(av_LAT1),4))
    dilute_drop_arr = np.zeros((len(av_LAT1),4))
    dense_av_arr = np.zeros((len(av_LAT1),4))
    area_fraction_arr = np.zeros((len(av_LAT1),2))
    thresh = 1.0
    dense_bulk1 = np.mean(av_LAT1[:,D/2:D,:,:])
    dense_bulk2 = np.mean(av_LAT2[:,D/2:D,:,:])
    lat_2d_shift_thresh_all = np.zeros((len(av_LAT1),L,L))
    lat_3d_shift = np.zeros((len(av_LAT1),D,L,L,3))
    lat_3d_dense_av = np.zeros((len(av_LAT1),D,9))

    ### Finds the phase, where the average density >= thresh ####
    for i in range(len(av_LAT1)):
        ### LATTICE OF TOTAL DENSITY TO CALCULATE WHERE THE PHASE IS
        lat_2d = np.mean(av_LAT2[i,0:5,:,:],axis=0) + np.mean(av_LAT1[i,0:5,:,:],axis=0) + np.mean(av_LAT_tether[i,0:5,:,:],axis=0)
        #### density in the bulk at iteration i #####
        dense_bulk_arr[i,0] = np.mean(np.mean(np.mean(av_LAT1[i,D/2:D,:,:],axis=1),axis=1))
        dense_bulk_arr[i,1] = np.mean(np.mean(np.mean(av_LAT2[i,D/2:D,:,:],axis=1),axis=1))
        ## Average over a grid ###
        lat_2d_grid = spatial_average(lat_2d)
        ## find c.o.m, and shift to center ##
        loc_cm = np.where(lat_2d_grid == np.max(lat_2d_grid))
        for k in range(0,L):
            for j in range(0,L):
                lat_2d_shift_thresh_all[i,k,j] = int(np.mean(lat_2d_grid[k,j] + (1-thresh)))
                for z in range(D):
                    lat_3d_shift[i,z,k,j,0] = av_LAT1[i,z,k,j]
                    lat_3d_shift[i,z,k,j,1] = av_LAT2[i,z,k,j]
                    lat_3d_shift[i,z,k,j,2] = av_LAT_tether[i,z,k,j]
                    
        pos_dense = np.where(lat_2d_shift_thresh_all[i] >= 1) 
        pos_dilute = np.where(lat_2d_shift_thresh_all[i] < 1)
        bulk_dense1 = np.mean(lat_3d_shift[i,D/2:D,pos_dense[0],pos_dense[1],0])
        bulk_dense2 = np.mean(lat_3d_shift[i,D/2:D,pos_dense[0],pos_dense[1],1])
        bulk_dilute1 = np.mean(lat_3d_shift[i,D/2:D,pos_dilute[0],pos_dilute[1],0])
        bulk_dilute2 = np.mean(lat_3d_shift[i,D/2:D,pos_dilute[0],pos_dilute[1],1]) 
        for z in range(D):
            lat_3d_dense_av[i,z,0] = np.mean(lat_3d_shift[i,z,pos_dense[0],pos_dense[1],0]) #- bulk_dense1
            lat_3d_dense_av[i,z,1] = np.mean(lat_3d_shift[i,z,pos_dense[0],pos_dense[1],1]) #- bulk_dense2
            lat_3d_dense_av[i,z,2] = np.mean(lat_3d_shift[i,z,pos_dilute[0],pos_dilute[1],0])# - bulk_dilute1
            lat_3d_dense_av[i,z,3] = np.mean(lat_3d_shift[i,z,pos_dilute[0],pos_dilute[1],1]) #- bulk_dilute2
            lat_3d_dense_av[i,z,4] = np.mean(lat_3d_shift[i,z,:,:,0]) #- bulk_dense1 - bulk_dilute1
            lat_3d_dense_av[i,z,5] = np.mean(lat_3d_shift[i,z,:,:,1]) #- bulk_dense2 - bulk_dilute2
            lat_3d_dense_av[i,z,6] = np.mean(lat_3d_shift[i,z,pos_dense[0],pos_dense[1],2]) 
            lat_3d_dense_av[i,z,7] = np.mean(lat_3d_shift[i,z,pos_dilute[0],pos_dilute[1],2])
            lat_3d_dense_av[i,z,8] = np.mean(lat_3d_shift[i,z,:,:,2])

        if len(pos_dense[0]) > 1:
            dCount += 1
            dense_drop_arr[i,0] = np.sum(lat_3d_dense_av[i,:,0])
            dilute_drop_arr[i,0] = np.sum(lat_3d_dense_av[i,:,2])
            dense_drop_arr[i,1] = np.sum(lat_3d_dense_av[i,:,1])
            dilute_drop_arr[i,1] = np.sum(lat_3d_dense_av[i,:,3])
            dense_drop_arr[i,2] = np.sum(lat_3d_dense_av[i,:,6])/float(Ltether)
            dilute_drop_arr[i,2] = np.sum(lat_3d_dense_av[i,:,7])/float(Ltether)

            dense_av_arr[i,0] = np.sum(lat_3d_dense_av[i,:,4])            
            dense_av_arr[i,1] = np.sum(lat_3d_dense_av[i,:,5])
            dense_av_arr[i,2] = np.sum(lat_3d_dense_av[i,:,8])
            
            area_fraction_arr[i,0] = len(pos_dense[0])/L**2.0
            area_fraction_arr[i,1] = len(pos_dense[1])/L**2.0
        else:
            gCount += 1
            dense_drop_arr[i,0]  = np.NaN
            dense_drop_arr[i,1]  = np.NaN
            dense_drop_arr[i,2] = np.NaN
            dense_drop_arr[i,3] = np.NaN

            dilute_drop_arr[i,0] = np.NaN
            dilute_drop_arr[i,1] = np.NaN
        
            dilute_drop_arr[i,2] = np.NaN
            dilute_drop_arr[i,3] = np.NaN


            dense_av_arr[i,0] = np.sum(lat_3d_dense_av[i,:,4])
            dense_av_arr[i,1] = np.sum(lat_3d_dense_av[i,:,5])
            dense_av_arr[i,2] = np.sum(lat_3d_dense_av[i,:,8])
            area_fraction_arr[i,1] = 0
            area_fraction_arr[i,0] = 0
        fh_data.write('%d\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\n'% (i,dense_drop_arr[i,0],dense_drop_arr[i,1], dense_drop_arr[i,2], dense_drop_arr[i,3],dilute_drop_arr[i,0],dilute_drop_arr[i,1],dilute_drop_arr[i,2],dilute_drop_arr[i,3],dense_av_arr[i,0],dense_av_arr[i,1],dense_av_arr[i,2],area_fraction_arr[i,0],area_fraction_arr[i,1],dense_bulk_arr[i,0],dense_bulk_arr[i,1] ))
        print i,dense_drop_arr[i,:],dilute_drop_arr[i,:],dense_av_arr[i],area_fraction_arr[i,0]
    if Plotting :
        ## 2D Density
        RT=DistanceLattice()
        r  = np.linspace(0.1,np.sqrt(2)*L/2,num=L)
        avArr = np.mean(lat_2d_shift_thresh_all,axis=0)
        avB1Arr = np.mean(lat_3d_shift[:,0,:,:,0],axis=0)
        avB2Arr = np.mean(lat_3d_shift[:,0,:,:,1],axis=0)
        avTArr = np.mean(lat_3d_shift[:,0,:,:,2],axis=0)

        fT = lambda r : avArr[(RT >= r-.5) & (RT < r+.5)].mean()
        f1 = lambda r : avB1Arr[(RT >= r-.5) & (RT < r+.5)].mean()
        f2 = lambda r : avB2Arr[(RT >= r-.5) & (RT < r+.5)].mean()
        fS = lambda r : avTArr[(RT >= r-.5) & (RT < r+.5)].mean()
        meanT = np.vectorize(fT)(r)
        mean1 = np.vectorize(f1)(r)
        mean2 = np.vectorize(f2)(r)
        meanS = np.vectorize(fS)(r)
        fig,ax = plt.subplots(2)
        ax[0].plot(r,meanT,color='purple')
        ax[0].scatter(r,meanT,60,edgecolor='k',color='purple')
        ax[0].plot(r,mean1,color='red')
        ax[0].scatter(r,mean1,60,edgecolor='k',color='red')
        ax[0].plot(r,mean2,color='blue')
        ax[0].scatter(r,mean2,60,edgecolor='k',color='blue')
        ax[0].plot(r,meanS,color='yellow')
        ax[0].scatter(r,meanS,60,edgecolor='k',color='yellow')
        ax[0].set_xlabel('Distance')
        ax[0].set_ylabel('Density')
        ts = np.arange(0,len(dense_drop_arr),1)
        ax[1].plot(ts,dense_drop_arr[:,0],c='r')
        ax[1].scatter(ts,dense_drop_arr[:,0],70,c='r',edgecolor='k')
        ax[1].plot(ts,dense_drop_arr[:,1],c='b')
        ax[1].scatter(ts,dense_drop_arr[:,1],70,c='b',edgecolor='k')
        plt.savefig(fn+'_test_rad_av.svg',transparent=True)

        surf_dense1 = np.mean(np.mean(lat_3d_shift[:,:,:,:,0],axis=0),axis=0)
        surf_dense2 = np.mean(np.mean(lat_3d_shift[:,:,:,:,1],axis=0),axis=0)
        surf_denseT = np.mean(np.mean(lat_3d_shift[:,:,:,:,2],axis=0),axis=0)
        surf_dense1yz = np.mean(np.mean(lat_3d_shift[:,:,:,:,0],axis=0),axis=1)
        surf_dense2yz = np.mean(np.mean(lat_3d_shift[:,:,:,:,1],axis=0),axis=1)
        surf_denseTyz = np.mean(np.mean(lat_3d_shift[:,:,:,:,2],axis=0),axis=1)
        
        fig,ax = plt.subplots(3,2)
        ax[0,0].imshow(surf_dense1,cmap='Reds')
        ax[1,0].imshow(surf_dense2,cmap='Blues')
        ax[2,0].imshow(surf_denseT,cmap='Purples')
        ax[0,1].imshow(surf_dense1yz,cmap='Reds')
        ax[1,1].imshow(surf_dense2yz,cmap='Blues')
        ax[2,1].imshow(surf_denseTyz,cmap='Purples')
        for i in range(len(ax)):
            for j in range(len(ax[i])):
                ax[i,j].set_xticks([])
                ax[i,j].set_yticks([])
        plt.savefig(fn+'2d_plots.svg',transparent=True)
        
        ### 3D Density
        ### Plotting Sum and Difference
        TotalDense = [np.nanmean(lat_3d_dense_av[:,k,0],axis=0) + np.nanmean(lat_3d_dense_av[:,k,1],axis=0) for k in range(D)]
        DiffDense = [np.nanmean(lat_3d_dense_av[:,k,0],axis=0) - np.nanmean(lat_3d_dense_av[:,k,1],axis=0) for k in range(D)]
        plt.figure()
        plt.plot(TotalDense,c='b', label ='Dense Total')
        plt.plot(DiffDense,c='r', label ='Dense Total')
        plt.scatter(np.arange(0,D,1),TotalDense,80,edgecolor='k',c='b')
        plt.scatter(np.arange(0,D,1),DiffDense,80,edgecolor='k',c='r')
        plt.xlabel('Distance From Membrane',fontsize= 30)
        plt.ylabel('Polymer Density',fontsize= 30)
        plt.legend()
        plt.tight_layout()
        plt.savefig(fn + '_dense_prof.svg',transparent=True)
 
    return 
def round_nearest(x, a):
    return round(x / a) * a


## Temporally average the lattice, and shift
def calc_averaged_lattice(poly_data,tethers):
    rng = range(len(poly_data)/4,len(poly_data))
    N = (len(poly_data) -len(poly_data)/4)/av_int 
    averaged_latticeb1 = np.zeros((N,D,L,L))
    averaged_latticeb2 = np.zeros((N,D,L,L))
    averaged_tether_lattice = np.zeros((N,D,L,L))
    count_int = 0
    for i in range(len(poly_data)/4,len(poly_data),av_int):
        count = 0
        for c in range(av_int):
            if count_int < len(averaged_latticeb1) and i+c < len(tethers) and i+c < len(poly_data):
                latticeb1,latticeb2,tether_lat = make_lattice(poly_data[i+c],tethers[i+c])
                averaged_latticeb1[count_int,:,:,:] += (latticeb1)/float(av_int)
                averaged_latticeb2[count_int,:,:,:] += (latticeb2)/float(av_int)
                averaged_tether_lattice[count_int,:,:,:] += (tether_lat)/float(av_int)
            count += 1
        count_int += 1
    return averaged_latticeb1,averaged_latticeb2,averaged_tether_lattice


if __name__ == "__main__":
    main()

