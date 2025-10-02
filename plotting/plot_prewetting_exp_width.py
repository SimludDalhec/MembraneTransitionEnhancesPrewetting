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
import scipy.optimize

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
grid_size =4
av_int = 2
Jb = 0
M = 0
Ti = 0
########################################################
### MAIN: READSLIST OF FILES, AND STORES DATA IN MASTER ARRAYS
### plots adsoprtoin vs Jbulk, mu_bulk
########################################################
def main():
    flist = sys.argv[1]
    fits = {}
    for line in open(flist,'r'):
        fpoly = line.strip().split()[0]
        fn = re.match('(.*)\.txt',fpoly).group(1)
        matches = re.match('Ti(.*)_Jb(.*)_u(.*)_u2(.*)_M(.*)_idx(.*)_speedup.*_pos\.txt',fpoly)
        if matches is not None:
            Ti,Jb,Mu,Mu2,M,idx = matches.groups(1)
        hPhi = float(Jb)
        p1pos,tether_arr = read_pos_file2(fpoly)
        fn_base = re.match('(.*)\.txt',fpoly).group()
        XZB1,XZB2,fit_params,bulk_conc,area_fraction,densities = get_adsorption(p1pos,tether_arr,fn_base)
        if float(Mu) not in fits.keys():
            fits[float(Mu)] = [[],[],[],[]]
            print(Mu,idx,area_fraction)
            
        fits[float(Mu)][0].append(bulk_conc)
        fits[float(Mu)][1].append(fit_params)
        fits[float(Mu)][2].append(densities)
        fits[float(Mu)][3].append(area_fraction)
        print('Debug',fits[float(Mu)])
    fWidth = 'prewet_width_Jb1.55_'+flist+'_component.txt'        
    fhw = open(fWidth,'w')
    mus = sorted(fits.keys())
    print(mus)
    widths= []
    for mu in mus:
        if len(fits[mu]) > 0:
            bulk_concs,fit_params,densitiess,areas= fits[mu]
            mean_bulk = np.mean(bulk_concs,axis=0)
            mean_fit = np.mean(fit_params,axis=0)
            widths.append(mean_fit[1])
            print(mu,widths)
            mean_densities = np.mean(densitiess,axis=0)
            mean_area = np.mean(areas,axis=0)
            print(densities,mean_densities)
            fhw.write('%0.7f\t%0.7f\t%0.7f\t%0.7f\t%0.7f\t%0.7f\t%0.7f\t%0.7f\t%0.7f\n' % (round(float(mu),4), mean_bulk,mean_fit[0],mean_fit[1],mean_fit[2],mean_densities[0],mean_densities[1],mean_densities[2],mean_area))            
    fhw.close()
    
    mus = sorted(fits.keys())
    plt.figure()
    plt.plot(mus,widths,color='k',linestyle='--')
    plt.scatter(mus,widths,100,edgecolor='k',color = 'green')
    plt.xlabel('Chemical Potential',fontsize=22)
    plt.ylabel('Thickness',fontsize=22)
    plt.tight_layout()
    plt.savefig('fitted_thickness.svg',transparent=True)
    return


################################
## Reads ising data and stores in nxLxL lattice
################################
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
    return p1pos,tether_arr


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
def get_adsorption(poly_data,tethers,fn):
    av_LAT1,av_LAT2,av_LAT_tether = calc_averaged_lattice(poly_data,tethers)
    dense_bulk_arr  = np.zeros((len(av_LAT1),2))
    dCount,gCount = 0,0
    dense_drop_arr = np.zeros((len(av_LAT1),2))
    dilute_drop_arr = np.zeros((len(av_LAT1),2))
    dense_av_arr = np.zeros((len(av_LAT1),2))
    area_fraction_arr = np.zeros((len(av_LAT1),2))
    thresh = 0.75
    dense_bulk1 = np.mean(av_LAT1[:,D/2:D,:,:])
    dense_bulk2 = np.mean(av_LAT2[:,D/2:D,:,:])
    lat_2d_shift_thresh_all = np.zeros((len(av_LAT1),L,L))
    lat_3d_shift = np.zeros((len(av_LAT1),D,L,L,3))
    lat_3d_dense_av = np.zeros((len(av_LAT1),D,8))

    ### Finds the phase, where the average density >= thresh ####
    for i in range(len(av_LAT1)):
        ### LATTICE OF TOTAL DENSITY TO CALCULATE WHERE THE PHASE IS
        lat_2d = np.mean(av_LAT2[i,0:5,:,:],axis=0) + np.mean(av_LAT1[i,0:5,:,:],axis=0) + np.mean(av_LAT_tether[i,0:5,:,:],axis=0)
        #### density in the bulk at iteration i #####
        dense_bulk_arr[i,0] = np.mean(av_LAT1[i,D/2:D,:,:])
        dense_bulk_arr[i,1] = np.mean(av_LAT2[i,D/2:D,:,:])
        ## Average over a grid ###
        lat_2d_grid = spatial_average(lat_2d)
        ## find c.o.m, and shift to center ##
        loc_cm = np.where(lat_2d_grid == np.max(lat_2d_grid))
        for k in range(0,L):
            for j in range(0,L):
                lat_2d_shift_thresh_all[i,(k-loc_cm[0][0] +L/2)%L ,(j-loc_cm[1][0]+L/2)%L] = int(np.mean(lat_2d_grid[k,j] + (1-thresh)))
                for z in range(D):
                    lat_3d_shift[i,z,(k-loc_cm[0][0] +L/2)%L ,(j-loc_cm[1][0]+L/2)%L,0] = av_LAT1[i,z,k,j]
                    lat_3d_shift[i,z,(k-loc_cm[0][0] +L/2)%L ,(j-loc_cm[1][0]+L/2)%L,1] = av_LAT2[i,z,k,j]
                    lat_3d_shift[i,z,(k-loc_cm[0][0] +L/2)%L ,(j-loc_cm[1][0]+L/2)%L,2] = av_LAT_tether[i,z,k,j]
        pos_dense = np.where(lat_2d_shift_thresh_all[i] >= 1) 
        pos_dilute = np.where(lat_2d_shift_thresh_all[i] < 1)
        for z in range(D):
            lat_3d_dense_av[i,z,0] = np.mean(lat_3d_shift[i,z,pos_dense[0],pos_dense[1],0]) - dense_bulk_arr[i,0]
            lat_3d_dense_av[i,z,1] = np.mean(lat_3d_shift[i,z,pos_dense[0],pos_dense[1],1]) - dense_bulk_arr[i,1]
            lat_3d_dense_av[i,z,2] = np.mean(lat_3d_shift[i,z,pos_dense[0],pos_dense[1],2])
            lat_3d_dense_av[i,z,3] = np.mean(lat_3d_shift[i,z,pos_dilute[0],pos_dilute[1],0]) - dense_bulk_arr[i,0]
            lat_3d_dense_av[i,z,4] = np.mean(lat_3d_shift[i,z,pos_dilute[0],pos_dilute[1],1]) - dense_bulk_arr[i,1]
            lat_3d_dense_av[i,z,5] = np.mean(lat_3d_shift[i,z,pos_dilute[0],pos_dilute[1],2]) - dense_bulk_arr[i,1]
            lat_3d_dense_av[i,z,6] = np.mean(lat_3d_shift[i,z,:,:,0]) - dense_bulk_arr[i,1]
            lat_3d_dense_av[i,z,7] = np.mean(lat_3d_shift[i,z,:,:,1]) - dense_bulk_arr[i,1]

        if len(pos_dense[0]) > 1:
            dCount += 1
            dense_drop_arr[i,0] = np.sum(lat_3d_dense_av[i,:,0])
            dilute_drop_arr[i,0] = np.sum(lat_3d_dense_av[i,:,2])
            dense_drop_arr[i,1] = np.sum(lat_3d_dense_av[i,:,1])
            dilute_drop_arr[i,1] = np.sum(lat_3d_dense_av[i,:,3])
            dense_av_arr[i,0] = np.sum(lat_3d_dense_av[i,:,6])
            dense_av_arr[i,1] = np.sum(lat_3d_dense_av[i,:,7])
            area_fraction_arr[i,0] = len(pos_dense[0])/L**2.0
            area_fraction_arr[i,1] = len(pos_dense[1])/L**2.0
        else:
            gCount += 1
            dense_drop_arr[i,0]  = np.NaN
            dense_drop_arr[i,1]  = np.NaN
            dilute_drop_arr[i,0] = np.NaN
            dilute_drop_arr[i,1] = np.NaN
            dense_av_arr[i,0] = np.sum(lat_3d_dense_av[i,:,6])
            dense_av_arr[i,1] = np.sum(lat_3d_dense_av[i,:,7])
            area_fraction_arr[i,1] = np.NaN
            area_fraction_arr[i,0] = np.NaN
        width1 = [lat_3d_dense_av[i,z,0] -  lat_3d_dense_av[i,z,2] for z in range(D)]    
        width2 = [lat_3d_dense_av[i,z,1] -  lat_3d_dense_av[i,z,3] for z in range(D)]

    ### 3D Density 
    bulk1Densep = np.nanmean(lat_3d_dense_av[:,:,0],axis=0)
    bulk2Densep = np.nanmean(lat_3d_dense_av[:,:,1],axis=0)
    tetherDensep = np.nanmean(lat_3d_dense_av[:,:,2],axis=0)
    bulk1Dilutep = np.nanmean(lat_3d_dense_av[:,:,3],axis=0)
    bulk2Dilutep = np.nanmean(lat_3d_dense_av[:,:,4],axis=0)
    tetherDilutep = np.nanmean(lat_3d_dense_av[:,:,5],axis=0)
    

    xBulk = np.arange(4,D,1)
    xBulkF = np.arange(0,D-4,1)
    denseBulk = bulk2Densep[xBulk]
    bulk_conc = (dense_bulk1+dense_bulk2)/2    
    phi_0 = denseBulk[0]
    avg_bulk1=np.mean(bulk1Densep[0:5])
    avg_bulk2=np.mean(bulk2Densep[0:5])
    avg_tether=np.mean(tetherDensep[0:5])
    densities = [avg_bulk1,avg_bulk2,avg_tether]
    
    params = None
    fit = None
    area_fraction = np.nanmean(area_fraction_arr[:,0])    
    if np.sum(np.isnan(denseBulk)) > 0 or area_fraction < 0.001:
        params,fit = ([0,np.inf,0],0)
        densities = [np.NaN,np.NaN,np.NaN]
        print "NO FIT"
    else:
        params,fit = scipy.optimize.curve_fit(expon, xBulkF, denseBulk, (phi_0, 1, 0))
        print "FIT PARAMS:", params, bulk_conc,area_fraction,densities

    XZB1 = None
    XZB2 = None
    return XZB1,XZB2,params,bulk_conc,area_fraction,densities

def expon(x, m, t, l):
    return m * np.exp(-t*x) + l

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

