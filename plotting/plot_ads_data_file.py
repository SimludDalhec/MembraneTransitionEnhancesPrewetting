#!/usr/bin/python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys,re,math
D,L=30,64
def read_file(f,ftype):
    count = 0
    termcount = 11
    if ftype == 1:
        termcount = 16
    for line in open(f,'r'):
        count+= 1
    data_arr = np.zeros((count,15))
    for line in open(f,'r'):
        terms=line.strip().split()
        termsR = np.zeros(len(terms))
        for i  in range(len(terms)):
            if terms[i] == 'nan':
                termsR[i] = np.NaN
            else:
                termsR[i] = terms[i]
        if len(termsR) == termcount:
            idx = int(terms[0])
            if ftype == 0:
                i,drop_ads1,drop_ad2,dil_ads1,dil_ads2,av_ads1,av_ads2,area_frac,area_frac2,av1,av2 = [float(t) for t in termsR]
            else:
                i,drop_ads1,drop_ad2,drop_ads_tether,drop_ising,dil_ads1,dil_ads2,dil_ads_tether,dil_ising,av_ads1,av_ads2,av_ads_tether,area_frac,area_frac2,av1,av2 = [float(t) for t in termsR]
                data_arr[idx,10] = drop_ads_tether
                data_arr[idx,11] = dil_ads_tether
                data_arr[idx,12] = drop_ising
                data_arr[idx,13] = dil_ising
                data_arr[idx,14] = av_ads_tether

            data_arr[idx,0] = drop_ads1
            data_arr[idx,1] = drop_ad2
            data_arr[idx,2] = dil_ads1
            data_arr[idx,3] = dil_ads2
            data_arr[idx,4] = av_ads1
            data_arr[idx,5] = av_ads2
            data_arr[idx,6] = area_frac
            data_arr[idx,7] = area_frac2
            data_arr[idx,8] = av1
            data_arr[idx,9] = av2
    return data_arr

def plot_ads(data,norm):
    Jbs = sorted(data.keys())
    plt.figure()
    for J in Jbs:
        dataF = 'av_adsorption_data_Ti' + str(J)+'.txt'
        fh = open(dataF,'w')
        Mus = sorted(data[J].keys())
        av_bulkAll = [np.mean(data[J][Mu][0:len(data[J][Mu]),8]) + np.mean(data[J][Mu][0:len(data[J][Mu]),9])for Mu in Mus]
        av_bulk1 = [np.mean(data[J][Mu][0:len(data[J][Mu]),8])for Mu in Mus]
        av_bulk2 = [np.mean(data[J][Mu][0:len(data[J][Mu]),9]) for Mu in Mus]
        av_ads2 = [(np.nanmean(data[J][Mus[k]][0:len(data[J][Mus[k]]),5]))for k in range(len(Mus))]  
        av_ads1 = [(np.nanmean(data[J][Mus[k]][0:len(data[J][Mus[k]]),4])) for k in range(len(Mus))] 
        av_err2 = [(np.nanstd(data[J][Mus[k]][0:len(data[J][Mus[k]]),5])) for k in range(len(Mus))] 
        av_err1 = [(np.nanstd(data[J][Mus[k]][0:len(data[J][Mus[k]]),4])) for k in range(len(Mus))]

        col = np.random.rand(3,)
        label=None
        if J == 2.0:
            label = "Single Component"
            col = 'orange'
        elif J > 0 and J < 10:
            label = "Near-Critical, T = "+str(round(J,3))
            col = 'purple'
        elif J == 0:
            label = "Solid"
            col = 'indianred'
        elif J == 10:
            label = "Single Component"
            col = 'orange'
        elif J == -1:
            label = "Bulk 1 only"
            col = 'darkgreen'
            plt.scatter(av_bulk1,av_ads1,80,edgecolor='k',color=col)
            plt.errorbar(av_bulk1,av_ads1,yerr=av_err1,color=col,label=label)
            fh.write("\t".join([str(k) for k in av_ads1]))        
        if J != -1:
            plt.scatter(av_bulk2,av_ads2,80,edgecolor='k',color=col)
            plt.errorbar(av_bulk2,av_ads2,yerr=av_err2,color=col,label=label)
            fh.write("\t".join([str(k) for k in av_ads2]))
        fh.close()
    plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel(r'$c_{bulk}$',fontsize=32)
    plt.ylabel('Average Adsorption',fontsize=32)
    plt.tight_layout()
    plt.gca().set_aspect(0.75)
    plt.xlim(0.00001,0.025)
    plt.legend()
    plt.savefig('av_ads_plot.svg',transparent=True)
    
    plt.figure()
    for J in Jbs:
        Mus = sorted(data[J].keys())
        Ldata = [len(data[J][Mu]) for Mu in Mus]
        av_bulkAll = [np.mean(data[J][Mu][0:len(data[J][Mu]),8]) + np.mean(data[J][Mu][0:len(data[J][Mu]),9])for Mu in Mus]
        av_bulk1 = [np.mean(data[J][Mu][0:len(data[J][Mu]),8]) for Mu in Mus]
        av_bulk2 = [ np.mean(data[J][Mu][0:len(data[J][Mu]),9]) for Mu in Mus]
        av_ads2 = [(np.nanmean(data[J][Mus[k]][0:len(data[J][Mus[k]]),5]))for k in range(len(Mus))] 
        av_ads1 = [(np.nanmean(data[J][Mus[k]][0:len(data[J][Mus[k]]),4]) ) for k in range(len(Mus))] 
        drop_ads = [(np.nanmean(data[J][Mus[k]][0:len(data[J][Mus[k]]),0])) for k in range(len( Mus))] 
        drop_ads2 = [(np.nanmean(data[J][Mus[k]][0:len(data[J][Mus[k]]),1])) for k in range(len(Mus))]
        drop_err2 = [np.nanstd(data[J][Mus[k]][:,1]) for k in range(len(Mus))]
        drop_err = [np.nanstd(data[J][Mus[k]][:,0]) for k in range(len(Mus))]
        area = [np.nanmean(data[J][Mus[k]][:,6],axis=0) for k in range(len(Mus))]
        print(Mus,av_bulk1)
        area_err = [np.nanstd(data[J][Mus[k]][:,6]) for Mu in Mus]
        col = np.random.rand(3,)
        label=None
        if J == 10:
            label = "Single Component"
            col = 'orange'
        elif J > 0 and J < 10:
            label = "Near-Critical, T = "+str(round(J,3))
            col = 'purple'
        elif J == 10:
            label = "No Tethers"
            col = 'steelblue'
        elif J == 0:
            label = "Solid"
            col = 'indianred'
        elif J == -1:
            label = "Bulk 1 only"
            col = 'darkgreen'
            plt.scatter(av_bulk1,drop_ads,80,edgecolor='k',color=col,label=label)
            plt.errorbar(av_bulk1,drop_ads,yerr=drop_err,color=col)
        if J != -1:
            plt.scatter(av_bulk2,drop_ads2,80,edgecolor='k',color=col,label=label)
            plt.errorbar(av_bulk2,drop_ads2,yerr=drop_err2,color=col)
        dataF = 'adsorption_data_Jb' + str(J)+'.txt'
        fh = open(dataF,'w')
        fh.write('Chemical Potential\t soluble_bulk1\t soluble_bulk2\t adsorption_av_bulk1\t adsorption_ads_bulk2\t adsorption_drop_ads1\t adsorption_drop_bulk2\t area_fraction\t drop_stdev_bulk1\t drop_stdev_bulk2\t area_stdev\n')
        for k in range(len(Mus)):
            fh.write('%0.5f\t%0.9f\t%0.9f\t%0.9f\t%0.9f\t%0.9f\t%0.9f\t%0.9f\t%0.9f\t%0.9f\t%0.9f\n' % (Mus[k], av_bulk1[k], av_bulk2[k], av_ads1[k], av_ads2[k],drop_ads[k],drop_ads2[k], area[k],drop_err[k],drop_err2[k],area_err[k]))
        fh.close()
    plt.legend()
    plt.xscale('log')
    plt.xlabel(r'$c_{bulk}$',fontsize=32)
    plt.ylabel('Adsorption',fontsize=32)
    plt.tight_layout()
    plt.gca().set_aspect(0.75)
    plt.xlim(0.00001,0.025)
    plt.savefig('ads_drop_plot.svg',transparent=True)

def plot_area(data):
    plt.figure()
    Tis = sorted(data.keys())
    for Ti in Tis:
        dataF = 'area_data_Ti' + str(Ti)+'.txt'
        fh = open(dataF,'w')
        Mus = sorted(data[Ti].keys())
        area = [np.mean(data[Ti][Mu][:,6]) for Mu in Mus]
        area_err = [np.nanstd(data[Ti][Mu][:,6]) for Mu in Mus]
        av_bulkAll = [np.mean(data[Ti][Mu][0:len(data[Ti][Mu]),8]) + np.mean(data[Ti][Mu][0:len(data[Ti][Mu]),9])for Mu in Mus]
        av_bulk1 = [np.mean(data[Ti][Mu][0:len(data[Ti][Mu]),8])for Mu in Mus]

        col = np.random.rand(3,)
        if Ti == 10:
            label = "Single Component"
            col = 'orange'
        elif Ti > 0 and Ti < 10:
            label = "Near-Critical, T = "+str(round(Ti,3))
            col = 'purple'
        elif Ti == 0:
            label = "Solid"
            col = 'indianred'
        elif Ti == 2:
            label = "No Tethers"
            col = 'steelblue'
        elif Ti == -1:
            label = "Blue only"
            col = 'darkgreen'
            plt.scatter(av_bulkAll,area,80,edgecolor='k',color=col)
            plt.errorbar(av_bulkAll,area,yerr=area_err,color=col)
            fh.write('\t'.join([str(a) for a in area]))
        if Ti != -1:
            plt.scatter(av_bulkAll,area,80,edgecolor='k',color=col)
            plt.errorbar(av_bulkAll,area,yerr=area_err,color=col)
            fh.write('\t'.join([str(a) for a in area]))
        fh.close()
    plt.legend()
    plt.xscale('log')
    plt.xlabel(r'$c_{bulk}$',fontsize=32)
    plt.ylabel('Area Fraction',fontsize=32)
    plt.gca().set_aspect(0.75)
    plt.ylim(-0.0025,0.09)
    plt.tight_layout()
    plt.xlim(0.00005,0.05)
    plt.savefig('area_plot.svg',transparent=True)


def plot_surface_data(data):
    Tis = sorted(data.keys())
    plt.figure()
    for Ti in Tis:
        if Ti > 0:
            mus = sorted(data[Ti].keys())
            av_bulk1 = [np.mean(data[Ti][Mu][0:len(data[Ti][Mu]),8])for Mu in mus]
            avbulk2 = [np.mean(data[Ti][Mu][0:len(data[Ti][Mu]),9]) for Mu in mus]

            tether_dense = [np.nanmean(data[Ti][mu][:,10]) for mu in mus]
            tether_dil = [np.nanmean(data[Ti][mu][:,11]) for mu in mus]
            tether_av = [np.nanmean(data[Ti][mu][:,14]) for mu in mus]
            ising_dense = [-np.nanmean(data[Ti][mu][:,12]) for mu in mus]
            ising_dil = [-np.nanmean(data[Ti][mu][:,13]) for mu in mus]
            label = 'Near Critical'
            col='purple'
            if Ti == 2.0:
                label = 'Single component'
                col = 'orange'
            plt.scatter(avbulk2,tether_dense,80,color=col,edgecolor='k')
            plt.plot(avbulk2,tether_dense,color=col)
            plt.scatter(avbulk2,tether_dil,80,color=col,edgecolor='w')
            plt.plot(avbulk2,tether_dil,color=col)
            
            plt.scatter(avbulk2,ising_dense,80,color='black',edgecolor='k')
            plt.plot(avbulk2,ising_dense,color='black')
            plt.scatter(avbulk2,ising_dil,80,color='black',edgecolor='w')
            plt.plot(avbulk2,ising_dil,color='black')
            dataF = 'surface_data_Ti' + str(Ti)+'.txt'
            fh = open(dataF,'w')
            for k in mus:
                avbulk1 = np.mean(data[Ti][k][0:len(data[Ti][k]),8])
                b1_dense = np.nanmean(data[Ti][k][:,0])
                b2_dense = np.nanmean(data[Ti][k][:,1])
                b1_dil = np.nanmean(data[Ti][k][:,2])
                b2_dil = np.nanmean(data[Ti][k][:,3])
                b1_av = np.nanmean(data[Ti][k][:,4])
                b2_av = np.nanmean(data[Ti][k][:,5])


                tether_dense = np.nanmean(data[Ti][k][:,10])
                tether_dil = np.nanmean(data[Ti][k][:,11])
                tether_av = np.nanmean(data[Ti][k][:,14])
                ising_dense = np.nanmean(data[Ti][k][:,12])
                ising_dil = np.nanmean(data[Ti][k][:,13])
                fh.write('%0.5f\t%0.f5\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\n' % ( k,avbulk1,b1_dense,b1_dil,b1_av,b2_dense,b2_dil,b2_av,tether_dense,tether_dil,tether_av,ising_dense,ising_dil))
            fh.close()
    plt.ylim(-1,5)
    plt.xscale('log')
    plt.savefig('test_surface_prop.svg',transparent=True)
    return
data_by_mu = {}
for fname in open(sys.argv[1],'r'):
    matches = re.match('Ti(.*)_Jb(.*)_u(.*)_u2(.*)_M(.*)_ads_full_data.txt',fname)
    matches2 = re.match('Ti(.*)_Jb(.*)_u(.*)_M(.*)_pos_ads.txt',fname)
    matches3 = re.match('noTether\/Ti(.*)_Jb(.*)_u(.*)_u2(.*)_M(.*)_pos_ads.txt',fname)
    ftype = 0
    if matches is not None:        
        Ti,Jb,mu,mu2,M = matches.groups()
        ftype = 1
    elif matches2 is not None:
        Ti,Jb,mu,M = matches2.groups()
        mu2 = mu    
        Ti = -1
    elif matches3 is not None:
        Ti,Jb,mu,mu2,M=matches3.groups()
        Ti = 10
    else:
        Jb,mu,cT = re.match('Jb(.*)_u(.*)_cT(.*)_pos_ads.txt',fname).groups()
        Ti = 0
    if float(Ti) not in data_by_mu.keys():
        data_by_mu[float(Ti)] = {}

    data = read_file(fname.strip(),ftype)
    print(Ti,mu,np.mean(data[1:len(data),6],axis=0))
    dCount,gCount = 0,0
    for k in range(len(data)):
        if math.isnan(data[k][0]):
            gCount += 1
        else:
            dCount += 1
    if np.mean(data[1:len(data),6],axis=0)  < 16/L**2.0:
        for k in range(len(data)):
            data[k][0] = data[k][4]
            data[k][1] = data[k][5]
            data[k][2] = data[k][4]
            data[k][3] = data[k][5]
            data[k][6] = np.NaN
            data[k][7] = np.NaN
            data[k][10] = data[k][14]
            data[k][11] = data[k][14]
            if float(Ti) == 2.0:
                data[k][12] = -1
                data[k][13] = -1
            else:
                data[k][12] = 0
                data[k][13] = 0
    if len(data) == 0:
        data = np.NaN*np.ones(15)
    else:
        data_by_mu[float(Ti)][float(mu)] = data[1:]

        print(len(data_by_mu[float(Ti)][float(mu)]),np.mean(data[:,6]))
Mus = sorted(data_by_mu[10].keys())
print("Mus",Mus)
norm = 1#np.nanmean(data_by_mu[2.0][Mus[0]][:][5])
plot_ads(data_by_mu,norm)
plot_area(data_by_mu)
plot_surface_data(data_by_mu)

