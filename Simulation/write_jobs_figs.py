#!?usr/bin/python
import sys
import os
import numpy as np

#### Generates parameters  and commands for all plots used in the main text and supplement #######
submit_script = sys.argv[1]
# @ u = -4.5, prewetting occurs ~ hphi =1,cbound = 0.06,Jbulk > 0.6


cbound = 0.04
s = 1
mu1Arr=[-4.5,-4.0,-3.1,-2.9,-2.8,-2.7]
############# Figure 1 ######################
Jbulk= 1.6
mem_comp_arr = [0]
jobfile = open('job_file_Fig1.sh','w')
for mem_comp in mem_comp_arr:
    for mu1 in mu1Arr:
        Tm = 10
        for fidx in range(1,5):
            jobfile.write('python %s %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %d %0.4f\n' % (submit_script,cbound,Jbulk,mu1,mu1,mem_comp,Tm,fidx,s))
jobfile.close()
                
JbulkArr = 1.6
mu1Arr = [-3.1,-2.9]
JmemArr = [1.25,0.9]
mem_comp_arr = [0,0.5]
jobfile = open('job_file_Fig2C.sh','w')
for mem_comp in mem_comp_arr:
    
    for mu1 in mu1Arr:
        if mem_comp != 0:
            for Tm in JmemArr:
                for fidx in range(0,10):
                    jobfile.write('python %s %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %d %0.4f\n' % (submit_script,cbound,Jbulk,mu1,mu1,mem_comp,Tm,fidx,s))                
        else:
            Tm = 10
            for fidx in range(0,10):
                jobfile.write('python %s %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %d %0.4f\n' % (submit_script,cbound,Jbulk,mu1,mu1,mem_comp,Tm,fidx,s))
jobfile.close()


                        
Jbulk= 1.6
mu1Arr = np.arange(-6,-1.9,0.1)
JmemArr = [1.05]
mem_comp_arr = [0,0.5]
jobfile = open('job_file_FigS2C.sh','w')
for mem_comp in mem_comp_arr:
    for mu1 in mu1Arr:
        if mem_comp != 0:
            for Tm in JmemArr:
                for fidx in range(0,5):
                    jobfile.write('python %s %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %d %0.4f\n' % (submit_script,cbound,Jbulk,mu1,mu1,mem_comp,Tm,fidx,s))                
        else:
            Tm = 10
            for fidx in range(0,5):
                jobfile.write('python %s %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %d %0.4f\n' % (submit_script,cbound,Jbulk,mu1,mu1,mem_comp,Tm,fidx,s))

jobfile.close()

#### Fig 4
Jbulk = 1.6
mu1Arr = np.arange(-6,-2.6,0.1)
JmemArr = [1.05]
mem_comp_arr = [0,0.5]
jobfile = open('job_file_FigS4.sh','w')
for mem_comp in mem_comp_arr:
    for mu1 in mu1Arr:
        if mem_comp != 0:
            for Tm in JmemArr:
                for fidx in range(0,5):
                    jobfile.write('python %s %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %d %0.4f\n' % (submit_script,cbound,Jbulk,mu1,mu1,mem_comp,Tm,fidx,s))                
        else:
            Tm = 10
            for fidx in range(0,5):
                jobfile.write('python %s %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %d %0.4f\n' % (submit_script,cbound,Jbulk,mu1,mu1,mem_comp,Tm,fidx,s))
jobfile.close()
