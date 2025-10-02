#!?usr/bin/python
import sys
import os
import numpy as np

submit_script = sys.argv[1]
# parameters for prewetting transition
# @ u = -4.5, prewetting occurs ~ hphi =1,cbound = 0.06,Jbulk > 0.6


cbound = 0.04
s = 1
#### Fig 4
Jbulk = 1.6
mu1Arr = np.arange(-6,-2.6,0.1)
JmemArr = [1.05]
mem_comp_arr = [0,0.5]
jobfile = open('job_file_FigS4.sh','w')
for mu1 in mu1Arr:
    for fidx in range(0,5):
        jobfile.write('python %s %0.3f %0.3f %0.3f %0.3f %d %0.4f\n' % (submit_script,cbound,Jbulk,mu1,mu1,fidx,s))                
