# * *************************************************************
# *  
# *   Soft Active Mater on Surfaces (SAMoS)
# *   
# *   Author: Rastko Sknepnek
# *  
# *   Division of Physics
# *   School of Engineering, Physics and Mathematics
# *   University of Dundee
# *   
# *   (c) 2013, 2014
# * 
# *   School of Science and Engineering
# *   School of Life Sciences 
# *   University of Dundee
# * 
# *   (c) 2015
# * 
# *   Author: Silke Henkes
# * 
# *   Department of Physics 
# *   Institute for Complex Systems and Mathematical Biology
# *   University of Aberdeen  
# * 
# *   (c) 2014, 2015
# *  
# *   This program cannot be used, copied, or modified without
# *   explicit written permission of the authors.
# * 
# * ***************************************************************

# Utility code for converting standard data output (dat) of the simulation
# into the VTP format suitable for visualization with ParaView

import matplotlib
matplotlib.use('Agg')


import sys
import argparse
import pickle
import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
import math as m
from datetime import *

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="input file (defect data columns format)")
parser.add_argument("-o", "--output", type=str, default='out.png', help="output name for the figure (PNG format)")
parser.add_argument("-d", "--data", type=str, default='out.dat', help="textual output")
parser.add_argument("-a", "--autocorrel", type=str, default='autocorrel.dat', help="name of the autocorrelation data file")
parser.add_argument("-H", "--histogram", type=str, default='hist.dat', help="name of the time histogram file")
parser.add_argument("-D", "--defname", type=str, default='def.dat', help="name of the defect number vs. time file")
parser.add_argument("-b", "--bins", type=int, default=25, help="number of bins in the histogram")
parser.add_argument("-s", "--skip", type=int, default=0, help="skip this many time steps")
parser.add_argument("-S", "--step", type=float, default=5000, help="Time step multiplier")
parser.add_argument("-R", "--radius", type=float, default=16.0, help="System's radius")
parser.add_argument("--v0", type=float, default=1.0, help="v0")
parser.add_argument("--tau", type=float, default=1.0, help="tau flip")
parser.add_argument("-J", "--J", type=float, default=1.0, help="J")
args = parser.parse_args()

v0 =  args.v0
R = args.radius
J = args.J
tau_flip = args.tau


print
print "\tActive Particles on Curved Spaces (APCS)"
print "\tConverts dat files to VTP files"
print 
print "\tRastko Sknepnek"
print "\tUniversity of Dundee"
print "\t(c) 2014"
print "\t----------------------------------------------"
print
print "\tInput file : ", args.input
print "\tOutput file : ", args.output
print "\tData file : ", args.data
print "\tHistogram file : ", args.histogram
print "\tAutocorrelation data file : ", args.autocorrel
print "\tSkip timesteps : ", args.skip
print "\tTime steps per snapshots : ", args.step
print "\tNumber of bins : ", args.bins
print "\tR : ", R
print "\tv0 : ", v0
print "\tJ : ", J
print "\ttau flip : ", args.tau

start = datetime.now()

ext = args.input.split('.')[-1]

if ext == 'p':
  inp = open(args.input,'rb')
  inp_data = pickle.load(inp)
  data = []
  for d in inp_data['defects_n']:
    if len(d) >= 4:
      line = [len(d)]
      count = 0
      for df in d:
        if count < 4:
          line.extend(df[1:])
        count += 1
      data.append(line)
else:
  inp = open(args.input,'r')
  data = map(lambda x: map(float, x.strip().split()),inp.readlines()[1:])

steps = args.step*np.array(range(len(data)))

datamat = np.matrix(data[args.skip:])
steps = steps[args.skip:] - args.skip*args.step

print "\tTotal number of samples : ", len(steps)
print

ndef = np.ravel(np.array(datamat[:,0]))
d1 = np.array(datamat[:,1:4]/lin.norm(datamat[0,1:4]))
d2 = np.array(datamat[:,4:7]/lin.norm(datamat[0,4:7]))
d3 = np.array(datamat[:,7:10]/lin.norm(datamat[0,7:10]))
d4 = np.array(datamat[:,10:14]/lin.norm(datamat[0,10:14]))


a12 = np.degrees(np.arccos(np.sum(d1.T*d2.T,0)))
a13 = np.degrees(np.arccos(np.sum(d1.T*d3.T,0)))
a14 = np.degrees(np.arccos(np.sum(d1.T*d4.T,0)))
a23 = np.degrees(np.arccos(np.sum(d2.T*d3.T,0)))
a24 = np.degrees(np.arccos(np.sum(d2.T*d4.T,0)))
a34 = np.degrees(np.arccos(np.sum(d3.T*d4.T,0)))
angles = (a12 + a13 + a14 + a23 + a24 + a34)/6.0
allangles = np.concatenate((a12,a13,a14,a23,a24,a34))

A = np.zeros((steps.size,6))
for i in xrange(steps.size):
  a = np.array(np.sort([a12[i],a13[i],a14[i],a23[i],a24[i],a34[i]]))
  A[i] = a

datout = open(args.data,'w')
datout.write('%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n' % (R, v0, np.mean(angles), np.std(angles), np.mean(A[:,0]), np.std(A[:,0]), np.mean(A[:,1]), np.std(A[:,1]), np.mean(A[:,2]), np.std(A[:,2]), np.mean(A[:,3]), np.std(A[:,3]),np.mean(A[:,4]), np.std(A[:,4]),np.mean(A[:,5]), np.std(A[:,5])))
print '<a> = ',np.mean(angles), '+/-', np.std(angles)
print '<a1> = ',np.mean(A[:,0]), '+/-', np.std(A[:,0])
print '<a2> = ',np.mean(A[:,1]), '+/-', np.std(A[:,1])
print '<a3> = ',np.mean(A[:,2]), '+/-', np.std(A[:,2])
print '<a4> = ',np.mean(A[:,3]), '+/-', np.std(A[:,3])
print '<a5> = ',np.mean(A[:,4]), '+/-', np.std(A[:,4])
print '<a6> = ',np.mean(A[:,5]), '+/-', np.std(A[:,5])
datout.close()

ang = angles - np.mean(angles)
ang_correl = np.zeros(ang.size)
for tau in xrange(ang.size):
  ang_correl[tau] = np.mean(ang[:ang.size-tau]*ang[tau:])
ang_correl /= ang_correl[0]


autodata = open(args.autocorrel,'w')
for (step,ac) in zip(steps,ang_correl):
  autodata.write('%f  %f\n' % (step,ac))
autodata.close()

plt.rc('text', usetex=True)
plt.xlim((0,max(steps[:ang_correl.size/2])))
plt.ylim((1.1*min(ang_correl),1.0*max(ang_correl)))
plt.title(r'Autocorrelation function R = '+str(R)+', $v_0$ = '+str(v0)+' J = '+str(J))
plt.plot(steps[:ang_correl.size/2],ang_correl[:ang_correl.size/2],'b-', linewidth=2.0)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.xlabel(r'$\tau$',fontsize=20)
plt.ylabel(r'autocorrelation function',fontsize=20)
autocorrelname = '.'.join(args.autocorrel.split('.')[:-1])
plt.savefig(autocorrelname+'.png',format='png')
plt.clf()


hist = np.histogram(angles,bins=args.bins)
bins = hist[1]
vals = hist[0].astype(np.float32)
histdata = open(args.histogram,'w')
for (b,v) in zip(bins,vals):
  histdata.write('%f  %f\n' % (b,v))
histdata.close()


plt.rc('text', usetex=True)
plt.xlim((0,180))
#plt.ylim((0,200))
plt.title(r'Histogram of all $\alpha$, R = '+str(R)+', $v_0$ = '+str(v0)+' J = '+str(J))
bins = np.linspace(0,180,90)
plt.hist(allangles,bins,color='green',alpha=0.8,normed=1)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.xlabel(r'$\alpha$',fontsize=20)
plt.ylabel(r'$P(\alpha)$',fontsize=20)
histname = '.'.join(args.histogram.split('.')[:-1])
plt.savefig(histname+'.png',format='png')
plt.clf()

plt.rc('text', usetex=True)
plt.xlim((0,180))
#plt.ylim((0,200))
plt.title(r'Histogram of $<\alpha>$, R = '+str(R)+', $v_0$ = '+str(v0)+' J = '+str(J))
bins = np.linspace(0,180,90)
plt.hist(angles,bins,color='green',alpha=0.8,normed=1)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.xlabel(r'$<\alpha>$',fontsize=20)
plt.ylabel(r'$P(<\alpha>)$',fontsize=20)
histname = '.'.join(args.histogram.split('.')[:-1])
plt.savefig('mean_'+histname+'.png',format='png')
plt.clf()


end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print

plt.rc('text', usetex=True)
plt.xlim((0,max(steps)))
plt.ylim((0,180.0))
plt.title(r'$\bar{\alpha} = $'+str(round(np.mean(angles),1))+'$\pm$'+str(round(np.std(angles),1))+' for R = '+str(R)+', $v_0$ = '+str(v0)+' J = '+str(J)+', $\\tau_{flip}$ = '+str(tau_flip))
plt.plot(steps,A[:,0],'r-', label=r'$\alpha_{1}$')
plt.plot(steps,A[:,1],'g-', label=r'$\alpha_{2}$')
plt.plot(steps,A[:,2],'b-', label=r'$\alpha_{3}$')
plt.plot(steps,A[:,3],'y-', label=r'$\alpha_{4}$')
plt.plot(steps,A[:,4],'c-', label=r'$\alpha_{5}$')
plt.plot(steps,A[:,5],'m-', label=r'$\alpha_{6}$')
plt.plot(steps,angles,'k-',linewidth=3, label=r'$\left<\alpha\right>$')
plt.tick_params(axis='both', which='major', labelsize=20)
plt.xlabel('time step',fontsize=20)
plt.ylabel(r'$\left<\alpha\right>({}^\circ)$',fontsize=20)
plt.legend(loc=4)
plt.savefig(args.output,format='png')
plt.clf()

plt.rc('text', usetex=True)
plt.xlim((0,max(steps)))
plt.ylim((0,20))
plt.title(r'Number of defects vs. time for R = '+str(R)+', $v_0$ = '+str(v0)+' J = '+str(J))
plt.plot(steps,ndef,'r-', linewidth=2)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.xlabel('time step',fontsize=20)
plt.ylabel(r'$N_{def}$',fontsize=20)
defname = '.'.join(args.defname.split('.')[:-1])
plt.savefig(defname+'.png',format='png')

#plt.show()


