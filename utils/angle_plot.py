# ################################################################
#  
#    Active Particles on Curved Spaces (APCS)
#    
#    Author: Rastko Sknepnek
#   
#    Division of Physics
#    School of Engineering, Physics and Mathematics
#    University of Dundee
#    
#    (c) 2013, 2014
#    
#    This program cannot be used, copied, or modified without
#    explicit permission of the author.
#  
# ################################################################

# Utility code for converting standard data output (dat) of the simulation
# into the VTP format suitable for visualization with ParaView

import sys
import argparse
from read_data import *
import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
import math as m
from datetime import *


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="input file (defect data columns format)")
parser.add_argument("-o", "--output", type=str, help="output directory")
parser.add_argument("-s", "--skip", type=int, default=0, help="skip this many time steps")
parser.add_argument("-S", "--step", type=float, default=5000, help="Time step multiplier")
args = parser.parse_args()

print
print "\tActive Particles on Curved Spaces (APCS)"
print "\tConverts dat files to VTP files"
print 
print "\tRastko Sknepnek"
print "\tUniversity of Dundee"
print "\t(c) 2014"
print "\t----------------------------------------------"
print
print "\tInput files : ", args.input
print "\tOutput files : ", args.output
print "\tSkip timesteps : ", args.skip
print

start = datetime.now()

inp = open(args.input,'r')
data = map(lambda x: map(float, x.strip().split()),inp.readlines()[1:])

steps = args.step*np.array(range(len(data)))

datamat = np.matrix(data[args.skip:])
steps = steps[args.skip:]


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

A = np.zeros((steps.size,6))
for i in xrange(steps.size):
  a = np.array(np.sort([a12[i],a13[i],a14[i],a23[i],a24[i],a34[i]]))
  A[i] = a


print '<a> = ',np.mean(angles), '+/-', np.std(angles)
print '<a1> = ',np.mean(A[:,0]), '+/-', np.std(A[:,0])
print '<a2> = ',np.mean(A[:,1]), '+/-', np.std(A[:,1])
print '<a3> = ',np.mean(A[:,2]), '+/-', np.std(A[:,2])
print '<a4> = ',np.mean(A[:,3]), '+/-', np.std(A[:,3])
print '<a5> = ',np.mean(A[:,4]), '+/-', np.std(A[:,4])
print '<a6> = ',np.mean(A[:,5]), '+/-', np.std(A[:,5])

end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print

plt.rc('text', usetex=True)
plt.xlim((0,max(steps)))
plt.plot(steps,A[:,0],'r-', label=r'$\alpha_{1}$')
plt.plot(steps,A[:,1],'g-', label=r'$\alpha_{2}$')
plt.plot(steps,A[:,2],'b-', label=r'$\alpha_{3}$')
plt.plot(steps,A[:,3],'y-', label=r'$\alpha_{4}$')
plt.plot(steps,A[:,4],'c-', label=r'$\alpha_{5}$')
plt.plot(steps,A[:,5],'m-', label=r'$\alpha_{6}$')
plt.plot(steps,angles,'k-',linewidth=3, label=r'$\left<\alpha\right>$')
plt.tick_params(axis='both', which='major', labelsize=20)
plt.xlabel('time step',fontsize=20)
plt.ylabel(r'$\left<\alpha_{min}\right>({}^\circ)$',fontsize=20)
plt.legend(loc=4)
plt.show()


