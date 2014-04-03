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

# Utility code for computing time dependence of the distance between 
# two arbitrary particles

import sys
import argparse
from read_data import *
from math import *
from glob import glob
import numpy as np
from datetime import *

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="base name of the input files")
parser.add_argument("-o", "--output", type=str, default='dist_p1_XXX_p2_YYY.dat', help="distance between particles XXX and YYY as a function of time")
parser.add_argument("-s", "--skip", type=int, default=0, help="skip this many samples")
parser.add_argument("-p", "--p1", type=int, default=0, help="id of the 1st particle to track")
parser.add_argument("-q", "--p2", type=int, default=1, help="id of the 2nd particle to track")
args = parser.parse_args()

print
print "\tActive Particles on Curved Spaces (APCS)"
print "\tInterparticle distance over time"
print 
print "\tRastko Sknepnek"
print "\tUniversity of Dundee"
print "\t(c) 2014"
print "\t----------------------------------------------"
print
print "\tInput files : ", args.input
print "\tOutput files : ", args.output
print "\tSkip frames : ", args.skip
print "\tId of 1st particle : ", args.p1
print "\tId of 2nd particle : ", args.p2
print

start = datetime.now()

files = sorted(glob(args.input+'*.dat'))[args.skip:]

out = open(args.output+'.dat','w')
p1 = args.p1
p2 = args.p2


out.write('# Interparticle distance between : %d and %d\n' % (p1,p2))
out.write('# timestep dist n1_dot_n2  v1_dot_v2\n')  

idx = 0
for f in files:
  print "Processing file : ", f
  data = ReadData(f)
  timestep = int(f.split('_')[-1].split('.')[0])
  
  if data.has_header:
    x, y, z = np.array(data.data[data.keys['x']]), np.array(data.data[data.keys['y']]), np.array(data.data[data.keys['z']])
    vx, vy, vz = np.array(data.data[data.keys['vx']]), np.array(data.data[data.keys['vy']]), np.array(data.data[data.keys['vz']])
    nx, ny, nz = np.array(data.data[data.keys['nx']]), np.array(data.data[data.keys['ny']]), np.array(data.data[data.keys['nz']])
    v = np.vstack((vx,vy,vz)).transpose()
    len_v = np.apply_along_axis(np.linalg.norm,1,v).reshape((len(v),1))
    x1, y1, z1 = x[p1], y[p1], z[p1]
    x2, y2, z2 = x[p2], y[p2], z[p2]
    d = sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
    n1x, n1y, n1z = nx[p1], ny[p1], nz[p1]
    n2x, n2y, n2z = nx[p2], ny[p2], nz[p2]
    n_dot_n = n1x*n2x + n1y*n2y + n1z*n2z
    v1x, v1y, v1z = vx[p1], vy[p1], vz[p1]
    v2x, v2y, v2z = vx[p2], vy[p2], vz[p2]
    v_dot_v = v1x*v2x + v1y*v2y + v1z*v2z
    out.write('%d %f %f %f\n' % (idx, d, n_dot_n, v_dot_v))
  idx += 1


out.close()

end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print
