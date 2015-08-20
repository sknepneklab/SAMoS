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

# Utility code for computing mean square displacement (MSD) on a sphere

import sys
import argparse
from read_data import *
from math import *
from glob import glob
import numpy as np
from datetime import *
import vtk
from vtk import *
 
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="base name of the input files")
parser.add_argument("-o", "--output", type=str, default='msd.dat', help="MSD data")
parser.add_argument("-s", "--skip", type=int, default=0, help="skip this many samples")
parser.add_argument("-r", "--radius", type=float, default=1.0, help="sphere radius")
parser.add_argument("-S", "--step", type=float, default=5000, help="Time step multiplier")
args = parser.parse_args()

print
print "\tActive Particles on Curved Spaces (APCS)"
print "\tTrajectory of each particle"
print 
print "\tRastko Sknepnek"
print "\tUniversity of Dundee"
print "\t(c) 2014"
print "\t----------------------------------------------"
print
print "\tInput files : ", args.input
print "\tOutput files : ", args.output
print "\tSkip frames : ", args.skip
print "\tSphere radius : ", args.radius
print "\tTime step multiplier : ", args.step
print

start = datetime.now()

files = sorted(glob(args.input+'*.dat'))[args.skip:]

if len(files) == 0:
  files = sorted(glob(args.input+'*.dat.gz'))[args.skip:]


init_data = ReadData(files[0])
x0 = np.array(init_data.data[init_data.keys['x']])
r = np.zeros((x0.size,3,len(files)))


dist_dat = []

i = 0
for f in files:
  print "Processing file : ", f
  data = ReadData(f)
  x, y, z = np.array(data.data[data.keys['x']]), np.array(data.data[data.keys['y']]), np.array(data.data[data.keys['z']])
  rr = np.vstack((x,y,z)).T
  r[:,:,i] = rr
  #dist = args.radius*np.arccos(np.einsum('ij,ij->i',r,r0)/args.radius**2)
  #dist_dat.append(dist**2)
  i += 1

print 'Computing MSD...'

all_dist = np.zeros((x0.size,len(files)))
all_dist_r = np.zeros((x0.size,len(files)))
for tau in xrange(1,len(files)):
  print 'Computing for tau : ', tau
  vals = np.nan_to_num(args.radius*np.arccos(np.einsum('ijk,ijk->ik',r[:,:,:-tau],r[:,:,tau:])/args.radius**2))
  dr = (r[:,:,:-tau]-r[:,:,tau:])**2
  ddr = np.sum(dr,axis=1)
  dist = np.mean(vals,axis=1)
  dist_r = np.mean(ddr,axis=1)
  all_dist[:,tau] = dist**2
  all_dist_r[:,tau] = dist_r
  

msd = np.mean(all_dist,axis=0)
msd_std = np.std(all_dist,axis=0)
msd_r = np.mean(all_dist_r,axis=0)
msd_r_std = np.std(all_dist_r,axis=0)

t = args.step*np.arange(len(files))

out = open(args.output,'w')

for tt,m,ms,m_r,m_rs in zip(t,msd,msd_std,msd_r,msd_r_std):
  out.write('%f  %f  %f  %f  %f\n' % (tt,m,ms,m_r,m_rs))

out.close()

end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print


