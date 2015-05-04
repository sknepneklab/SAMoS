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

# Utility code for computing time correlation function for
# velocity (or director) in plane for a series of velocity/director files.

import sys
import argparse
from read_data import *
from math import *
from glob import glob
from numpy import histogram, cross, mean, std, dot
from numpy.linalg import norm
from datetime import *




parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", type=str, default='time_correl.dat', help="contains time correlation data")
parser.add_argument("-i", "--input", type=str, help="base name of the input files")
parser.add_argument("-s", "--skip", type=int, default=0, help="skip this many samples")
parser.add_argument("-d", "--director", action="store_true", help="compute time correlations of the director (velocity is default)")
args = parser.parse_args()

print
print "\tActive Particles on Curved Spaces (APCS)"
print "\tComputes time correlation function of velocity"
print "\tor director field on plane"
print 
print "\tRastko Sknepnek"
print "\tUniversity of Dundee"
print "\t(c) 2014"
print "\t----------------------------------------------"
print
print "\tInput files : ", args.input
print "\tOutput files : ", args.output
print "\tSkip frames : ", args.skip
if args.director:
  print "\tCompute time correlation of the director field."
print

start = datetime.now()

files = sorted(glob(args.input+'*.dat'))[args.skip:]

if args.director:
  qx, qy, qz = 'nx', 'ny', 'nz'
else:
  qx, qy, qz = 'vx', 'vy', 'vz'

vals = []
idx = 0
for f in files:
  print "Processing file : ", f
  data = ReadData(f)
  frame = []
  if data.has_header:
    for (ax, ay, az) in zip(data.data[data.keys[qx]],data.data[data.keys[qy]],data.data[data.keys[qz]]):
      frame.append([ax,ay,az])
    vals.append([idx,frame])  
  idx += 1

time_file = open(args.output,'w')

print "Computing averages..."

for tau in range(0,len(vals)-1):
  print "Processing tau : ", tau
  avg_tau = 0.0
  for t in range(0,len(vals)-tau):
    frame1 = vals[t][1]
    frame2 = vals[t+tau][1]
    avg = 0.0
    for (v1,v2) in zip(frame1,frame2):
      avg += dot(v1,v2)
    avg /= len(frame1)
    avg_tau += avg
  time_file.write('%d   %f\n' % (tau,avg_tau/(len(vals)-tau)))
    
time_file.close()
  


end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print
