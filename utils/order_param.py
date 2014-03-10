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

# Utility code for computing (time) average of the order parameter
# for a series of director files.

import sys
import argparse
from read_data import *
from math import *
from glob import glob
import numpy as np
from datetime import *

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", type=str, default='hist.dir', help="velocity histogram")
parser.add_argument("-i", "--input", type=str, help="base name of the input files")
parser.add_argument("-s", "--skip", type=int, default=0, help="skip this many samples")
parser.add_argument("-t", "--time", type=str, default='time_seq.dat', help="record time sequence of the order parameter")
parser.add_argument("-T", "--type", type=str, default='sphere', help="system type (plane or sphere)")
args = parser.parse_args()

print
print "\tActive Particles on Curved Spaces (APCS)"
print "\tComputes order parameter for spherical systems"
print 
print "\tRastko Sknepnek"
print "\tUniversity of Dundee"
print "\t(c) 2014"
print "\t----------------------------------------------"
print
print "\tInput files : ", args.input
print "\tOutput files : ", args.output
print "\tSkip frames : ", args.skip
print "\tStore OP time sequence in : ", args.time
print "\tSystem type : ", args.type
print

start = datetime.now()

files = sorted(glob(args.input+'*.dat'))[args.skip:]

time_file = open(args.time,'w')

op = []
idx = 0
for f in files:
  print "Processing file : ", f
  data = ReadData(f)
  timestep = int(f.split('_')[-1].split('.')[0])
  
  if data.has_header:
    vx, vy, vz = np.array(data.data[data.keys['vx']]), np.array(data.data[data.keys['vy']]), np.array(data.data[data.keys['vz']])
    v = np.vstack((vx,vy,vz)).transpose()
    if args.type == "sphere":
      x, y, z = np.array(data.data[data.keys['x']]), np.array(data.data[data.keys['y']]), np.array(data.data[data.keys['z']])
      r = np.vstack((x,y,z)).transpose()
      len_r = np.apply_along_axis(np.linalg.norm,1,r).reshape((len(r),1))
      len_v = np.apply_along_axis(np.linalg.norm,1,v).reshape((len(v),1))
      r = r/len_r
      v = v/len_v
      loc_op = np.cross(r,v)
      mean_loc_op = np.linalg.norm(np.apply_along_axis(sum,0,loc_op))/data.N
    else:
      mean_loc_op = np.linalg.norm(np.apply_along_axis(sum,0,v))/data.N
    time_file.write('%d  %f\n' % (timestep,mean_loc_op))
    op.append(mean_loc_op)

time_file.close()

out = open(args.output,'w')
out.write('#  op  std_error\n')
out.write('%f  %f\n' % (np.mean(op), np.std((op))/sqrt(len(op)-1)))
out.close()

end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print


