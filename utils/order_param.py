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
from numpy import histogram, cross, mean, std
from numpy.linalg import norm
from datetime import *

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", type=str, default='hist.dir', help="velocity histogram")
parser.add_argument("-i", "--input", type=str, help="base name of the input files")
parser.add_argument("-s", "--skip", type=int, default=0, help="skip this many samples")
parser.add_argument("-t", "--time", type=str, default='time_seq.dat', help="record time sequence of the order parameter")
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
print

start = datetime.now()

files = sorted(glob(args.input+'*.dat'))[args.skip:]

time_file = open(args.time,'w')

op = []
idx = 0
for f in files:
  print "Processing file : ", f
  data = ReadData(f)
  
  if data.has_header:
    v = [0,0,0]
    for (x, y, z, nx, ny, nz) in zip(data.data[data.keys['x']],data.data[data.keys['y']],data.data[data.keys['z']],data.data[data.keys['nx']],data.data[data.keys['ny']],data.data[data.keys['nz']]):
      ln = sqrt(x*x+y*y+z*z)
      lnn = sqrt(nx*nx + ny*ny + nz*nz)
      v = map(lambda x, y: x + y, v, cross([x/ln,y/ln,z/ln],[nx/lnn,ny/lnn,nz/lnn]))
    time_file.write('%d  %f  %f  %f  %f\n' % (idx, norm(v)/data.N, v[0]/data.N, v[1]/data.N, v[2]/data.N))
  idx += 1
  op.append(norm(v)/data.N)

time_file.close()

out = open(args.output,'w')
out.write('#  op  std_error\n')
out.write('%f  %f\n' % (mean(op), std((op))/sqrt(len(op)-1)))
out.close()

end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print


