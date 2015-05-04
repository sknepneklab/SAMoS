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
#    (c) 2013
#    
#    This program cannot be used, copied, or modified without
#    explicit permission of the author.
#  
# ################################################################


# Computes average velocity histogram

import sys
import argparse
from read_data import *
from math import *
from glob import glob
from numpy import histogram
from datetime import *

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--min", type=float, default=0.0, help="minimum velocity")
parser.add_argument("-m", "--max", type=float, default=3.0, help="maximum velocity")
parser.add_argument("-b", "--bin", type=int, default=50, help="number of bins")
parser.add_argument("-o", "--output", type=str, default='hist.dir', help="velocity histogram")
parser.add_argument("-i", "--input", type=str, help="base name of the input files")
parser.add_argument("-s", "--skip", type=int, default=0, help="skip this many samples")
args = parser.parse_args()

print
print "\tActive Particles on Curved Spaces (APCS)"
print "\tAverage velocity distribution"
print 
print "\tRastko Sknepnek"
print "\tUniversity of Dundee"
print "\t(c) 2014"
print "\t----------------------------------------------"
print
print "\tMinmum velocity : ", args.min
print "\tMaximum velocity : ", args.max
print "\tNumber of bins : ", args.bin
print

start = datetime.now()

files = sorted(glob(args.input+'*.dat'))[args.skip:]

vel = []
for f in files:
  data = ReadData(f)

  if data.has_header:
    for (vx, vy, vz) in zip(data.data[data.keys['vx']],data.data[data.keys['vy']],data.data[data.keys['vz']]):
      vel.append(sqrt(vx*vx + vy*vy + vz*vz))

hist = histogram(vel,bins=args.bin,range=(args.min,args.max))

out = open(args.output,'w')
out.write('#  v  num\n')
for (n,v) in zip(hist[0],hist[1]):
  out.write('%f  %d\n' % (v, n))
out.close()

end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print
 