import sys
import argparse
from read_data import *
from math import *
from glob import glob
import numpy as np
from datetime import *
from vtktools import *

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

# Utility code for computing distances between seed (0th) particle 
# and the rest of the particles

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="input file in ACPS dat format")
parser.add_argument("-o", "--output", type=str, default='dist.dat', help="list of all distances")
args = parser.parse_args()

print
print "\tActive Particles on Curved Spaces (APCS)"
print "\tDistance to seed particle"
print 
print "\tRastko Sknepnek"
print "\tUniversity of Dundee"
print "\t(c) 2014"
print "\t----------------------------------------------"
print
print "\tInput files : ", args.input
print "\tOutput files : ", args.output
print

start = datetime.now()


data = ReadData(args.input)

x, y, z = np.array(data.data[data.keys['x']]), np.array(data.data[data.keys['y']]), np.array(data.data[data.keys['z']])
x0, y0, z0 = x[0], y[0], z[0]

dist = np.sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2) 

out = open(args.output,'w')

for d in dist:
  out.write('%f\n' % d)

out.close()


end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print
