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

# Utility code for computing trajectory of a given particle

import sys
import argparse
from read_data import *
from math import *
from glob import glob
import numpy as np
from datetime import *
from vtktools import *

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="base name of the input files")
parser.add_argument("-o", "--output", type=str, default='trajXXX.dat', help="particle trajectory (XXXX is particle number)")
parser.add_argument("-s", "--skip", type=int, default=0, help="skip this many samples")
parser.add_argument("-n", "--id", type=int, default=0, help="id of the particle to track")
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
print

start = datetime.now()

files = sorted(glob(args.input+'*.dat'))[args.skip:]

out = open(args.output+'.dat','w')
pid = args.id


xx = []
yy = []
zz = []
vvx = []
vvy = []
vvz = []
nnx = []
nny = []
nnz = []
op = []
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
    xx.append(x[pid])
    yy.append(y[pid])
    zz.append(z[pid])
    vvx.append(vx[pid])
    vvy.append(vy[pid])
    vvz.append(vz[pid])
    nnx.append(nx[pid])
    nny.append(ny[pid])
    nnz.append(nz[pid])
    out.write('%f %f %f %f  %f  %f  %f\n' % (x[pid],y[pid],z[pid],len_v[pid],vx[pid],vy[pid],vz[pid]))

out.close()

vtk_writer = VTK_XML_Serial_Unstructured()

vtk_writer.snapshot(args.output+'.vtu',xx,yy,zz,vx=vvx,vy=vvy,vz=vvz,nx=nnx,ny=nny,nz=nnz)


end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print


