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
# into the VTK format suitable for visualization with ParaView

import sys
import argparse
from read_data import *
import numpy as np
from datetime import *

from vtktools import *

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="input file")
parser.add_argument("-o", "--output", type=str, help="output file")
parser.add_argument("-d", "--distance", type=str, default=None, help="file with interparticle distances")
args = parser.parse_args()

print
print "\tActive Particles on Curved Spaces (APCS)"
print "\tConverts dat files to VTK files"
print 
print "\tRastko Sknepnek"
print "\tUniversity of Dundee"
print "\t(c) 2014"
print "\t----------------------------------------------"
print
print "\tInput files : ", args.input
print "\tOutput files : ", args.output
if args.distance != None:
  print "\tReading distances to 0th particle from : ", args.distance
print

start = datetime.now()

vtk_writer = VTK_XML_Serial_Unstructured()
has_v = False
has_n = False

data = ReadData(args.input)

dist = []
if args.distance != None:
  distinp = open(args.distance,'r')
  dist = distinp.readlines()
  dist = map(lambda x: x.strip(),dist)
  dist = map(float,dist)
  distinp.close()

if not (data.keys.has_key('x') and data.keys.has_key('y') and data.keys.has_key('z')):
  raise "Particle coordinate not specified in the input data."

x = np.array(data.data[data.keys['x']])
y = np.array(data.data[data.keys['y']])
z = np.array(data.data[data.keys['z']])

data_to_print = {}

if (data.keys.has_key('vx') or data.keys.has_key('vy') or data.keys.has_key('vz')):
  vx = np.array(data.data[data.keys['vx']])
  vy = np.array(data.data[data.keys['vy']])
  vz = np.array(data.data[data.keys['vz']])
  data_to_print['vx'] = vx
  data_to_print['vy'] = vy
  data_to_print['vz'] = vz
  has_v = True

if (data.keys.has_key('nx') or data.keys.has_key('ny') or data.keys.has_key('nz')):
  nx = np.array(data.data[data.keys['nx']])
  ny = np.array(data.data[data.keys['ny']])
  nz = np.array(data.data[data.keys['nz']])
  data_to_print['nx'] = nx
  data_to_print['ny'] = ny
  data_to_print['nz'] = nz
  has_n = True

r = np.ones(len(x))  
  
if has_v and has_n:
  if len(dist) > 0:
    vtk_writer.snapshot(args.output, x, y, z, vx=vx, vy=vy, vz=vz, nx=nx, ny=ny, nz=nz,radii=r,dist=dist)
  else:
    vtk_writer.snapshot(args.output, x, y, z, vx=vx, vy=vy, vz=vz, nx=nx, ny=ny, nz=nz,radii=r)
elif has_v:
  vtk_writer.snapshot(args.output, x, y, z, vx=vx, vy=vy, vz=vz,radii=r)
elif has_n:
  vtk_writer.snapshot(args.output, x, y, z, nx=nx, ny=ny, nz=nz,radii=r)
else:
  vtk_writer.snapshot(args.output, x, y, z, radii=r)

end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print
  
