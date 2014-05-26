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
from datetime import *
import vtk


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="input file")
parser.add_argument("-o", "--output", type=str, help="output file")
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
print

start = datetime.now()

Points = vtk.vtkPoints()


has_v = False
has_n = False

data = ReadData(args.input)

if not (data.keys.has_key('x') and data.keys.has_key('y') and data.keys.has_key('z')):
  raise "Particle coordinate not specified in the input data."

x = np.array(data.data[data.keys['x']])
y = np.array(data.data[data.keys['y']])
z = np.array(data.data[data.keys['z']])


if (data.keys.has_key('vx') or data.keys.has_key('vy') or data.keys.has_key('vz')):
  vx = np.array(data.data[data.keys['vx']])
  vy = np.array(data.data[data.keys['vy']])
  vz = np.array(data.data[data.keys['vz']])
  has_v = True

if (data.keys.has_key('nx') or data.keys.has_key('ny') or data.keys.has_key('nz')):
  nx = np.array(data.data[data.keys['nx']])
  ny = np.array(data.data[data.keys['ny']])
  nz = np.array(data.data[data.keys['nz']])
  has_n = True

r = np.ones(len(x))  

Radii = vtk.vtkDoubleArray()
Radii.SetNumberOfComponents(1)
Radii.SetName('Radius')

if has_v:
  Velocities = vtk.vtkDoubleArray()
  Velocities.SetNumberOfComponents(3)
  Velocities.SetName("Velocity")

if has_n:
  Directors = vtk.vtkDoubleArray()
  Directors.SetNumberOfComponents(3)
  Directors.SetName("Directors")

for (xx,yy,zz,rr) in zip(x,y,z,r):
  Points.InsertNextPoint(xx,yy,zz)
  Radii.InsertNextValue(rr)
  
if has_v:
  for (vvx,vvy,vvz) in zip(vx,vy,vz):
    Velocities.InsertNextTuple3(vvx,vvy,vvz)

if has_n:
  for (nnx,nny,nnz) in zip(nx,ny,nz):
    Directors.InsertNextTuple3(nnx,nny,nnz)

polydata = vtk.vtkPolyData()
polydata.SetPoints(Points)

polydata.GetPointData().AddArray(Radii)

if has_v:
  polydata.GetPointData().AddArray(Velocities)

if has_n:
  polydata.GetPointData().AddArray(Directors)

polydata.Modified()
writer = vtk.vtkXMLPolyDataWriter()
writer.SetFileName(args.output+'.vtp')
writer.SetInputData(polydata)
writer.SetDataModeToAscii()
writer.Write()


end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print
  
