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

# Utility code for computing trajectory of a given particle

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

if len(files) == 0:
  files = sorted(glob(args.input+'*.dat.gz'))[args.skip:]

out = open(args.output+'.dat','w')
pid = args.id

Points = vtk.vtkPoints()
Lines = vtk.vtkCellArray()


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

Velocities = vtk.vtkDoubleArray()
Velocities.SetNumberOfComponents(3)
Velocities.SetName("Velocity")

Directors = vtk.vtkDoubleArray()
Directors.SetNumberOfComponents(3)
Directors.SetName("Directors")

for (x,y,z,vx,vy,vz,nx,ny,nz) in zip(xx,yy,zz,vvx,vvy,vvz,nnx,nny,nnz):
  Points.InsertNextPoint(x,y,z)
  Velocities.InsertNextTuple3(vx,vy,vz)
  Directors.InsertNextTuple3(nx,ny,nz)

Line = vtk.vtkLine()
for i in range(len(xx)-1):
  Line.GetPointIds().SetId(0,i)
  Line.GetPointIds().SetId(1,i+1)
  Lines.InsertNextCell(Line)

polydata = vtk.vtkPolyData()
polydata.SetPoints(Points)
polydata.SetLines(Lines)

polydata.GetPointData().AddArray(Velocities)
polydata.GetPointData().AddArray(Directors)

polydata.Modified()
writer = vtk.vtkXMLPolyDataWriter()
writer.SetFileName(args.output+'.vtp')
writer.SetInputData(polydata)
writer.SetDataModeToAscii()
writer.Write()

out.close()

#vtk_writer = VTK_XML_Serial_Unstructured()

#vtk_writer.snapshot(args.output+'.vtu',xx,yy,zz,vx=vvx,vy=vvy,vz=vvz,nx=nnx,ny=nny,nz=nnz)


end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print


