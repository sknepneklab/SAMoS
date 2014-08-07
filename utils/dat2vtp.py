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
import numpy.linalg as lin
from datetime import *
from glob import glob
from scipy.spatial import ConvexHull
import math as m
import vtk


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="input file (base name)")
parser.add_argument("-o", "--output", type=str, help="output directory")
parser.add_argument("-s", "--skip", type=int, default=0, help="skip this many samples")
parser.add_argument("--connected", action='store_true', default=False, help="Include Delaunay triangulation data")
parser.add_argument("--nematic", action='store_true', default=False, help="Shift n vectors such that particle is in the middle of director.")
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
print "\tSkip frames : ", args.skip
print

start = datetime.now()

files = sorted(glob(args.input+'*.dat'))[args.skip:]

for f in files:
  print "Processing file : ", f

  Points = vtk.vtkPoints()


  has_v = False
  has_n = False

  data = ReadData(f)

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
    # Negative directors to mimic neamtic (silly but should work)
    if args.nematic:
      NDirectors = vtk.vtkDoubleArray()
      NDirectors.SetNumberOfComponents(3)
      NDirectors.SetName("NDirectors")

  for (xx,yy,zz,rr) in zip(x,y,z,r):
    Points.InsertNextPoint(xx,yy,zz)
    Radii.InsertNextValue(rr)
    
  if has_v:
    for (vvx,vvy,vvz) in zip(vx,vy,vz):
      Velocities.InsertNextTuple3(vvx,vvy,vvz)

  if has_n:
    for (nnx,nny,nnz) in zip(nx,ny,nz):
      if args.nematic:
        Directors.InsertNextTuple3(0.5*nnx,0.5*nny,0.5*nnz)
        NDirectors.InsertNextTuple3(-0.5*nnx,-0.5*nny,-0.5*nnz)
      else:
        Directors.InsertNextTuple3(nnx,nny,nnz)        

  if args.connected:
    Lines = vtk.vtkCellArray()
    Line = vtk.vtkLine()
    Lengths = vtk.vtkDoubleArray()
    Lengths.SetNumberOfComponents(1)
    Lengths.SetName('Length')
    NNeighs = vtk.vtkDoubleArray()
    NNeighs.SetNumberOfComponents(1)
    NNeighs.SetName('NNeigh')
    Areas = vtk.vtkDoubleArray()
    Areas.SetNumberOfComponents(1)
    Areas.SetName('Area')
    Faces = vtk.vtkCellArray()
    Polygon = vtk.vtkPolygon()
    points = np.column_stack((x,y,z)) 
    hull = ConvexHull(points)
    edges = []
    nneighs = [0 for i in xrange(len(x))]
    for h in hull.simplices:
      i, j, k = h
      if not sorted([i,j]) in edges: edges.append(sorted([i,j]))
      if not sorted([i,k]) in edges: edges.append(sorted([i,k]))
      if not sorted([j,k]) in edges: edges.append(sorted([j,k]))
      #a1 = points[j]-points[i]
      #a2 = points[k]-points[i]
      #area = 0.5*lin.norm(np.cross(a1,a2))
      #Areas.InsertNextValue(area)
      #Polygon.GetPointIds().SetNumberOfIds(3)
      #Polygon.GetPointIds().SetId(0, i)
      #Polygon.GetPointIds().SetId(1, j)
      #Polygon.GetPointIds().SetId(2, k)
      #Faces.InsertNextCell(Polygon)
    for (i,j) in edges:
      Line.GetPointIds().SetId(0,i)
      Line.GetPointIds().SetId(1,j)
      Lines.InsertNextCell(Line)
      nneighs[i] += 1
      nneighs[j] += 1
      dx, dy, dz = x[i]-x[j], y[i]-y[j], z[i]-z[j]
      Lengths.InsertNextValue(m.sqrt(dx*dx + dy*dy + dz*dz))
    for n in nneighs:
      NNeighs.InsertNextValue(n)

  polydata = vtk.vtkPolyData()
  polydata.SetPoints(Points)
  if args.connected:
    polydata.GetPointData().AddArray(NNeighs)
    polydata.SetLines(Lines)
    polydata.GetCellData().AddArray(Lengths)
    #polydata.SetPolys(Faces)
    #polydata.GetCellData().AddArray(Areas)

  polydata.GetPointData().AddArray(Radii)

  if has_v:
    polydata.GetPointData().AddArray(Velocities)

  if has_n:
    polydata.GetPointData().AddArray(Directors)
    if args.nematic:
      polydata.GetPointData().AddArray(NDirectors)

  polydata.Modified()
  writer = vtk.vtkXMLPolyDataWriter()
  outname = '.'.join(f.split('.')[:-1])
  writer.SetFileName(args.output+'/'+outname+'.vtp')
  writer.SetInputData(polydata)
  writer.SetDataModeToAscii()
  writer.Write()
  
end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print
  
