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
parser.add_argument("-E", "--end_sample", type=int, default=None, help="last sample")
parser.add_argument("-S", "--step", type=int, default=1, help="step between samples")
parser.add_argument("-C", "--contact", type=str, default=None, help="contact network data file")
parser.add_argument("-e", "--exclude", type=float, default=None, help="exclude all contact line that are longer than this value")
parser.add_argument("--connected", action='store_true', default=False, help="Include Delaunay triangulation data")
parser.add_argument("--nematic", action='store_true', default=False, help="Shift n vectors such that particle is in the middle of director.")
parser.add_argument("-l", "--length", type=float, default=1.0, help="rod length")
parser.add_argument("-b", "--bonds", type=str, default=None, help="bond file")
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
if args.contact != None:
  print "\tContact network data file : ", args.contact
if args.exclude != None:
  print "\tExclude all contact lines that are longer than : ", args.exclude 
print

start = datetime.now()

if args.end_sample == None:
  files = sorted(glob(args.input+'*.dat'))[args.skip::args.step]
  if len(files) == 0:
    files = sorted(glob(args.input+'*.dat.gz'))[args.skip::args.step]
else:
  files = sorted(glob(args.input+'*.dat'))[args.skip:args.end_sample:args.step]
  if len(files) == 0:
    files = sorted(glob(args.input+'*.dat.gz'))[args.skip:args.end_sample:args.step]

if args.contact != None:
  cont_files = sorted(glob(args.contact+'*.con'))[args.skip:]  
  if len(files) != len(cont_files):
    print "There has to be same number of data and contact files."
    sys.exit(1)

# read bonds
bonds = []
if args.bonds != None:
  with open(args.bonds,'r') as bond_file:
    lines = bond_file.readlines()
    #print lines
    #lines = lines.split('\n')
    lines = map(lambda x: x.strip(), lines)
    for line in lines:
      b = line.split()
      bonds.append((int(b[2]),int(b[3])))
    
    
u=0
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
  
  Lx, Ly, Lz = np.max(x) - np.min(x), np.max(y) - np.min(y), np.max(z) - np.min(z)

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

  if (data.keys.has_key('radius')):
    r = np.array(data.data[data.keys['radius']])
  else:
    r = np.ones(len(x))  
  
  if (data.keys.has_key('type')):
    tp = np.array(data.data[data.keys['type']])
  else:
    tp = np.ones(len(x))
    
  if (data.keys.has_key('flag')):
    flag = np.array(data.data[data.keys['flag']])
  else:
    flag = np.arange(len(x))
  

  Radii = vtk.vtkDoubleArray()
  Radii.SetNumberOfComponents(1)
  Radii.SetName('Radius')
  
  Types = vtk.vtkDoubleArray()
  Types.SetNumberOfComponents(1)
  Types.SetName('Type')
  
  Flags = vtk.vtkDoubleArray()
  Flags.SetNumberOfComponents(1)
  Flags.SetName('Flag')

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

  for (xx,yy,zz,rr,t,f) in zip(x,y,z,r,tp,flag):
    Points.InsertNextPoint(xx,yy,zz)
    Radii.InsertNextValue(rr)
    Types.InsertNextValue(t)
    Flags.InsertNextValue(f)
    
  if has_v:
    for (vvx,vvy,vvz) in zip(vx,vy,vz):
      Velocities.InsertNextTuple3(vvx,vvy,vvz)

  if has_n:
    rod_len = args.length
    for (nnx,nny,nnz) in zip(nx,ny,nz):
      if args.nematic:
        Directors.InsertNextTuple3(0.5*rod_len*nnx,0.5*rod_len*nny,0.5*rod_len*nnz)
        NDirectors.InsertNextTuple3(-0.5*rod_len*nnx,-0.5*rod_len*nny,-0.5*rod_len*nnz)
      else:
        Directors.InsertNextTuple3(rod_len*nnx,rod_len*nny,rod_len*nnz)        

  if args.contact != None:
    if args.connected:
      print "Error! You cannot use --connected flag with the contact network data"
      sys.exit(1)
    Lines = vtk.vtkCellArray()
    Line = vtk.vtkLine()
    contact = open(cont_files[u],'r')
    con_lines = contact.readlines()
    con_lines = map(lambda x: x.strip().split(),con_lines)
    edges = []
    for line in con_lines:
      if not line[0] == '#':
        i, j = int(line[1]), int(line[2])
        if args.exclude != None:
          dr = np.sqrt((x[i]-x[j])**2 + (y[i]-y[j])**2 + (z[i]-z[j])**2)
          if dr < args.exclude:
            edges.append((i,j))
        else:
          edges.append((i,j))
    contact.close()
    for (i,j) in edges:
      Line.GetPointIds().SetId(0,i)
      Line.GetPointIds().SetId(1,j)
      Lines.InsertNextCell(Line)

  if args.bonds != None:
    Lines = vtk.vtkCellArray()
    Line = vtk.vtkLine()
    Lengths = vtk.vtkDoubleArray()
    Lengths.SetNumberOfComponents(1)
    Lengths.SetName('BondLength')
    for (i,j) in bonds:
      dr = np.sqrt((x[i]-x[j])**2 + (y[i]-y[j])**2 + (z[i]-z[j])**2)
      if dr < 0.5*max([Lx,Ly,Lz]):
        Line.GetPointIds().SetId(0,i)
        Line.GetPointIds().SetId(1,j)
        Lines.InsertNextCell(Line)
        Lengths.InsertNextValue(dr)

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

  if args.contact != None:
    polydata.SetLines(Lines)

  if args.bonds != None:
    polydata.SetLines(Lines)
    polydata.GetCellData().AddArray(Lengths)

  polydata.GetPointData().AddArray(Radii)
  polydata.GetPointData().AddArray(Types)
  polydata.GetPointData().AddArray(Flags)

  if has_v:
    polydata.GetPointData().AddArray(Velocities)
    
  if has_n:
    polydata.GetPointData().AddArray(Directors)
    if args.nematic:
      polydata.GetPointData().AddArray(NDirectors)

  polydata.Modified()
  writer = vtk.vtkXMLPolyDataWriter()
  #outname = '.'.join(f.split('.')[:-1])
  #print outname
  outname='frame_%09d' % u
  u+=1
  writer.SetFileName(args.output+'/'+outname+'.vtp')
  if vtk.VTK_MAJOR_VERSION <= 5:
    writer.SetInput(polydata)
  else:
    writer.SetInputData(polydata)
  writer.SetDataModeToAscii()
  writer.Write()
  
end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print
  
