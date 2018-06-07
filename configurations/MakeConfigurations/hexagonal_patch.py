# * ***************************************************************************
# *
# *  Copyright (C) 2013-2016 University of Dundee
# *  All rights reserved. 
# *
# *  This file is part of SAMoS (Soft Active Matter on Surfaces) program.
# *
# *  SAMoS is free software; you can redistribute it and/or modify
# *  it under the terms of the GNU General Public License as published by
# *  the Free Software Foundation; either version 2 of the License, or
# *
# *  (at your option) any later version.
# *  SAMoS is distributed in the hope that it will be useful,
# *  but WITHOUT ANY WARRANTY; without even the implied warranty of
# *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# *  GNU General Public License for more details.
# *
# *  You should have received a copy of the GNU General Public License
# *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
# *
# * ****************************************************************************/

# Utility code for building a regular hexagonal patch

from datetime import *
from math import *
import numpy as np
import argparse
from random import uniform 

# Base class for Python particle

class Particle:
  
  def __init__(self, idx, tp=1, R=1.0, l=1.0):
    self.idx = idx
    self.tp = tp
    self.R = R
    self.r = [0.0,0.0,0.0]
    self.v = [0.0,0.0,0.0]
    self.f = [0.0,0.0,0.0]
    self.n = [0.0,0.0,0.0]
    self.nvx = [0.0,0.0,1.0]
    self.area = 0.0
    self.boundary = 0
    

class TriangularPatch:
  
  def __init__(self, n, a, v):
    self.n = n
    self.a = a
    self.__generate_pos()
    self.__generate_vel(v)
    self.__generate_director()
    
  def __generate_pos(self):
    self.particles = []
    self.boundary = []
    vert_fact = 0.5*self.a*np.sqrt(3)
    idx = 0
    for i in range(self.n):
      y = vert_fact*i
      for j in range(2*self.n-1-i):
        x = 0.5*(-(2*self.n-i-2)*self.a) + j*self.a
        p = Particle(idx)
        p.r = [x, y, 0.0]
        p.area = vert_fact*self.a
        self.particles.append(p)
        if (i == self.n-1) and (not idx in self.boundary): self.boundary.append(idx)
        if (j == 0 or j == (2*self.n-1-i-1)) and (not idx in self.boundary): self.boundary.append(idx)
        idx += 1
        if (y != 0.0):
          p = Particle(idx)
          p.r = [x, -y, 0.0]
          p.area = vert_fact*self.a
          self.particles.append(p)
          if (j == 0 or j == (2*self.n-1-i-1)) and (not idx in self.boundary): self.boundary.append(idx)
          if (i == self.n-1) and (not idx in self.boundary): self.boundary.append(idx)
          idx += 1

  def __generate_vel(self, vav=1.0):
    for p in self.particles:
      phi = uniform(0,2*pi)
      p.v = [vav*cos(phi), vav*sin(phi), 0.0]
  
  def __generate_director(self):
    for p in self.particles:
      phi = uniform(0,2*pi)
      p.n = [cos(phi), sin(phi), 0.0]
  
  def scale(self, sx, sy):
    for p in self.particles:
      p.r[0] *= sx 
      p.r[1] *= sy 

  def write(self,outfile):
    gentime = datetime.now()
    out = open(outfile,'w')
    out.write('# Total of %d particles\n' % len(self.particles))
    out.write('# Generated on : %s\n' % str(gentime))
    out.write('keys: id  type radius  x   y   z   vx   vy   vz   nx   ny   nz nvx nvy nvz area boundary\n')
    for p in self.particles:
      x, y, z = p.r
      vx, vy, vz = p.v
      nx, ny, nz = p.n
      nvx, nvy, nvz = p.nvx
      if p.idx in self.boundary: bnd = 1
      else: bnd = 0
      out.write('%d  %d  %f %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f %d\n' % (p.idx,p.tp,p.R,x,y,z,vx,vy,vz,nx,ny,nz,nvx,nvy,nvz,p.area,bnd))
    out.close()

  def write_boundary(self, outfile):
    angles = np.zeros(len(self.boundary))
    for i in range(angles.size):
      p = self.particles[self.boundary[i]]
      angles[i] = np.arctan2(p.r[1],p.r[0])
    with open(outfile,'w') as out:
      idxs = np.argsort(angles)
      for i in range(idxs.size-1):
        out.write('{:3d} {:5d} {:5d}\n'.format(i, self.boundary[idxs[i]], self.boundary[idxs[i+1]]))
      out.write('{:3d} {:5d} {:5d}\n'.format(idxs.size-1, self.boundary[idxs[-1]], self.boundary[idxs[0]]))

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--n", type=int, default=3, help="number of layers")
parser.add_argument("-a", "--a",  type=float, default=1.0, help="lattice constant")
parser.add_argument("-o", "--output", type=str, default='out.dat', help="output file")
parser.add_argument("-v", "--vavr", type=float, default=1.0, help="average velocity")
parser.add_argument("-x", "--sx", type=float, default=1.0, help="scale in x direction")
parser.add_argument("-y", "--sy", type=float, default=1.0, help="scale in y direction")
parser.add_argument("-b", "--boundary", type=str, default=None, help="boundary particles")
args = parser.parse_args()

print
print "\tSoft Active Matter on Surfaces (SAMoS)"
print "\tBuilding a regular hexagonal patch on a plane"
print 
print "\tRastko Sknepnek"
print "\tUniversity of Dundee"
print "\t(c) 2018"
print "\t----------------------------------------------"
print 
print "\tLattice constant : ", args.a
print "\tNumber of layers : ", args.n
print "\tAverage velocity : ", args.vavr
print "\tOutput file : ", args.output
if args.boundary != None:
  print "\tBoundary file : ", args.boundary
print 

start = datetime.now()

t = TriangularPatch(args.n, args.a, args.vavr)
t.scale(args.sx, args.sy)
t.write(args.output)
if args.boundary != None:
  t.write_boundary(args.boundary)

end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print
