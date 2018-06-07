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

# Utility code for building random initial configuration on xy plane

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
    self.l = l
    self.omega = 0.0


class TriangularPlane:
  
  def __init__(self, length, nside, N, v):
    self.L = (length,length)
    self.N = N
    self.nside=nside
    self.__generate_pos()
    self.__generate_vel(v)
    self.__generate_director()
    
  def __generate_pos(self):
    self.particles = [Particle(i) for i in range(self.N)]
    # Triangular lattice now. Is nside x nside, first put it on a tilted polygon compatible with the lattice
    # Apply periodic boundary conditions if necessary
    dx = self.L[0]/self.nside
    dy = self.L[1]/self.nside*np.sin(np.pi/3)
    for i in range(self.nside):
      for j in range(self.nside):
	x = i*dx + 0.5*j*dx
	if x > self.L[0]:
          x = x - self.L[0]
	y = j*dy
	z = 0
	k = nside*j+i
	self.particles[k].r = [x,y,z]


  def __generate_vel(self,vav=1.0):
    for p in self.particles:
      phi = uniform(0,2*pi)
      p.v = [vav*cos(phi),vav*sin(phi),0.0]
  
  def __generate_director(self):
    for p in self.particles:
      phi = uniform(0,2*pi)
      p.n = [cos(phi),sin(phi),0.0]
  
  def write(self,outfile):
    gentime = datetime.now()
    out = open(outfile,'w')
    out.write('# Total of %d particles\n' % self.N)
    out.write('# Generated on : %s\n' % str(gentime))
    out.write('# id  type radius  x   y   z   vx   vy   vz   nx   ny   nz \n')
    for p in self.particles:
      x, y, z = p.r
      vx, vy, vz = p.v
      nx, ny, nz = p.n
      out.write('%d  %d  %f %f  %f  %f  %f  %f  %f  %f  %f  %f \n' % (p.idx,p.tp,p.R,x,y,z,vx,vy,vz,nx,ny,nz))
    out.close()
    
parser = argparse.ArgumentParser()
parser.add_argument("-L", "--length", type=float, default=10.0, help="box length")
parser.add_argument("-f", "--phi",  type=float, default=0.5, help="packing fraction")
parser.add_argument("-o", "--output", type=str, default='out.dat', help="output file")
parser.add_argument("-v", "--vavr", type=float, default=1.0, help="average velocity")
args = parser.parse_args()

print
print "\tActive Particles on Curved Spaces (APCS)"
print "\tBuilding a triangular configuration (xy plane)"
print 
print "\tSilke Henkes"
print "\tUniversity of Aberdeen"
print "\t(c) 2016"
print "\t----------------------------------------------"
print 
print "\tL : ", args.length
print "\tPacking fraction : ", args.phi
print "\tAverage velocity : ", args.vavr
print "\tOutput file : ", args.output
print 

start = datetime.now()

# Create number of particles. phi = N pi r^2/L^2, so N = phi L^2/ pi r_0^2
# Round to the nearest square
r0=1.0
N0 = args.phi*args.length**2/(np.pi*r0**2)
print "Raw number of particles " + str(N0)
nside = int(round(np.sqrt(N0)))
print "Number of particles a side " + str(nside)
N=int(nside**2)
phireal = N*np.pi*r0**2/args.length**2
print "Number of particles " + str(N) + " and real phi " + str(phireal)

p = TriangularPlane(args.length,nside,N,args.vavr)
p.write(args.output)

end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print
