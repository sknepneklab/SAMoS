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
from random import uniform 
from math import *
import argparse


from particle import *

class Plane:
  
  def __init__(self, Lx, Ly, N, v):
    self.L = (Lx,Ly)
    self.N = N
    self.particles = [Particle(i) for i in range(N)]
    self.__generate_pos()
    self.__generate_vel(v)
    self.__generate_director()
    
  def __generate_pos(self):
    for i in range(self.N):
      x = uniform(-0.5*self.L[0],0.5*self.L[0])
      y = uniform(-0.5*self.L[1],0.5*self.L[1])
      z = 0
      self.particles[i].r = [x,y,z]

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
    out.write('# id  type radius  x   y   z   vx   vy   vz   nx   ny   nz  omega\n')
    for p in self.particles:
      x, y, z = p.r
      vx, vy, vz = p.v
      nx, ny, nz = p.n
      out.write('%d  %d  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n' % (p.idx,p.tp,p.R,x,y,z,vx,vy,vz,nx,ny,nz,p.omega))
    out.close()
    
class Circle:
  def __init__(self, R, N, v, poly,boundary):
    self.R = R
    self.N = N
    self.poly = poly
    self.boundary=boundary
    self.Nbound = int(1.5*2.0*pi*R/2.0) # Dense boundary wall, linear packing fraction 1.5
    self.particles = [Particle(i,1) for i in range(N)]
    self.boundparticles = [Particle(i,2) for i in range(N,N+self.Nbound)]
    self.__generate_posinside()
    self.__generate_posbound()
    self.__generate_vel(v)
    self.__generate_director()
    
  def __generate_posinside(self):
    for i in range(self.N):
      rpart=4*self.R**2
      while rpart>(self.R-1.0)**2:
        x = uniform(-self.R+1.0,self.R-1.0)
        y = uniform(-self.R+1.0,self.R-1.0)
        z = 0
        rpart=x**2+y**2
      self.particles[i].r = [x,y,z]
      self.particles[i].R = uniform(1-0.5*self.poly,1+0.5*self.poly)
  
  def __generate_posbound(self):
	for i in range(self.Nbound):
	  # fractional position along circumference (angle)
	  sval = 1.0*i/self.Nbound*2*pi
	  self.boundparticles[i].r=[self.R*cos(sval),self.R*sin(sval),0]
	  # introduce a little noise there; no ratchets here
	  self.boundparticles[i].R = uniform(1-0.5*self.poly,1+0.5*self.poly)
	  
  def __generate_vel(self,vav=1.0):
    for p in self.particles:
      phi = uniform(0,2*pi)
      p.v = [vav*cos(phi),vav*sin(phi),0.0]
  
  def __generate_director(self):
    for p in self.particles:
      phi = uniform(0,2*pi)
      p.n = [cos(phi),sin(phi),0.0]
    # This is actually intriguing: make them point inwards, outwards or random
    for i in range(self.Nbound):
	  if self.boundary==2:
	    sval = 1.0*i/self.Nbound*2*pi+pi
	  elif self.boundary==1:
	    sval = 1.0*i/self.Nbound*2*pi
	  else:
		sval = uniform(0,2*pi)
	  self.boundparticles[i].n=[cos(sval),sin(sval),0]
	  
  def write(self,outfile):
    gentime = datetime.now()
    out = open(outfile,'w')
    out.write('# Total of %s particles\n' % str(self.N+self.Nbound))
    out.write('# Generated on : %s\n' % str(gentime))
    out.write('# id  type radius  x   y   z   vx   vy   vz   nx   ny   nz  omega\n')
    for p in self.particles:
      x, y, z = p.r
      vx, vy, vz = p.v
      nx, ny, nz = p.n
      out.write('%d  %d  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n' % (p.idx,p.tp,p.R,x,y,z,vx,vy,vz,nx,ny,nz,p.omega))
    for p in self.boundparticles:
      x, y, z = p.r
      vx, vy, vz = p.v
      nx, ny, nz = p.n
      out.write('%d  %d  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n' % (p.idx,p.tp,p.R,x,y,z,vx,vy,vz,nx,ny,nz,p.omega))
    out.close()
	  
parser = argparse.ArgumentParser()
#parser.add_argument("-r", "--R", type=float, default=10.0, help="glued circle radius")
parser.add_argument("-n", "--N", type=int, default=1000, help="number of particles")
parser.add_argument("-p", "--poly", type=float, default=0.0, help="polydispersity fraction")
parser.add_argument("-f", "--phi",  type=float, default=1.0, help="packing fraction")
parser.add_argument("-o", "--output", type=str, default='out.dat', help="output file")
parser.add_argument("-v", "--vavr", type=float, default=0.1, help="average velocity")
parser.add_argument("-b", "--boundary", type=int, default=0, help="boundary n: 0 random 1 outward 2 inward")
args = parser.parse_args()

print
print "\tActive Particles on Curved Spaces (APCS)"
print "\tBuilding a glued circle random configuration on the plane"
print 
print "\tRastko Sknepnek"
print "\tUniversity of Dundee"
print "\t(c) 2013"
print 
print "\tSilke Henkes"
print "\tUniversity of Aberdeen"
print "\t(c) 2014"
print "\t----------------------------------------------"
print 
#print "\tRadius : ", args.R
print "\tpolydispersity: ", args.poly
#N = int(round(args.R**2*args.phi/4.0)) # An individual particle has area pi R^2 = 4 pi
print "\tPacking fraction : ", args.phi
print "\tNumber of particles : ", args.N
R = sqrt(args.N/args.phi)
print "\tRadius : ", R
print "\tAverage velocity : ", args.vavr
print "\tBoundary type (0 random 1 outward 2 inward) : ", args.boundary
print "\tOutput file : ", args.output
print 

start = datetime.now()
#def __init__(self, R, N, v, poly):
p = Circle(R, args.N, args.vavr,args.poly,args.boundary)
p.write(args.output)

end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print