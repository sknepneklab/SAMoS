# *****************************************************************************
# *
# *  This Python script is a part of tha analysis of the data published in 
# *  the paper: "Universal motion patterns in confluent cell monolayers"
# *  by Silke Henkes, Kaja Kostanjevec, J. Martin Collinson, Rastko Sknepnek, 
# *  and Eric Bertin, Jounral name, vol, page (2019).
# *
# *  Please refer to the document Computational_summary.pdf for a detailed
# *  description of the tasks performed by this script.
# * 
# *****************************************************************************

# Utility code for building random initial configuration on xy plane

from datetime import *
import numpy.random
from random import uniform 
from math import *
import argparse

# Base class for Python particle

class Particle:
  
  def __init__(self, idx, a0, tp=1, R=1.0, l=1.0):
    self.idx = idx
    self.tp = tp
    self.R = R
    self.r = [0.0,0.0,0.0]
    self.v = [0.0,0.0,0.0]
    self.f = [0.0,0.0,0.0]
    self.n = [0.0,0.0,0.0]
    self.l = l
    self.normal = [0.0,0.0,1.0]
    self.area = a0

class Circle:
  def __init__(self, R, N, v, poly,boundary,bdens,a0):
    self.a0 = a0
    self.radius = sqrt (a0/pi)
    self.R = R
    self.N = N
    self.poly = poly
    self.boundary=boundary
    self.Nbound = int(bdens*2.0*pi*R/(2.0*self.radius)) # Dense boundary wall, linear packing fraction 2.0 (for stiff glued boundaries)
    area = self.a0*numpy.random.uniform(1-0.5*self.poly,1+0.5*self.poly,self.N+self.Nbound)
    self.particles = [Particle(i,area[i],1,  sqrt (area[i]/pi)) for i in range(N)]
    self.boundparticles = [Particle(i,area[i],2,  sqrt (area[i]/pi)) for i in range(N,N+self.Nbound)]
    self.__generate_posinside()
    self.__generate_posbound()
    self.__generate_vel(v)
    self.__generate_director()
    
  def __generate_posinside(self):
    for i in range(self.N):
      rpart=4*self.R**2
      while rpart>(self.R-self.radius)**2:
        x = uniform(-self.R+self.radius,self.R-self.radius)
        y = uniform(-self.R+self.radius,self.R-self.radius)
        z = 0
        rpart=x**2+y**2
      self.particles[i].r = [x,y,z]
  
  def __generate_posbound(self):
	for i in range(self.Nbound):
	  # fractional position along circumference (angle)
	  sval = 1.0*i/self.Nbound*2*pi
	  self.boundparticles[i].r=[self.R*cos(sval),self.R*sin(sval),0]
	  # introduce a little noise there; no ratchets here
	  #self.boundparticles[i].R = uniform(1-0.5*self.poly,1+0.5*self.poly)
	  
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
    out.write('keys: id type radius  x   y   z   vx   vy   vz   nx   ny   nz  nvx   nvy   nvz   area  boundary \n')
    for p in self.particles:
      x, y, z = p.r
      vx, vy, vz = p.v
      nx, ny, nz = p.n
      nvx, nvy, nvz = p.normal
      out.write('%d  %d  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %d\n' % (p.idx,p.tp,p.R,x,y,z,vx,vy,vz,nx,ny,nz,nvx,nvy,nvz,p.area,0))
    for p in self.boundparticles:
      x, y, z = p.r
      vx, vy, vz = p.v
      nx, ny, nz = p.n
      out.write('%d  %d  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %d\n' % (p.idx,p.tp,p.R,x,y,z,vx,vy,vz,nx,ny,nz,nvx,nvy,nvz,p.area,1))
    
    # write boundary file. Careful, this has to be clockwise and close on itself on the last one
  def writeBoundary(self,outfile):
    gentime = datetime.now()
    out = open(outfile,'w')
    out.write('# Generated on : %s\n' % str(gentime))
    out.write('# Total of %s boundary particles\n' % str(self.Nbound))
    out.write('# label id1 id2\n')
    # first line: first and last
    out.write("%d %d %d\n" % (0,self.N,self.N+self.Nbound-1))
    for k in range(1,self.Nbound):
        out.write("%d %d %d\n" % (k,self.N+self.Nbound-k,self.N+self.Nbound-k-1))
    out.close()
	  
parser = argparse.ArgumentParser()
#parser.add_argument("-r", "--R", type=float, default=10.0, help="glued circle radius")
parser.add_argument("-n", "--N", type=int, default=1500, help="number of particles")
parser.add_argument("-p", "--poly", type=float, default=0.0, help="polydispersity fraction")
parser.add_argument("-f", "--phi",  type=float, default=1.0, help="packing fraction")
parser.add_argument("-a", "--area",  type=float, default=382, help="cell area")
parser.add_argument("-o", "--output", type=str, default='out.dat', help="output file")
parser.add_argument("-i", "--bfile", type=str, default='boundary.dat', help="boundary output file")
parser.add_argument("-v", "--vavr", type=float, default=0.1, help="average velocity")
parser.add_argument("-b", "--boundary", type=int, default=0, help="boundary n: 0 random 1 outward 2 inward")
parser.add_argument("-d", "--bdens", type=float, default=1.2, help="packing fraction at the boundary")
args = parser.parse_args()

print
print "\tSoft Active Matter on Surfaces"
print "\tBuilding a glued circle random configuration on the plane as cell initial configuration"
print 
#print "\tRadius : ", args.R
print "\tpolydispersity: ", args.poly
#N = int(round(args.R**2*args.phi/4.0)) # An individual particle has area pi R^2 = 4 pi
print "\tPacking fraction : ", args.phi
print "\tNumber of particles : ", args.N
print "\tArea of particles : ", args.area
R = sqrt(args.N*args.area/(args.phi*pi))
print "\tRadius : ", R
print "\tAverage velocity : ", args.vavr
print "\tBoundary type (0 random 1 outward 2 inward) : ", args.boundary
print "\tOutput file : ", args.output
print 

start = datetime.now()
p = Circle(R, args.N, args.vavr,args.poly,args.boundary,args.bdens,args.area)
p.write(args.output)
p.writeBoundary(args.bfile)

end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print