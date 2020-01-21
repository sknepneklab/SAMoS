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
from random import uniform 
from math import *
import argparse


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
    

class Plane:
  
  def __init__(self, Lx, Ly, N, poly, ratio, v):
    self.L = (Lx,Ly)
    self.N = N
    self.ratio=ratio
    self.particles = [Particle(i) for i in range(N)]
    self.__generate_pos()
    self.__generate_vel(v)
    self.__generate_director()
    self.set_radius(poly)
    self.set_type()
      
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
  
  def set_radius(self,poly):
    for i in range(self.N):
      self.particles[i].R = uniform(1-0.5*poly,1+0.5*poly)
  
  def set_type(self):
	ntracer=int(self.ratio*self.N)
	for i in range(self.N-ntracer):
	  self.particles[i].tp = 1
	for i in range(self.N-ntracer,self.N):
	  self.particles[i].tp = 2
  
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
    
parser = argparse.ArgumentParser()
parser.add_argument("-x", "--lx", type=float, default=10.0, help="box length in x direction")
parser.add_argument("-y", "--ly", type=float, default=10.0, help="box length in y direction")
parser.add_argument("-f", "--phi",  type=float, default=0.5, help="packing fraction")
parser.add_argument("-o", "--output", type=str, default='out.dat', help="output file")
parser.add_argument("-v", "--vavr", type=float, default=1.0, help="average velocity")
parser.add_argument("-p", "--poly", type=float, default=0.0, help="polydispersity fraction")
parser.add_argument("-r", "--ratio", type=float, default=0.0, help="fraction of tracer particles")
args = parser.parse_args()

#N = int(round(1.0/pi*args.lx*args.ly*args.phi/(args.eta*args.a1**2+(1-args.eta)*args.a2**2)))
N = int(round(1.0/pi*args.lx*args.ly*args.phi))

print
print "\tSoft Active Matter on Surfaces (SAMoS)"
print "\tBuilding of a random flat configuration (xy plane)"
print "\t----------------------------------------------"
print 
print "\tLx : ", args.lx
print "\tLy : ", args.ly
print "\tPacking fraction : ", args.phi
print "\tNumber of particles : ", N
print "\tAverage velocity : ", args.vavr
print "\tOutput file : ", args.output
print "\tPolydispersity :", args.poly
print "\tFraction of tracer particles :", args.ratio
print 

start = datetime.now()

p = Plane(args.lx, args.ly, N, args.poly, args.ratio, args.vavr)
p.write(args.output)

end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print