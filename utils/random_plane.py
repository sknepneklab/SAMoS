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
#    (c) 2013
#    
#    This program cannot be used, copied, or modified without
#    explicit permission of the author.
#  
# ################################################################

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
      phi = 0.0 #uniform(0,2*pi)
      p.n = [cos(phi),sin(phi),0.0]
  
  def set_radius(self,radii):
    for i in xrange(len(radii)):
      self.particles[i].R = radii[i]
  
  def set_type(self,types):
    for i in xrange(len(types)):
      self.particles[i].tp = types[i]
  
  def set_lens(self,lens):
    for i in xrange(len(lens)):
      self.particles[i].l = lens[i]
  
  def write(self,outfile):
    gentime = datetime.now()
    out = open(outfile,'w')
    out.write('# Total of %d particles\n' % self.N)
    out.write('# Generated on : %s\n' % str(gentime))
    out.write('# id  type radius  x   y   z   vx   vy   vz   nx   ny   nz  omega l\n')
    for p in self.particles:
      x, y, z = p.r
      vx, vy, vz = p.v
      nx, ny, nz = p.n
      out.write('%d  %d  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n' % (p.idx,p.tp,p.R,x,y,z,vx,vy,vz,nx,ny,nz,p.omega,p.l))
    out.close()
    
parser = argparse.ArgumentParser()
parser.add_argument("-x", "--lx", type=float, default=10.0, help="box length in x direction")
parser.add_argument("-y", "--ly", type=float, default=10.0, help="box length in y direction")
parser.add_argument("-f", "--phi",  type=float, default=0.5, help="packing fraction")
parser.add_argument("-o", "--output", type=str, default='out.dat', help="output file")
parser.add_argument("-v", "--vavr", type=float, default=1.0, help="average velocity")
parser.add_argument("--a1",  type=float, default=1.0, help="radius of particles of type 1")
parser.add_argument("--a2",  type=float, default=0.5, help="radius of particles of type 2")
parser.add_argument("--eta",  type=float, default=0.5, help="fraction of particles of type 1")
parser.add_argument("--l1",  type=float, default=2.0, help="length of rod of type 1")
parser.add_argument("--l2",  type=float, default=1.0, help="length of rod of type 2")
args = parser.parse_args()

N = int(round(1.0/pi*args.lx*args.ly*args.phi/(args.eta*args.a1**2+(1-args.eta)*args.a2**2)))

print
print "\tActive Particles on Curved Spaces (APCS)"
print "\tBuilding of a random flat configuration (xy plane)"
print 
print "\tRastko Sknepnek"
print "\tUniversity of Dundee"
print "\t(c) 2013"
print "\t----------------------------------------------"
print 
print "\tLx : ", args.lx
print "\tLy : ", args.ly
print "\tPacking fraction : ", args.phi
print "\tNumber of particles : ", N
print "\tAverage velocity : ", args.vavr
print "\tOutput file : ", args.output
print "\tFraction of particles of type 1 : ", args.eta
print "\tRadius of particles of type 1 : ", args.a1
print "\tRadius of particles of type 2 : ", args.a2
print 

start = datetime.now()

p = Plane(args.lx, args.ly, N, args.vavr)
radii = []
types = []
lens = []
for i in xrange(N):
  if i < args.eta*N: 
    radii.append(args.a1)
    types.append(1)
    lens.append(args.l1)
  else: 
    radii.append(args.a2)
    types.append(2)
    lens.append(args.l2)
p.set_radius(radii)
p.set_type(types)
p.set_lens(lens)
p.write(args.output)

end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print