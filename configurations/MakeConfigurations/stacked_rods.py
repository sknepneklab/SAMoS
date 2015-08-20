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

class PlaneRods:
  
  def __init__(self, Lx, Ly, Nx,Ny, v, random_orient = True):
    self.L = (Lx,Ly)
    self.N = (Nx,Ny)
    self.__generate_lattice()
    self.__generate_vel(v)
    self.__generate_director(random_orient)
    

  def __generate_lattice(self):
    self.particles = [Particle(i) for i in range(self.N[0]*self.N[1])]
    idx = 0
    dx=self.L[0]/self.N[0]
    dy=self.L[1]/self.N[1]
    for i in range(self.N[0]):
      x = -0.5*self.L[0] + (0.5+i)*dx
      for j in range(self.N[1]):
        y = -0.5*self.L[1] + (0.5+j)*dy
        self.particles[idx].r = [x,y,0]
        idx += 1

  def __generate_vel(self,vav=1.0):
    for p in self.particles:
      phi = uniform(0,2*pi)
      p.v = [vav*cos(phi),vav*sin(phi),0.0]
  
  def __generate_director(self, random_orinet):
    for p in self.particles:
      if random_orient:
        if uniform(-0.5,0.5) < 0.0:
          phi = 0.0 
        else:
          phi = pi
      else:
        phi = 0.0
        
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
    ntot=self.N[0]*self.N[1]
    out.write('# Total of %d particles\n' % ntot )
    out.write('# System sizes lx = ' + str(self.L[0]) + ' ly = ' + str(self.L[1]) + '\n')
    out.write('# Generated on : %s\n' % str(gentime))
    out.write('# id  type radius  x   y   z   vx   vy   vz   nx   ny   nz  omega l\n')
    for p in self.particles:
      x, y, z = p.r
      vx, vy, vz = p.v
      nx, ny, nz = p.n
      out.write('%d  %d  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n' % (p.idx,p.tp,p.R,x,y,z,vx,vy,vz,nx,ny,nz,p.omega,p.l))
    out.close()
    
parser = argparse.ArgumentParser()
parser.add_argument("-L", "--laverage", type=float, default=10.0, help="system size")
parser.add_argument("-f", "--phi",  type=float, default=0.5, help="packing fraction")
parser.add_argument("-r", "--sigma",  type=float, default=1.0, help="rod cap radius / half width")
parser.add_argument("-l", "--length",  type=float, default=2.0, help="rod length (minus caps)")
parser.add_argument("-o", "--output", type=str, default='out.dat', help="output file")
parser.add_argument("-v", "--vavr", type=float, default=1.0, help="average velocity")
parser.add_argument("-O", "--orient", type=bool, default=False, help="random orientation of rods")
args = parser.parse_args()



print
print "\tActive Particles on Curved Spaces (APCS)"
print "\tBuilding of a random flat configuration (xy plane)"
print 
print "\tRastko Sknepnek"
print "\tUniversity of Dundee"
print "\t(c) 2013, 2014, 2015"
print "\t----------------------------------------------"
print 


start = datetime.now()


N = int(args.phi*args.laverage**2/(2*args.length*args.sigma+pi*args.sigma**2))
print N
aspect = (args.length+2*args.sigma)/(2.0*args.sigma)
print aspect
nx=int(sqrt(N)/aspect)
ny=int(sqrt(N)*aspect)
lx=aspect*args.laverage*nx/sqrt(N)
ly=args.laverage*ny/(aspect*sqrt(N))
print nx
print ny
print lx
print ly


random_orient = args.orient

p = PlaneRods(lx, ly, nx,ny, args.vavr, random_orient=random_orient)

radii = []
types = []
lens = []
for i in xrange(len(p.particles)): 
  radii.append(args.sigma)
  types.append(1)
  lens.append(args.length)
p.set_radius(radii)
p.set_type(types)
p.set_lens(lens)
p.write(args.output)

end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print
