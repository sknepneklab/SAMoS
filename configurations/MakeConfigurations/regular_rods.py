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

# Utility code for building regular rods initial configuration on xy plane

from datetime import *
from random import uniform, randint 
from math import *
import argparse


from particle import *

class Plane:
  
  def __init__(self, Lx, Ly, N, lx, ly, sigma, l):
    self.L = (Lx,Ly)
    self.N = N
    self.lx = lx
    self.ly = ly
    self.sigma = sigma
    self.l = l
    self.__generate()
        
  def __generate(self):
    self.particles = []
    i = 0
    n = 0
    add_n = True
    while add_n:
      x = -0.5*self.L[0] + (n+0.5)*self.lx
      if x > 0.5*(self.L[0]-self.lx)+1e-3:
        add_n = False
      else:
        m = 0
        add_m = True
        while add_m:
          y = -0.5*self.L[1] + (m+0.5)*self.ly
          if y > 0.5*(self.L[1]-self.ly)+1e-3:
            add_m = False
          else:
            self.particles.append(Particle(i))
            self.particles[i].r = [x,y,0.0]
            self.particles[i].n = [2.0*(randint(0,1)-0.5),0.0,0.0]
            self.particles[i].v = [0.0,0.0,0.0]
            self.particles[i].R = 0.5*self.sigma
            self.particles[i].l = self.l
            self.particles[i].omega = 0.0
            i += 1
            m += 1
        n += 1
  
  def write(self,outfile):
    gentime = datetime.now()
    out = open(outfile,'w')
    out.write('# Total of %d particles\n' % len(self.particles))
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
parser.add_argument("-a", "--radius", type=float, default=0.5, help="rod radius")
parser.add_argument("-l", "--length", type=float, default=2.0, help="rod length")
args = parser.parse_args()

area = args.lx*args.ly
sigma = 2.0*args.radius
Arod = sigma*(args.length+0.25*pi*sigma)
N = int(round(area*args.phi/Arod))
p = (args.length+sigma)/sigma
lx = sqrt(Arod*p/args.phi)
ly = sqrt(Arod/(p*args.phi))

print p, lx, ly

print
print "\tActive Particles on Curved Spaces (APCS)"
print "\tBuilding of a random flat configuration (xy plane)"
print 
print "\tRastko Sknepnek"
print "\tUniversity of Dundee"
print "\t(c) 2013, 2014, 2015"
print "\t----------------------------------------------"
print 
print "\tLx : ", args.lx
print "\tLy : ", args.ly
print "\tPacking fraction : ", args.phi
print "\tNumber of particles : ", N
print "\tOutput file : ", args.output
print "\tRod radius : ", args.radius
print "\tRod length : ", args.length
print 

start = datetime.now()

random_orinet = True
#if args.l1 != 2.0 or args.l2 != 1.0:
  #random_orinet = True

p = Plane(args.lx, args.ly, N, lx, ly, sigma, args.length)
p.write(args.output)

print "Actual packing fraction for this box : ", len(p.particles)*Arod/area

end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print
