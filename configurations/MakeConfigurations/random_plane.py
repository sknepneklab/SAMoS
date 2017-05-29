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
from scipy.integrate import dblquad
import argparse


from particle import *

def bump(x,y,*args):
    A, a, b = args
    fact = -2 * A * exp(-x*x/(a*a) - y*y/(b*b))
    hx = x*fact/(a*a)
    hy = y*fact/(b*b)
    return sqrt(1.0 + hx*hx + hy*hy)


class Plane:
  
  def __init__(self, Lx, Ly, N, v, l=None, random_orinet = True):
    self.L = (Lx,Ly)
    self.N = N
    if l == None:
      self.__generate_pos()
    else:
      self.__generate_lattice(l)
    self.__generate_vel(v)
    self.__generate_director(random_orinet)
    
  def __generate_pos(self):
    self.particles = [Particle(i) for i in range(self.N)]
    for i in range(self.N):
      x = uniform(-0.5*self.L[0],0.5*self.L[0])
      y = uniform(-0.5*self.L[1],0.5*self.L[1])
      z = 0
      self.particles[i].r = [x,y,z]

  def __generate_lattice(self,l):
    Nx = int(floor(sqrt(self.N/(l+1))))
    Ny = (l+1)*Nx
    self.particles = [Particle(i) for i in range(Nx*Ny)]
    idx = 0
    for i in range(Nx):
      x = -0.5*self.L[0] + (0.5+i)*(l+1)
      for j in range(Ny):
        y = -0.5*self.L[1] + j
        self.particles[idx].r = [x,y,0]
        idx += 1

  def __generate_vel(self,vav=1.0):
    for p in self.particles:
      phi = uniform(0,2*pi)
      p.v = [vav*cos(phi),vav*sin(phi),0.0]
  
  def __generate_director(self, random_orinet):
    for p in self.particles:
      if random_orinet:
        phi = uniform(0,2*pi)
      else:
        if uniform(-0.5,0.5) < 0.0:
          phi = 0.0 
        else:
          phi = pi
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
      
  def set_mass(self,masses):
    for i in xrange(len(masses)):
      self.particles[i].mass = masses[i]
  
  def gaussian_bump(self,A=1.0,a=1.0,b=1.0):
    for i in xrange(len(self.particles)):
      x, y, z = self.particles[i].r
      z = A*exp(-x*x/(a*a) - y*y/(b*b))
      self.particles[i].r = [x,y,z]
  
  def write(self,outfile):
    gentime = datetime.now()
    out = open(outfile,'w')
    out.write('# Total of %d particles\n' % self.N)
    out.write('# Generated on : %s\n' % str(gentime))
    #out.write('# id  type radius  x   y   z   vx   vy   vz   nx   ny   nz  omega l\n')
    out.write('keys: id  type radius mass x   y   z   vx   vy   vz   nx   ny   nz  length\n')
    for p in self.particles:
      x, y, z = p.r
      vx, vy, vz = p.v
      nx, ny, nz = p.n
      out.write('%d  %d  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n' % (p.idx,p.tp,p.R,p.mass,x,y,z,vx,vy,vz,nx,ny,nz,p.l))
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
parser.add_argument("--l1",  type=float, default=2.0, help="length of rod of type 1 (in terms of particle radius)")
parser.add_argument("--l2",  type=float, default=1.0, help="length of rod of type 2")
parser.add_argument("--amplitude",  type=float, default=1.0, help="Gaussian bump amplitude")
parser.add_argument("-a", "--ga",  type=float, default=1.0, help="Gaussian bump x width")
parser.add_argument("-b", "--gb",  type=float, default=1.0, help="Gaussian bump y width")
parser.add_argument("-l","--lattice", action='store_true', help="make lattice")
parser.add_argument("-g","--gaussian", action='store_true', help="make Gaussian bump")
parser.add_argument("-r","--rods", action='store_true', help="make rods")
parser.add_argument("-m", "--mass",  type=float, default=1.0, help="particle mass (same for all particles)")
args = parser.parse_args()

if args.gaussian:
    V = dblquad(bump,-0.5*args.lx,0.5*args.lx,lambda x: -0.5*args.ly, lambda x: 0.5*args.ly,args=(args.amplitude,args.ga,args.gb))[0]
    print 'Area : ', V
else:
    V = args.lx*args.ly

if args.rods:
    N = int(round(0.5*args.phi*V/((args.eta*args.a1*(args.l1 + 0.5*pi*args.a1)+(1-args.eta)*args.a2**2*(args.l2 + 0.5*pi*args.a2)))))
else:
    N = int(round(1.0/pi*V*args.phi/(args.eta*args.a1**2+(1-args.eta)*args.a2**2)))

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
print "\tAverage velocity : ", args.vavr
print "\tOutput file : ", args.output
print "\tFraction of particles of type 1 : ", args.eta
print "\tRadius of particles of type 1 : ", args.a1
print "\tRadius of particles of type 2 : ", args.a2
print "\tLength of particles of type 1 (in units of radius) : ", args.l1
print "\tLength of particles of type 2 (in units of radius) : ", args.l2
print 

start = datetime.now()

random_orinet = True
#if args.l1 != 2.0 or args.l2 != 1.0:
  #random_orinet = True

if args.lattice:
  p = Plane(args.lx, args.ly, N, args.vavr, args.l1, random_orinet=random_orinet)
else:
  p = Plane(args.lx, args.ly, N, args.vavr, random_orinet=random_orinet)
radii = []
types = []
lens = []
masses = []
for i in xrange(len(p.particles)):
  if i < args.eta*N: 
    radii.append(args.a1)
    types.append(1)
    lens.append(args.l1)
    masses.append(args.mass)
  else: 
    radii.append(args.a2)
    types.append(2)
    lens.append(args.l2)
    masses.append(args.mass)
p.set_radius(radii)
p.set_type(types)
p.set_lens(lens)
p.set_mass(masses)
if args.gaussian:
  p.gaussian_bump(args.amplitude,args.ga,args.gb)
p.write(args.output)

end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print
