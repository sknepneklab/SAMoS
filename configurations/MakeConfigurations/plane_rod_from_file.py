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

# Utility code for building initial configuration on xy plane
# where particle coordinates and oriantion vectors are read from a file.

from datetime import *
from random import uniform 
from math import *
from scipy.integrate import dblquad
import numpy as np
import argparse


from particle import *


class Plane:
  
  def __init__(self,inp, v):
    self.inp = inp
    self.__read_pos_dir()
    self.__generate_vel(v)
        
  def __read_pos_dir(self):
    data = np.loadtxt(self.inp)
    min_x, max_x = np.min(data[:,0]), np.max(data[:,0])
    min_y, max_y = np.min(data[:,1]), np.max(data[:,1])
    self.Lx, self.Ly = max_x - min_x, max_y - min_y
    print "Lx = ", self.Lx
    print "Ly = ", self.Ly
    self.N = data.shape[0]
    self.particles = [Particle(i) for i in xrange(self.N)]
    for i in range(self.N):
      x, y, nx, ny = data[i]
      self.particles[i].r = [x-0.5*self.Lx,y-0.5*self.Ly,0.0]
      len_n = np.sqrt(nx*nx + ny*ny)
      self.particles[i].n = [nx/len_n,ny/len_n,0.0]

  def __generate_vel(self,vav=1.0):
    for p in self.particles:
      phi = uniform(0,2*pi)
      p.v = [vav*cos(phi),vav*sin(phi),0.0]
  
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
parser.add_argument("-i", "--input", type=str, help="inout file name")
parser.add_argument("-o", "--output", type=str, default='out.dat', help="output file")
parser.add_argument("-v", "--vavr", type=float, default=1.0, help="average velocity")
parser.add_argument("-p", "--phi", type=float, default=0.5, help="packing fraction")
parser.add_argument("-l", "--laspect",  type=float, default=4.0, help="aspect ratio")
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
print "\tinput file : ", args.input
print "\tAverage velocity : ", args.vavr
print "\tOutput file : ", args.output
start = datetime.now()

random_orinet = True
#if args.l1 != 2.0 or args.l2 != 1.0:
  #random_orinet = True

p = Plane(args.input, args.vavr)
# unscaled particle area based on the aspect ratio -- rescale for proper packing fraction
part_area = 1.0*(2.0*(args.laspect-1.0) + pi)
print part_area
part_real = p.Lx*p.Ly*args.phi/p.N
print part_real
rescale=np.sqrt(part_real/part_area)
print rescale
radii = []
types = []
lens = []
for i in xrange(len(p.particles)):
    radii.append(rescale)
    types.append(1)
    lens.append(rescale*2.0*args.laspect)
p.set_radius(radii)
p.set_type(types)
p.set_lens(lens)
p.write(args.output)

end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print
