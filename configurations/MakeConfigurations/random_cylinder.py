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

# Utility code for building random initial configuration on a sphere

from datetime import *
from random import uniform 
from math import *
import argparse
import numpy as np

from particle import *

def projection_matrix(axis):
  n = axis/np.sqrt(np.dot(axis,axis))
  return (np.identity(3) - np.outer(n,n))


def rotation_matrix(axis,theta):
  n = axis/np.sqrt(np.dot(axis,axis))
  a = np.cos(theta/2)
  b,c,d = -n*np.sin(theta/2)
  return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                  [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                  [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])

class Cylinder:
  
  def __init__(self, R, H, N, v, A, n):
    self.R = R
    self.H = H
    self.N = N
    self.A = A
    self.n = n
    self.particles = [Particle(i) for i in range(N)]
    self.__generate_pos()
    self.__generate_vel(v)
    self.__generate_director()
    
  def __generate_pos(self):
    for i in range(self.N):
      z = uniform(-0.5*self.H,0.5*self.H)
      phi = uniform(0.0, 2*pi)
      R = self.R + self.A*sin(2.0*pi/self.H*self.n*z)
      x = R*cos(phi)
      y = R*sin(phi)
      self.particles[i].r = [x,y,z]

  def __generate_vel(self,vav=1.0):
    for p in self.particles:
      zgrad = -4.0*self.A*self.n*pi/self.H*cos(2.0*pi/self.H*self.n*p.r[2])*(self.R+self.A*sin(2.0*pi/self.H*self.n*p.r[2]))
      axis = np.array([p.r[0],p.r[1],zgrad])
      v = np.array([uniform(-0.5,0.5), uniform(-0.5,0.5), uniform(-0.5,0.5)])
      v = np.dot(projection_matrix(axis),v)
      vlen = sqrt(sum(v**2))
      v *= vav/vlen
      theta = uniform(0,2*pi)
      v = np.dot(rotation_matrix(axis,theta),v)
      p.v = v
 
  def __generate_director(self):
    for p in self.particles:
      zgrad = -4.0*self.A*self.n*pi/self.H*cos(2.0*pi/self.H*self.n*p.r[2])*(self.R+self.A*sin(2.0*pi/self.H*self.n*p.r[2]))
      axis = np.array([p.r[0],p.r[1],zgrad])
      n = np.array([uniform(-0.5,0.5), uniform(-0.5,0.5), uniform(-0.5,0.5)])
      n = np.dot(projection_matrix(axis),n)
      nlen = sqrt(sum(n**2))
      n *= 1.0/nlen
      theta = uniform(0,2*pi)
      n = np.dot(rotation_matrix(axis,theta),n)
      p.n = n
 
  def write(self,outfile):
    gentime = datetime.now()
    out = open(outfile,'w')
    out.write('# Total of %d particles\n' % self.N)
    out.write('# Generated on : %s\n' % str(gentime))
    out.write('# id  type radius  x   y   z   vx   vy   vz   nx   ny   nz   omega\n')
    for p in self.particles:
      x, y, z = p.r
      vx, vy, vz = p.v
      nx, ny, nz = p.n
      out.write('%d  %d  %f %f  %f  %f  %f  %f  %f  %f  %f  %f   %f\n' % (p.idx,p.tp,p.R,x,y,z,vx,vy,vz,nx,ny,nz,p.omega))
    out.close()
    


parser = argparse.ArgumentParser()
parser.add_argument("-R", "--radius", type=float, default=10.0, help="cylinder radius")
parser.add_argument("-H", "--height", type=float, default=20.0, help="cylinder height")
parser.add_argument("-A", "--amplitude", type=float, default=0.0, help="amplitude for a wavy cylinder")
parser.add_argument("-n", "--n", type=float, default=1.0, help="number of nodes (for a wavy cylinder)")
parser.add_argument("-f", "--phi",  type=float, default=0.5, help="packing fraction")
parser.add_argument("-o", "--output", type=str, default='out.dat', help="output file")
parser.add_argument("-v", "--vavr", type=float, default=1.0, help="average velocity")
args = parser.parse_args()

print
print "\tActive Particles on Curved Spaces (APCS)"
print "\tBuilding of a random spherical configuration"
print 
print "\tRastko Sknepnek"
print "\tUniversity of Dundee"
print "\t(c) 2013"
print "\t----------------------------------------------"
print 
print "\tRadius : ", args.radius
print "\tHeight : ", args.height
print "\tPacking fraction : ", args.phi
print "\tAverage velocity : ", args.vavr
print "\tAmplitude (wavy cylinder) : ", args.amplitude
print "\tNumber of nodes (wavy cylinder) : ", args.n
N=int(round(2.0*args.radius*args.height*args.phi))
print "\tTotal number of particles : ", N
print "\tOutput file : ", args.output
print 

start = datetime.now()


c = Cylinder(args.radius, args.height, N, args.vavr, args.amplitude, args.n)
c.write(args.output)

end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print


