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
# *   (c) 2015, 2016
# * 
# *   Author: Silke Henkes
# * 
# *   Department of Physics 
# *   Institute for Complex Systems and Mathematical Biology
# *   University of Aberdeen  
# * 
# *   (c) 2014, 2015, 2016
# *  
# *   This program cannot be used, copied, or modified without
# *   explicit written permission of the authors.
# * 
# * ***************************************************************

# Utility code for building random initial configuration on a torus

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

class Torus:
  
  def __init__(self, R, r, N, v):
    self.R = R
    self.r = r
    self.N = N
    self.particles = [Particle(i) for i in range(N)]
    self.__generate_pos()
    self.__generate_vel(v)
    self.__generate_director()
    
  def __generate_pos(self):
    for i in range(self.N):
      phi = uniform(0.0, 2*pi)
      theta = uniform(0.0, 2*pi)
      x = (self.R + self.r*cos(theta))*cos(phi)
      y = (self.R + self.r*cos(theta))*sin(phi)
      z = self.r*sin(theta)
      self.particles[i].r = [x,y,z]

  def __generate_vel(self,vav=1.0):
    for p in self.particles:
      x, y, z = p.r
      fact1 = sqrt(x*x + y*y)
      fact2 = (self.R - fact1)/fact1
      axis = np.array([-2*x*fact2,-2*y*fact2,2*z])
      v = np.array([uniform(-0.5,0.5), uniform(-0.5,0.5), uniform(-0.5,0.5)])
      v = np.dot(projection_matrix(axis),v)
      vlen = sqrt(sum(v**2))
      v *= vav/vlen
      theta = uniform(0,2*pi)
      v = np.dot(rotation_matrix(axis,theta),v)
      p.v = v
 
  def __generate_director(self):
    for p in self.particles:
      x, y, z = p.r
      fact1 = sqrt(x*x + y*y)
      fact2 = (self.R - fact1)/fact1
      axis = np.array([-2*x*fact2,-2*y*fact2,2*z])
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
parser.add_argument("-R", "--Radius", type=float, default=10.0, help="large radius of torus")
parser.add_argument("-r", "--radius", type=float, default=4.0, help="small radius of torus")
parser.add_argument("-n", "--density",  type=float, default=0.5, help="number density")
parser.add_argument("-o", "--output", type=str, default='out.dat', help="output file")
parser.add_argument("-v", "--vavr", type=float, default=1.0, help="average velocity")
args = parser.parse_args()

print
print "\tSoft Active Matter on Surfaces (SAMoS)"
print "\tBuilding of a random configuration on torus"
print 
print "\tRastko Sknepnek"
print "\tUniversity of Dundee"
print "\t(c) 2016"
print "\t----------------------------------------------"
print 
print "\tLarge radius : ", args.Radius
print "\tSmall radius : ", args.radius
print "\tNumber density : ", args.density
print "\tAverage velocity : ", args.vavr
A = 4*pi**2*args.Radius*args.radius
N=int(round(A*args.density))
print "\tTotal number of particles : ", N
print "\tOutput file : ", args.output
print 

start = datetime.now()


c = Torus(args.Radius, args.radius, N, args.vavr)
c.write(args.output)

end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print


