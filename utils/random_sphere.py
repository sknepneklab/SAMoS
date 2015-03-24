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

class Sphere:
  
  def __init__(self, R, N, v):
    self.R = R
    self.N = N
    self.particles = [Particle(i) for i in range(N)]
    self.__generate_pos()
    self.__generate_vel(v)
    self.__generate_director()
    
  def __generate_pos(self):
    for i in range(self.N):
      u = uniform(-1,1)
      t = uniform(0,2*pi)
      x = sqrt(1-u*u)*cos(t)
      y = sqrt(1-u*u)*sin(t)
      z = u
      self.particles[i].r = [self.R*x,self.R*y,self.R*z]

  def __generate_vel(self,vav=1.0):
    for p in self.particles:
      axis = np.array(p.r)
      v = np.array([uniform(-0.5,0.5), uniform(-0.5,0.5), uniform(-0.5,0.5)])
      v = np.dot(projection_matrix(axis),v)
      vlen = sqrt(sum(v**2))
      v *= vav/vlen
      theta = uniform(0,2*pi)
      v = np.dot(rotation_matrix(axis,theta),v)
      p.v = v
 
  def __generate_director(self):
    for p in self.particles:
      axis = np.array(p.r)
      n = np.array([uniform(-0.5,0.5), uniform(-0.5,0.5), uniform(-0.5,0.5)])
      n = np.dot(projection_matrix(axis),n)
      nlen = sqrt(sum(n**2))
      n *= 1.0/nlen
      theta = uniform(0,2*pi)
      n = np.dot(rotation_matrix(axis,theta),n)
      p.n = n
 
  def set_lens(self,lens):
    for i in xrange(len(lens)):
      self.particles[i].l = lens[i]
      
  def write(self,outfile):
    gentime = datetime.now()
    out = open(outfile,'w')
    out.write('# Total of %d particles\n' % self.N)
    out.write('# Generated on : %s\n' % str(gentime))
    out.write('# id  type radius  x   y   z   vx   vy   vz   nx   ny   nz  omega  l\n')
    for p in self.particles:
      x, y, z = p.r
      vx, vy, vz = p.v
      nx, ny, nz = p.n
      out.write('%d  %d  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n' % (p.idx,p.tp,p.R,x,y,z,vx,vy,vz,nx,ny,nz,p.omega,p.l))
    out.close()
    


parser = argparse.ArgumentParser()
parser.add_argument("-R", "--radius", type=float, default=10.0, help="sphere radius")
parser.add_argument("-f", "--phi",  type=float, default=0.5, help="packing fraction")
parser.add_argument("-o", "--output", type=str, default='out.dat', help="output file")
parser.add_argument("-v", "--vavr", type=float, default=1.0, help="average velocity")
parser.add_argument("-l", "--length", type=float, default=2.0, help="rod length")
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
print "\tPacking fraction : ", args.phi
print "\tAverage velocity : ", args.vavr
N=int(round(4.0*args.radius**2*args.phi))
print "\tTotal number of particles : ", N
print "\tOutput file : ", args.output
print 

start = datetime.now()


s = Sphere(args.radius, N, args.vavr)
lens = [args.length for i in range(N)]
s.set_lens(lens)
s.write(args.output)

end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print


