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

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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

class SphereCornea:

	def __init__(self, R, N, alpha,poly, v,boundary):
		self.R = R
		self.N = N
		self.alpha=alpha # half the cone angle of the cornea
		self.poly = poly
		self.boundary=boundary # as before
		self.Nbound = int(1.5*2.0*pi*self.R*sin(self.alpha)/2.0) # Dense boundary wall, linear packing fraction 1.5
		self.particles = [Particle(i,0) for i in range(N)]
		self.boundparticles = [Particle(i,1) for i in range(N,N+self.Nbound)]
		self.boundparticles.extend([Particle(i,2) for i in range(N+self.Nbound,N+2*self.Nbound)])
		self.__generate_posinside()
		self.__generate_posbound()
		self.__generate_posbound2()
		self.__generate_vel(v)
		self.__generate_director()


	def __generate_posbound(self):
		for i in range(self.Nbound):
			# fractional position along circumference (angle)
			sval = 1.0*i/self.Nbound*2*pi
			self.boundparticles[i].r=[self.R*cos(sval)*sin(self.alpha),self.R*sin(sval)*sin(self.alpha),self.R*cos(self.alpha)]
			# introduce a little noise there; no ratchets here
			self.boundparticles[i].R = uniform(1-0.5*self.poly,1+0.5*self.poly)
		
	# 2nd layer to prevent escapees
	def __generate_posbound2(self):
		dang=0.5/self.R
		for i in range(self.Nbound,2*self.Nbound):
			# fractional position along circumference (angle)
			sval = (1.0*i+0.5)/self.Nbound*2*pi
			self.boundparticles[i].r=[self.R*cos(sval)*sin(self.alpha+dang),self.R*sin(sval)*sin(self.alpha+dang),self.R*cos(self.alpha+dang)]
			# introduce a little noise there; no ratchets here
			self.boundparticles[i].R = uniform(1-0.5*self.poly,1+0.5*self.poly)
		
	def __generate_posinside(self):
		for i in range(self.N):
			u = uniform(cos(self.alpha),1)
			t = uniform(0,2*pi)
			x = sqrt(1-u*u)*cos(t)
			y = sqrt(1-u*u)*sin(t)
			z = u
			self.particles[i].r = [self.R*x,self.R*y,self.R*z]
			self.particles[i].R = uniform(1-0.5*self.poly,1+0.5*self.poly)

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
		# This is actually intriguing: make them point inwards, outwards or random
		for i in range(self.Nbound):
			if self.boundary==2:
				sval = 1.0*i/self.Nbound*2*pi+pi
				self.boundparticles[i].n=[cos(sval)*cos(self.alpha),sin(sval)*cos(self.alpha),sin(self.alpha)]
			elif self.boundary==1:
				sval = 1.0*i/self.Nbound*2*pi
				self.boundparticles[i].n=[cos(sval)*cos(self.alpha),sin(sval)*cos(self.alpha),-sin(self.alpha)]
			else:
				sval = uniform(0,2*pi)
				self.boundparticles[i].n=[cos(sval)*cos(self.alpha),sin(sval)*cos(self.alpha),sin(self.alpha)]

		
	def write(self,outfile):
		gentime = datetime.now()
		out = open(outfile,'w')
		out.write('# Total of %s particles\n' % str(self.N+self.Nbound))
		out.write('# Generated on : %s\n' % str(gentime))
		out.write('# id  type radius  x   y   z   vx   vy   vz   nx   ny   nz  omega\n')
		for p in self.particles:
			x, y, z = p.r
			vx, vy, vz = p.v
			nx, ny, nz = p.n
			out.write('%d  %d  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n' % (p.idx,p.tp,p.R,x,y,z,vx,vy,vz,nx,ny,nz,p.omega))
		for p in self.boundparticles:
			x, y, z = p.r
			vx, vy, vz = p.v
			nx, ny, nz = p.n
			out.write('%d  %d  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n' % (p.idx,p.tp,p.R,x,y,z,vx,vy,vz,nx,ny,nz,p.omega))
		out.close()
		
	def getPos(self,whichgroup):
		i=0
		rval=0.0
		if whichgroup==0:
			rval=np.empty((self.N,3))
			for p in self.particles:
				rval[i,:]=p.r
				i+=1
		elif whichgroup==1:
			rval=np.empty((self.Nbound,3))
			for i in range(self.Nbound):
				rval[i,:]=self.boundparticles[i].r 
		elif whichgroup==2:
			rval=np.empty((self.Nbound,3))
			for i in range(self.Nbound):
				rval[i,:]=self.boundparticles[self.Nbound+i].r 
		return rval


parser = argparse.ArgumentParser()
parser.add_argument("-R", "--radius", type=float, default=10.0, help="sphere radius")
parser.add_argument("-a", "--alpha", type=float, default=pi/2, help="cone angle")
parser.add_argument("-f", "--phi",  type=float, default=0.5, help="packing fraction")
parser.add_argument("-p", "--poly", type=float, default=0.0, help="polydispersity fraction")
parser.add_argument("-o", "--output", type=str, default='out.dat', help="output file")
parser.add_argument("-v", "--vavr", type=float, default=0.0, help="average velocity")
parser.add_argument("-b", "--boundary", type=int, default=0, help="boundary n: 0 random 1 outward 2 inward")
parser.add_argument("-vb", "--verbose", type=bool, default=False, help="plot configuration?")
args = parser.parse_args()

print
print "\tActive Particles on Curved Spaces (APCS)"
print "\tBuilding of an initial cornea configuration"
print 
print "\tSilke Henkes, Rastko Sknepnek"
print "\tUniversity of Aberdeen and University of Dundee"
print "\t(c) 2013, 2015"
print "\t----------------------------------------------"
print 
print "\tRadius : ", args.radius
print "\tPacking fraction : ", args.phi
print "\tAverage velocity : ", args.vavr
N=int(round(4.0*args.radius**2*args.phi*(1.0-cos(args.alpha))))
print "\tTotal number of particles : ", N
print "\tOutput file : ", args.output
print 

start = datetime.now()


s = SphereCornea(args.radius, N, args.alpha,args.poly, args.vavr,args.boundary)
s.write(args.output)

end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print

if args.verbose:
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	rinside=s.getPos(0)
	rstem=s.getPos(1)
	rblock=s.getPos(2)
	ax.scatter(rinside[:,0], rinside[:,1], rinside[:,2], zdir='z', c='b',s=20)
	ax.scatter(rstem[:,0], rstem[:,1], rstem[:,2], zdir='z', c='r',s=20)
	ax.scatter(rblock[:,0], rblock[:,1], rblock[:,2], zdir='z', c='g',s=20)
	plt.show()


