# ################################################################
#  
#    Active Particles on Curved Spaces (APCS)
#    
#   Author: Silke Henkes
#   
#    ICSMB, Department of Physics
#    University of Aberdeen
#
#    Author: Rastko Sknepnek
#   
#    Division of Physics
#    School of Engineering, Physics and Mathematics
#    University of Dundee
#    
#    (c) 2013, 2014, 2015
#    
#    This program cannot be used, copied, or modified without
#    explicit permission of the author.
#  
# ################################################################

# Main Python library for computing:
# angular profiles for a rotating band state
# Defects for both nematic and polar symmetry states
# Tesselations of the packing using the Ball-Blumenfeld / dual to the Maxwell-Cremona
# Visualization output for Paraview

# Structure:
# External utility function to read parameters

from read_data import *
#from read_param import *
from op import *
#from inertia import *
from glob import glob
from datetime import *
from random import uniform 
from math import *
import numpy as np
import argparse
import scipy.spatial.distance as sd

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
#from matplotlib import rc
import matplotlib
from mpl_toolkits.mplot3d import Axes3D

import vtk

# setting global plotting parameters
#matplotlib.rcParams['text.usetex'] = 'true'
matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['axes.linewidth'] = 2
matplotlib.rcParams['xtick.major.size'] = 8
matplotlib.rcParams['ytick.major.size'] = 8
#matplotlib.rcParams['font.size']=40.0
#matplotlib.rcParams['legend.fontsize']=22.0
matplotlib.rcParams['font.size']=28
matplotlib.rcParams['legend.fontsize']=14


cdict = {'red':   [(0.0,  0.25, 0.25),
           (0.3,  1.0, 1.0),
                   (0.5,  0.4, 0.0),
                   (1.0,  0.0, 0.0)],

         'green': [(0.0,  0.0, 0.0),
           (0.25,  0.0, 0.5),
                   (0.5, 1.0, 1.0),
                   (0.75,  0.5, 0.0),
                   (1.0,  0.0, 0.0)],

         'blue':  [(0.0,  0.0, 0.0),
                   (0.5,  0.0, 0.0),
                   (0.7, 1.0, 1.0),
                   (1.0,  0.25, 0.25)]}

# Default geometry class: 3 dimensional unconstrained space
class Geometry:
	def __init__(self,manifold):
		try:
			self.manifold=manifold
			print "Created new geometry " + manifold
		except:
			pass
		
	def RotateMatrixVectorial(self,axis,theta):
		axlen=np.sqrt(axis[:,0]**2+axis[:,1]**2+axis[:,2]**2)
		#print axlen
		axis[:,0]=axis[:,0]/axlen
		axis[:,1]=axis[:,1]/axlen
		axis[:,2]=axis[:,2]/axlen
		a=np.cos(theta/2)
		b=-axis[:,0]*np.sin(theta/2)
		c=-axis[:,1]*np.sin(theta/2)
		d=-axis[:,2]*np.sin(theta/2)
		rotmat=np.empty((len(axis[:,0]),3,3))
		rotmat[:,0,0]=a*a+b*b-c*c-d*d
		rotmat[:,0,1]=2*(b*c-a*d)
		rotmat[:,0,2]=2*(b*d+a*c)
		rotmat[:,1,0]=2*(b*c+a*d)
		rotmat[:,1,1]=a*a+c*c-b*b-d*d
		rotmat[:,1,2]=2*(c*d-a*b)
		rotmat[:,2,0]=2*(b*d-a*c)
		rotmat[:,2,1]=2*(c*d+a*b)
		rotmat[:,2,2]=a*a+d*d-b*b-c*c
		return rotmat

	def RotateVectorial(self,v,n,phi):
		vrot=np.empty(np.shape(v))
		np.shape(vrot)
		rotmat=rotate_matrix_vectorial(n,phi)
		np.shape(rotmat)
		vrot[:,0]=rotmat[:,0,0]*v[:,0]+rotmat[:,0,1]*v[:,1]+rotmat[:,0,2]*v[:,2]
		vrot[:,1]=rotmat[:,1,0]*v[:,0]+rotmat[:,1,1]*v[:,1]+rotmat[:,1,2]*v[:,2]
		vrot[:,2]=rotmat[:,2,0]*v[:,0]+rotmat[:,2,1]*v[:,1]+rotmat[:,2,2]*v[:,2]
		return vrot

	# Default: parallel transport just transports parallel; i.e. it does nothing in a flat geometry    
	def ParallelTransport(self,r1,r2,a2):
		return a2
		
	def ParallelTransportSingle(self,r1,r2,a2):
		return a2
		
	# Default: just the cartesian distance
	def GeodesicDistance(r1,r2):
		return np.sqrt(np.sum((r2-r1)**2))
        
# Spherical geometry   
class GeometrySphere(Geometry):
	def __init__(self,param):
		self.R=param.r
		print "Created new geometry sphere with radius " + str(R)
		super(GeometrySphere,self).__init__('sphere')
		
	# Fully vectorial version of parallel transport
	# 1.determine the cross product of the origins
	# 2.compute the magnitude of all the origin and cross vectors
	# 3.Compute the dot product of the origins
	# 4.The rotation axis is the direction of the cross product
	# 5.The rotation angle is the angle between the origin vectors, extracted from the dot product
	def ParallelTransport(self,r1,r2,a2):
		r2_x_r1=np.cross(r2,r1)
		#len_r2_x_r1=np.sqrt(r2_x_r1[:,0]**2+r2_x_r1[:,1]**2+r2_x_r1[:,2]**2)
		len_r2_x_r1=np.sqrt(np.sum(r2_x_r1**2,axis=1)) 
		#lenr1=np.sqrt(r1[:,0]**2+r1[:,1]**2+r1[:,2]**2)
		lenr1=np.sqrt(np.sum(r1**2,axis=1))
		#lenr2=np.sqrt(r2[:,0]**2+r2[:,1]**2+r2[:,2]**2)
		lenr2=np.sqrt(np.sum(r2**2,axis=1))
		dot_r1r2=r1[:,0]*r2[:,0]+r1[:,1]*r2[:,1]+r1[:,2]*r2[:,2]
		n=np.empty(np.shape(r1))
		n = r2_x_r1/len_r2_x_r1
		#n[:,0] = r2_x_r1[:,0]/len_r2_x_r1
		#n[:,1] = r2_x_r1[:,1]/len_r2_x_r1
		#n[:,2] = r2_x_r1[:,2]/len_r2_x_r1
		phi = np.arccos(dot_r1r2/(lenr1*lenr2))
		a2trans=rotate_vectorial(a2,n,-phi)
		return a2trans
		
	# same thing for one vector and a set (i.e. a particle and its neigbours)
	def ParallelTransportSingle(self,r1,r2,a2):
		r2_x_r1=np.cross(r2,r1)
		len_r2_x_r1=np.sqrt(np.sum(r2_x_r1**2,axis=1)) 
		lenr1=np.sqrt(np.sum(r1**2,axis=1))
		lenr2=np.sqrt(np.sum(r2**2,axis=1))
		dot_r1r2=np.dot(r1,r2)
		n=np.empty(np.shape(r1))
		n = r2_x_r1/len_r2_x_r1
		phi = np.arccos(dot_r1r2/(lenr1*lenr2))
		a2trans=rotate_vectorial(a2,n,-phi)
		return a2trans
		
	# Determine the Euler angles, essentially. Find theta and phi for each particle,
	def TangentBundle(self,rval):
		rhat=((rval).transpose()/(np.sqrt(np.sum(rval**2,axis=1))).transpose()).transpose()
		# Angle theta with the z axis. arccos is between 0 and pi, so that's ok already
		theta=np.arccos(rhat[:,2])
		# From the euler angles: rx = sin theta cos phi
		# Choosing correct quadrant through the sign of ry=sin theta sin phi
		phi=np.sign(rhat[:,1]/(np.sin(theta)))*np.arccos(rhat[:,0]/(np.sin(theta)))
		# The other two of our trio of local coordinate vectors
		etheta = np.empty(np.shape(rval))
		etheta[:,0]=np.cos(theta)*np.cos(phi)
		etheta[:,1]=np.cos(theta)*np.sin(phi)
		etheta[:,2]=-np.sin(theta)
		ephi=np.empty(np.shape(rval))
		ephi[:,0]=-np.sin(phi)
		ephi[:,1]=np.cos(phi)
		ephi[:,2]=0
		return theta,phi,etheta,ephi
        
# Plane with periodic boundary conditions. By default, the plane is along x and y
class GeometryPeriodicPlane(Geometry):
	def __init__(self,param):
		self.Lx=param.lx
		self.Ly=param.ly
		print "Created new geometry periodic plane with Lx = " + str(self.Lx) + " and Ly = " +str(self.Ly)
		super(GeometryPeriodicPlane,self).__init__('plane')
		
	def TangentBundle(self,rval):
		x=rval[:,0]
		y=rval[:,1]
		ex[:,0]=1.0*np.ones(shape(rval))
		ex[:,1]=1.0*np.zeros(shape(rval))
		ex[:,2]=1.0*np.zeros(shape(rval))
		ey[:,0]=1.0*np.zeros(shape(rval))
		ey[:,1]=1.0*np.ones(shape(rval))
		ey[:,2]=1.0*np.zeros(shape(rval))
		return x,y,ex,ey
		
	# Just the cartesian distance in the plane, modulo periodic boundary conditions
	def GeodesicDistance(r1,r2):
		drx=r2[0]-r1[0]
		drx-=self.Lx*round(drx)
		dry=r2[1]-r1[1]
		dry-=self.Ly*round(dry)
		return np.sqrt(drx**2+dry**2)

class GeometryTube(Geometry):
	def __init__(self,param):
		self.R=param.const_param['r']
		print "Created new geometry tube with radius = " + str(self.R)
		super(GeometryTube,self).__init__('tube')
		print "ERROR: Tube geometry has not yet been implemented! Geometry will default to 3d space."
		
class GeometryPeanut(Geometry):
	def __init__(self,param):
		super(GeometryPeanut,self).__init__('peanut')
		print "ERROR: Peanut geometry has not yet been implemented! Geometry will default to 3d space."
		
class GeometryHourglass(Geometry):
	def __init__(self,param):
		self.R=param.const_param['R']
		self.A=param.const_param['A']
		print "Created new geometry tube with radius = " + str(self.R) + " and amplitude " +str(self.A)
		super(GeometryHourglass,self).__init__('hourglass')
		print "ERROR: Hourglass geometry has not yet been implemented! Geometry will default to 3d space."
		
class Configuration:
	def __init__(self,param,filename):
		self.param=param
		# Read the local data
		geometries={'sphere':GeometrySphere,'plane':GeometryPeriodicPlane,'none':Geometry,'tube':GeometryTube,'peanut':GeometryPeanut,'hourglass':GeometryHourglass}
		print "Processing file : ", filename
		data = ReadData(filename)
		x, y, z = np.array(data.data[data.keys['x']]), np.array(data.data[data.keys['y']]), np.array(data.data[data.keys['z']])
		vx, vy, vz = np.array(data.data[data.keys['vx']]), np.array(data.data[data.keys['vy']]), np.array(data.data[data.keys['vz']])
		nx, ny, nz = np.array(data.data[data.keys['nx']]), np.array(data.data[data.keys['ny']]), np.array(data.data[data.keys['nz']])
		self.monodisperse=False
		if not data.keys.has_key('radius'): 
			# MISSING: read them in from the initial configuration
			self.radius = np.array([1.0 for i in range(self.param.N)])
			self.monodisperse=True
			self.sigma=1.0
		else: 
			self.radius = data.data[data.keys['radius']]	
		if data.keys.has_key('type'):
			self.ptype = data.data[data.keys['type']]
		if data.keys.has_key('flag'):
			self.flag = data.data[data.keys['flag']]
		self.rval = np.column_stack((x,y,z))
		self.vval = np.column_stack((vx,vy,vz))
		self.nval = np.column_stack((nx,ny,nz))
		# Create the right geometry environment (TBC):
		self.geom=geometries[param.constraint](param)
		
	# Tangent bundle: Coordinates in an appropriate coordinate system defined on the manifold
	# and the coordinate vectors themselves, for all the particles
	def getTangentBundle(self):
		self.x1,self.x2,self.e1,self.e2=self.geom.TangentBundle(self.rval)
		return self.x1,self.x2,self.e1,self.e2
	
	def rotateFrame(self,axis,rot_angle):
		self.rval = self.geom.RotateVectorial(self.rval,axis0,-rot_angle)
		self.vval = self.geom.RotateVectorial(self.vval,axis0,-rot_angle)
		self.nval = self.geom.RotateVectorial(self.nval,axis0,-rot_angle)
		self.nval=((self.nval).transpose()/(np.sqrt(np.sum(self.nval**2,axis=1))).transpose()).transpose()
		self.vel = np.sqrt(self.vval[:,0]**2 + self.vval[:,1]**2 + self.vval[:,2]**2)
        

	def compute_energy_and_pressure(self):
		eng = np.zeros(self.param.N)
		press = np.zeros(self.param.N)
		stress = np.zeros((self.param.N,3,3))
		# Establish how many types of particles we have
		# If it's only one, use the current setup
		if self.param.ntypes==1:
			if self.param.potential=='soft':
				if self.monodisperse:
					dmax=4*self.sigma**2
				for i in range(len(r)):
				#for i in range(10):
					#dist=np.sum((r-r[i,:])**2,axis=1)
					dist=self.geom.GeodesicDistance(r,r[i,:])
					if self.monodisperse: 
						neighbours=[index for index,value in enumerate(dist) if value <dmax]
					else:
						neighbours=[index for index,value in enumerate(dist) if value < (radius[i]+radius[j])**2]
					neighbours.remove(i)
					dr=np.sqrt(dist[neighbours])
					diff=radius[i]+radius[j]-dr
					fact = 0.5*self.param.pot_param['k']*diff
					eng_val = fact*diff
					press_val = fact*dr
					# Stress (force moment) has to be element by element) r_a F_b = -k r_a dist_b 
					drvec=r[neighbours,:]-r[i,:]
					Fvec=k*((diff/dr).transpose()*(drvec).transpose()).transpose()
					for u in range(3):
						for v in range(3):
							stress[neighbours,u,v]+=0.5*drvec[:,u]*Fvec[:,v]
					eng[neighbours]+=eng_val
					press[neighbours]+=press_val
			elif self.param.potential=='morse':
				# We are morse by hand right now ...
				D=self.param.pot_param['D']
				re=self.param.pot_param['re']
				a=self.param.pot.param['a']
				if self.monodisperse:
					dmax=16*self.sigma**2
				for i in range(self.param.N):
				#for i in range(10):
					#dist=np.sum((r-r[i,:])**2,axis=1)
					dist=self.geom.GeodesicDistance(r,r[i,:])
					if self.monodisperse: 
						neighbours=[index for index,value in enumerate(dist) if value <dmax]
					else:
						neighbours=[index for index,value in enumerate(dist) if value < (2*radius[i]+2*radius[j])**2]
					neighbours.remove(i)
					dr=np.sqrt(dist[neighbours])
					eng_val=D*(1-np.exp(-a*(dr-re)))**2
					fnorm=-2*a*D*np.exp(-a*(dr-re))*(1-np.exp(-a*(dr-re)))
					drvec=r[neighbours,:]-r[i,:]
					Fvec=((fnorm/dr).transpose()*(drvec).transpose()).transpose()
					press_val=fnorm*dr
					for u in range(3):
						for v in range(3):
							stress[neighbours,u,v]+=0.5*drvec[:,u]*Fvec[:,v]
					eng[neighbours]+=eng_val
					press[neighbours]+=press_val
			elif self.param.potential=='gaussian':
				print "Warning! Gaussian interaction has not yet been implemented! Returning zero energy and stresses"
			else:
				print "Warning! Unknown interaction type! Returning zero energy and stresses"
		else:
			print "Warning! Multiple types of particles interactin have not yet been implemented! Returning zero energy and stresses"
		return [eng, press, stress]
	
	# Flat case statistics (or other geometry statistics, if desired)
	def getStats(self,debug=False):
		ez = np.array([0,0,1])  # lab frame z-axis
		# The order parameter with v_0 still in it. Normalize in final polish
		orderparV=np.sum(vval,axis=0)/len(vval)
		orderpar=np.sum(nval,axis=0)/len(nval)
		print orderpar
		print orderparV
		direction = orderpar/np.linalg.norm(orderpar)
		directionV = orderparV/np.linalg.norm(orderparV)
		axisorth= np.cross(direction,directionV)
		axisval=np.linalg.norm(axisorth)
		alpha=np.arcsin(axisval)
		axisorth=axisorth/axisval
		axisnorm=np.cross(ez,directionV)
		axisnorm/=np.linalg.norm(axisnorm)
		
		print directionV
		print axisorth
		
		vel = np.sqrt(self.vval[:,0]**2 + self.vval[:,1]**2 + self.vval[:,2]**2)
		velnorm=((self.vval).transpose()/(vel).transpose()).transpose()
		
		eng, press,stress = self.compute_energy_and_pressure()
		print np.shape(stress)
		# Project the stresses into the e,theta,phi components. The rr component hast to be 0, and the r cross components
		# belong to the projection. So they are not all that interesting. 
		# We want the theta theta, theta phi, phi theta ant phi phi components (implicitly testing symmetries ...)
		# I give up on the notation. Stress is (N,3,3), the axes are (N,3). We want e_i sigma_ij e_j
		s_tt=np.sum(axisnorm*np.einsum('kij,j->ki',stress,axisnorm),axis=1)
		s_tp=np.sum(axisnorm*np.einsum('kij,j->ki',stress,directionV),axis=1)
		s_pt=np.sum(directionV*np.einsum('kij,j->ki',stress,axisnorm),axis=1)
		s_pp=np.sum(directionV*np.einsum('kij,j->ki',stress,directionV),axis=1)
		print np.shape(s_tt)
		# Mean density really makes no sense? Determined by the initial conditions in periodic boundary conditions.
		# I do not wish to set up artificial bins in a translationally invariant system
		vel_av=np.mean(vel)
		eng_av=np.mean(eng)
		press_av=np.mean(press)
		s_tt_av=np.mean(s_tt)
		s_tp_av=np.mean(s_tp)
		s_pt_av=np.mean(s_pt)
		s_pp_av=np.mean(s_pp)
		
		# Debugging output
		if debug==True:
			fig = plt.figure()
			ax = fig.add_subplot(111, projection='3d')
			ax.scatter(rval[:,0], rval[:,1], rval[:,2], zdir='z', c='b')
			
		return [vel_av,eng_av,press_av,s_tt_av,s_tp_av,s_pt_av,s_pp_av,alpha,direction,directionV,orderpar,orderparV]

 
# Profiles for the spherical case. 
# Warning: Implicitly spherical geometry. Results in other geometries will be garbage (should not crash though)
class Profiles:
	def __init__(self,conf,geom):
		self.conf=conf
		self.geom=geom
		#self.rval=self.conf.rval
		ez = np.array([0,0,1])  # lab frame z-axis
		# Simply get the axis as the mean crossproduct or r and v; assuming alignment. This should also not flip.
		self.direction=np.sum(np.cross(conf.rval,conf.vval),axis=0)
		self.orderpar=direction/len(conf.rval)
		print orderpar
		self.direction = self.direction/np.linalg.norm(self.direction)
		axis = np.cross(self.direction,ez)
		axis = axis/np.linalg.norm(axis)
		rot_angle = acos(np.dot(self.direction,ez))
		axis0 = np.empty(np.shape(self.conf.rval))
		axis0[:,0] = axis[0]
		axis0[:,1] = axis[1]
		axis0[:,2] = axis[2]
		self.conf.rotateFrame(axis0,-rot_angle)
		
		#rhat=((rval).transpose()/(np.sqrt(np.sum(rval**2,axis=1))).transpose()).transpose()
		#self.vel = np.sqrt(self.conf.vval[:,0]**2 + self.conf.vval[:,1]**2 + self.conf.vval[:,2]**2)
		velnorm=((self.conf.vval).transpose()/(vel).transpose()).transpose()
  
		self.theta,self.phi,self.etheta,self.ephi=geom.TangentBundle(self,rval)
		# Alpha, the angle between the local polarity and the equator; here represented by ephi
		self.alpha=-np.arcsin(np.sum(nval*etheta, axis=1))
		# Same thing for the velocity
		# No - add pi/2 to get something that does not add up to zero 
		#alpha_v=np.arccos(np.sum(velnorm*etheta, axis=1))
  
	def getProfiles(self,debug)
		eng, press,stress = conf.compute_energy_and_pressure()
		# Project the stresses into the e,theta,phi components. The rr component hast to be 0, and the r cross components
		# belong to the projection. So they are not all that interesting. 
		# We want the theta theta, theta phi, phi theta ant phi phi components (implicitly testing symmetries ...)
		# I give up on the notation. Stress is (N,3,3), the axes are (N,3). We want e_i sigma_ij e_j
		s_tt=np.sum(self.etheta*np.einsum('kij,kj->ki',stress,self.etheta),axis=1)
		s_tp=np.sum(self.etheta*np.einsum('...ij,...j->...i',stress,self.ephi),axis=1)
		s_pt=np.sum(self.ephi*np.einsum('...ij,...j->...i',stress,self.etheta),axis=1)
		s_pp=np.sum(self.ephi*np.einsum('...ij,...j->...i',stress,self.ephi),axis=1)
		
		# Setting up the binning. I changed this to go from -pi/2 to pi/2 consistently. This maybe makes less pretty pictures,
		# but the edges are going to be a lot cleaner. Also only one bin to handle accross multiple v0/J.
		# Can always rebin to less resolution if necessary
		# Position angle with the z axis
		theta_bin=np.linspace(0,np.pi,nbin+1)
		dtheta=theta_bin[1]-theta_bin[0]
		theta_out=theta_bin[:nbin]+dtheta/2-np.pi/2
		
		rho_profile, bin_edges = np.histogram(self.theta, bins=theta_bin,density=True)
		isdata=[index for index,value in enumerate(rho_profile) if (value >0)]
		normz=2*np.pi*self.conf.radius*abs(np.cos(theta_out))
		rho_profile[isdata]=rho_profile[isdata]/normz[isdata]
		rho_profile/=np.mean(rho_profile)
		vel_profile=np.zeros(np.shape(rho_profile))
		eng_profile=np.zeros(np.shape(rho_profile))
		press_profile=np.zeros(np.shape(rho_profile))
		s_tt_profile=np.zeros(np.shape(rho_profile))
		s_tp_profile=np.zeros(np.shape(rho_profile))
		s_pt_profile=np.zeros(np.shape(rho_profile))
		s_pp_profile=np.zeros(np.shape(rho_profile))
		alpha_profile=np.zeros(np.shape(rho_profile))
		alpha_v_profile=np.zeros(np.shape(rho_profile))
		for idx in range(nbin):
			inbin=[index for index,value in enumerate(self.theta) if (value >= theta_bin[idx]  and value<=theta_bin[idx+1])]
			#print len(inbin)
			if len(inbin)>0:
			vel_profile[idx]=np.mean(self.conf.vel[inbin])
			eng_profile[idx]=np.mean(eng[inbin])
			press_profile[idx]=np.mean(press[inbin])
			s_tt_profile[idx]=np.mean(s_tt[inbin])
			s_tp_profile[idx]=np.mean(s_tp[inbin])
			s_pt_profile[idx]=np.mean(s_pt[inbin])
			s_pp_profile[idx]=np.mean(s_pp[inbin])
			alpha_profile[idx]=np.mean(self.alpha[inbin])
			#alpha_v_profile[idx]=np.mean(alpha_v[inbin])
		
		# Debugging output
		if debug==True:
			fig = plt.figure()
			ax = fig.add_subplot(111, projection='3d')
			ax.scatter(self.conf.rval[:,0], self.conf.rval[:,1], self.conf.rval[:,2], zdir='z', c='b')
			
		return [theta_out,rho_profile,vel_profile,eng_profile,press_profile,s_tt_profile,s_tp_profile,s_pt_profile,s_pp_profile,alpha_profile,alpha_v_profile,self.direction,self.orderpar]
        

class Tesselation:
    
	def __init__(self,conf):
		self.conf=conf
		self.rval=self.conf.rval
		
	def findLoop(self,dmax):
		neighList=[]
		self.Ival=[]
		self.Jval=[]
		Inei=[]
		count=0
		# Identify all neighbours and add them to a list. Keep i->j and j->i separate
		# The label is in neighList, the particle numbers are in Ival and Jval
		if self.conf.monodisperse:
			if self.conf.potential=='soft':
				dmax=4*self.sigma**2
			elif self.conf.potential=='morse':
				dmax=16*self.sigma**2
			else:
				dmax=4*self.sigma**2
				print "Warning: unimplemented potential, defaulting to maximum contact distance 2"
		if self.conf.potential=='morse':
			re=self.conf.param.pot_param['re']
		for i in range(len(self.rval)):
			#dist=np.sum((self.rval-self.rval[i,:])**2,axis=1)
			dist=self.conf.geom.GeodesicDistance(self.rval,self.rval[i,:])
			if self.conf.monodisperse:		
				neighbours=[index for index,value in enumerate(dist) if value <dmax]
			else:
				if self.conf.potential=='soft':
					neighbours=[index for index,value in enumerate(dist) if value < (self.conf.radius[i]+self.conf.radius[index])**2]
				elif self.conf.potential=='morse':	
					neighbours=[index for index,value in enumerate(dist) if value < (re*self.conf.radius[i]+re*self.conf.radius[index])**2]
				else:
					neighbours=[index for index,value in enumerate(dist) if value < (self.conf.radius[i]+self.conf.radius[index])**2]
					print "Warning: unimplemented potential, defaulting to maximum contact distance 2"
			neighbours.remove(i)
			neighList.extend([u for u in range(count,count+len(neighbours))])
			self.Ival.extend([i for k in range(len(neighbours))])
			self.Jval.extend(neighbours)
			Inei.append([u for u in range(count,count+len(neighbours))])
			count+=len(neighbours)
		# Identify loops based on the neighbour list. Kick out any (one-way) contacts that have occured so far
		Jarray=np.array(self.Jval)
		self.LoopList=[]
		# The dual: which loops belong to which particle
		self.ParList=[[] for k in range(len(self.rval))]
		self.LoopCen=[]
		l=0
		while len(neighList)>0:
			idx=neighList[0]
			idxkeep=idx
			#print idx
			idx0=[]
			#llist0=[]
			llist=[]
			goneround=False
			while goneround==False:  
				# Sort neighbours counterclockwise according to their local angle  
				dr0hat=self.rval[self.Jval[idx],:]-self.rval[self.Ival[idx],:]
				dr0hat/=np.sqrt(np.sum(dr0hat**2))
				jnei0=Inei[self.Jval[idx]]
				jnei=list(Jarray[jnei0])  
		
				drvec=self.rval[jnei,:]-self.rval[self.Jval[idx],:]
				drhat=((drvec).transpose()/(np.sqrt(np.sum(drvec**2,axis=1))).transpose()).transpose()
				cbeta=np.einsum('kj,j->k',drhat,self.conf.e2[self.Jval[idx],:])
				sbeta=np.einsum('kj,j->k',drhat,self.conf.e1[self.Jval[idx],:])
				cbeta0=np.dot(dr0hat,self.conf.e1[self.Jval[idx],:])
				sbeta0=np.dot(dr0hat,self.conf.e2[self.Jval[idx],:])
			
				# arccos returns between 0 and pi. Just multiply by the sign of the sine
				beta=np.arccos(cbeta)*np.sign(sbeta)
				# Determine the angles from the contact (read backwards) to the others, and pick the largest, modulo 2pi
				beta0=np.arccos(cbeta0)*np.sign(sbeta0)-np.pi
				dbeta=beta-beta0
				dbeta-=2*np.pi*np.round((dbeta-np.pi)/(2*np.pi))
				# and throwing out the particle itself
				itself=jnei.index(self.Ival[idx])
				dbeta[itself]=-1
				cnt=np.argmax(dbeta)
			
				idx=jnei0[cnt]
				goneround = idx in idx0
				if goneround==False:
					idx0.append(idx)
					llist.append(Jarray[idx])
					self.ParList[Jarray[idx]].append(l)
			#print idx0
			#print llist
			#print len(neighList)
			for v in idx0:
				try:
					neighList.remove(v)
				except ValueError:
					pass
			# There may be rare isolated cases (rattlers?) where the first contact itself is not part of the eventual loop.
			# This causes problems, because the loop identified after that has been removed.
			# Remove the original contact, in case it hasn't
			try:
				#print idxkeep
				neighList.remove(idxkeep)
			except ValueError:
				pass
			looppos=rval[llist]
			self.LoopCen.append([np.mean(looppos[:,0]), np.mean(looppos[:,1]),np.mean(looppos[:,2])])
			self.LoopList.append(llist)
			l+=1
		return self.LoopList,self.Ival,self.Jval
      
	# Much prettier: a loop that is too big (as measured by the mean square distance of the distances to the particles)
	# Deconstruct it into lots of little loops (virtual ones), with defined centers
	def makeEdges(self,rmax):   
		for l0 in range(len(self.LoopList)):
			llist=self.LoopList[l0]
			looppos=self.rval[llist]
			dlvec=looppos-self.LoopCen[l0]
			isLong=np.sqrt(np.sum(np.sum(dlvec**2,axis=1)))/len(llist)
			if len(llist)>5:
				print llist
				print isLong
			if isLong>rmax:
				print "Loop " + str(l0) + " with particles " + str(llist) + " is too big! "
				for k in range(len(llist)):
					kside=k-1
					if kside<0:
						kside=len(llist)-1
					# Attempting to catch the inward pointing loops: the have to be global boundary ~sqrt(N)
					if len(llist)<0.5*np.sqrt(len(self.rval)):
						newcen=0.5*(self.rval[llist[k]]+self.rval[llist[kside]])-conf.param.sigma*dlvec[k,:]/np.sqrt(np.sum(dlvec[k,:]**2))
					else:
						newcen=0.5*(self.rval[llist[k]]+self.rval[llist[kside]])+conf.param.sigma*dlvec[k,:]/np.sqrt(np.sum(dlvec[k,:]**2))
					self.LoopCen.append(newcen)
					try:
						ParList[llist[k]].remove(l0)
					except ValueError:
						pass
					self.ParList[llist[k]].append(l)
					try:
						ParList[llist[kside]].remove(l0)
					except ValueError:
						pass
					self.ParList[llist[kside]].append(l)
					l+=1
        
	# While we are at it, we can construct the dual tesselation here.
	# All that's missing is to order the patches for the particles counterclockwise
	def OrderPatches(self):
		LoopCen1=np.array(self.LoopCen)
		for i in range(len(self.rval)):
			parray=np.array(self.ParList[i])
			drvec=LoopCen1[self.ParList[i]]-self.rval[i,:]
			# Optionally Take care of irregularities (in the form of too long bonds) here. These happen at the edges of connected stuff
			# The tesselation is correct, it's just not what we want
			drlen=np.sqrt(np.sum(drvec**2,axis=1))
			#drvec=rval[jnei,:]-rval[Jval[idx],:]
			drhat=((drvec).transpose()/(drlen).transpose()).transpose()
			cbeta=np.einsum('kj,j->k',drhat,self.conf.e1[i,:])
			sbeta=np.einsum('kj,j->k',drhat,self.conf.e2[i,:])
			# arccos returns between 0 and pi. Just multiply by the sign of the sine
			beta=np.arccos(cbeta)*np.sign(sbeta)
			# sort by angle and put back in ParList
			lorder=np.argsort(beta)
			ParList[i]=parray[lorder] 
			# Use the new ParList structure where loops belong to particles are stored
		return self.LoopList,self.LoopCen,self.ParList,self.Ival,self.Jval
    
class Defects:
	def __init__(self,tess,conf):
		self.LoopList=tess.LoopList
		self.conf=conf
		
		self.numdefect_n=0
		self.numdefect_v=0
		# Defect storage, up to 100
		# For n and velocity
		self.defects_n=np.zeros((100,4))
		self.defects_v=np.zeros((100,4))
		self.printnow=False
        
	def getDefects(self,symtype): 
		# Generalized algorithm for defects of any type
		# Count the defect charge. Times two, to use integers and easier if statements
		for u in range(len(self.LoopList)):
			thisLoop=self.LoopList[u]
			if symtype=='oldnematic':
				ndefect,vdefect=self.getDefectsGoldenfeld(thisLoop)
			elif symtype=='polar':
				ndefect,vdefect=self.getDefectsPolar(thisLoop)
			elif symtype=='nematic':
				ndefect,vdefect=self.getDefectsNematic(thisLoop)
			else:
				print "Unknown alignment symmetry type! Not tracking defects!"
				ndefect=0.0
				vdefect=0.0
			if abs(ndefect)>0:
				if self.numdefect_n<100:
					print "Found Defect in orientation field!"
					print ndefect
					# Construct the geometric centre of the defect
					rmhat=np.sum(self.conf.rval[thisLoop],axis=0)
					rmhat/=np.sqrt(np.sum(rmhat**2))
					# Charge of the defect
					self.defects_n[self.numdefect_n,0]=ndefect
					# Coordinates of the defect
					self.defects_n[self.numdefect_n,1:]=radius*rmhat
					self.numdefect_n+=1
			if abs(vdefect)>0:
				if self.numdefect_v<100:
					print "Found Defect in velocity field!"
					print vdefect
					# Construct the geometric centre of the defect
					rmhat=np.sum(self.conf.rval[thisLoop],axis=0)
					rmhat/=np.sqrt(np.sum(rmhat**2))
					# Charge of the defect
					self.defects_v[self.numdefect_v,0]=vdefect
					# Coordinates of the defect
					self.defects_v[self.numdefect_v,1:]=radius*rmhat
					self.numdefect_v+=1

		#print defects
		print 'Number of orientation field defects: ' + str(self.numdefect_n)
		print 'Number of velocity field defects: ' + str(self.numdefect_v)
		return self.defects_n, self.defects_v,self.numdefect_n,self.numdefect_v
        
	def getDefectsNematic(self,thisLoop): 
		# Generalized algorithm for defects of any type
		# Count the defect charge. Times two, to use integers and easier if statements
		# nval
		thetatot=0
		t0=thisloop[0]
		ctheta=1
		for t in thisLoop[1:-1]:
			ctheta=np.dot(self.conf.nval[t,:],np.sign(ctheta)*self.conf.nval[t0,:])
			stheta=np.dot(self.conf.rhat[t,:],np.cross(self.conf.nval[t,:],self.conf.nval[t0,:]))
			theta=np.arccos(ctheta)*np.sign(stheta)
			thetatot+=theta
			t0=t
		ndefect=0.5*int(round(thetatot/(np.pi)))
		# vhat
		thetatot=0
		t0=thisloop[0]
		ctheta=1
		for t in thisLoop[1:-1]:
			ctheta=np.dot(self.conf.vhat[t,:],np.sign(ctheta)*self.conf.vhat[t0,:])
			stheta=np.dot(self.conf.rhat[t,:],np.cross(self.conf.nval[t,:],self.conf.vhat[t0,:]))
			theta=np.arccos(ctheta)*np.sign(stheta)
			thetatot+=theta
			t0=t
		vdefect=0.5*int(round(thetatot/(np.pi)))
		return ndefect,vdefect
			
            
	def getDefectsGoldenfeld(self,thisLoop): 
		# Should already be ordered counterclockwise
		# Following a version of the Goldenfeld algorithm, with nx,ny,nz as is playing the role of the order parameter. The sphere is in cartesian space
		# The old nematic algorithm, based on the hemispheres            
		# The polarization vector nval
		ctheta=1
		coord=[]
		coord.append(nval[thisLoop[0],:])
		for t in range(1,len(thisLoop)):
			ctheta=np.dot(self.conf.nval[thisLoop[t],:],np.sign(ctheta)*self.conf.nval[thisLoop[t-1],:])
			# Nematic: append the order parameter, rotated through the *smaller* angle
			coord.append(np.sign(ctheta)*self.conf.nval[thisLoop[t],:])
			# Find out if the last point and the starting point are in the same hemisphere. 
		cdefect=np.dot(coord[t],coord[0])
		if cdefect<0:
			ndefect=0.5
		else:
			ndefect=0.0
		# The normalized velocity vector vhat
		ctheta=1
		coord=[]
		coord.append(vhat[thisLoop[0],:])
		for t in range(1,len(thisLoop)):
			ctheta=np.dot(self.conf.vhat[thisLoop[t],:],np.sign(ctheta)*self.conf.vhat[thisLoop[t-1],:])
			# Nematic: append the order parameter, rotated through the *smaller* angle
			coord.append(np.sign(ctheta)*self.conf.vhat[thisLoop[t],:])
			# Find out if the last point and the starting point are in the same hemisphere. 
		cdefect=np.dot(coord[t],coord[0])
		if cdefect<0:
			vdefect=0.5
		else:
			vdefect=0.0
		return ndefect,vdefect
       
	def getDefectsPolar(self,thisLoop):
		# Generalized algorithm for defects of any type
		# Count the defect charge. 
		# Should already be ordered counterclockwise
		# nval
		thetatot=0
		t0=thisLoop[-1]
		for t in thisLoop[0:len(thisLoop)]:
			ctheta=np.dot(self.conf.nval[t,:],self.conf.nval[t0,:])    
			stheta=np.dot(self.conf.rhat[t,:],np.cross(self.conf.nval[t,:],self.conf.nval[t0,:]))
			theta=np.arccos(ctheta)*np.sign(stheta)
			thetatot+=theta
			t0=t
		# Classify according to defects
		# For a polar one, we can only have integer defects
		ndefect=int(round(thetatot/(2*np.pi)))
		# vhat
		thetatot=0
		t0=thisLoop[-1]
		for t in thisLoop[0:len(thisLoop)]:
			ctheta=np.dot(self.conf.vhat[t,:],self.conf.vhat[t0,:])    
			stheta=np.dot(self.conf.rhat[t,:],np.cross(self.conf.vhat[t,:],self.conf.vhat[t0,:]))
			theta=np.arccos(ctheta)*np.sign(stheta)
			thetatot+=theta
			t0=t
			#if ctheta<0:
				#print "candidate: t t0 ctheta stheta theta thetatot"
				#print t, t0, ctheta, stheta, theta, thetatot
				#printnow=True
		# Classify according to defects
		# For a polar one, we can only have integer defects
		vdefect=int(round(thetatot/(2*np.pi)))
		#if printnow:
			#print thetatot
			#print thisLoop
		return ndefect,vdefect            

	def PlotDefects(self):
		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')
		ax.scatter(self.conf.rval[:,0], self.conf.rval[:,1], self.conf.rval[:,2], zdir='z', c='b',s=4)
		ax.scatter(defects_n[:,1], defects_n[:,2], defects_n[:,3], zdir='z', c='r',s=50)
		ax.scatter(defects_v[:,1], defects_v[:,2], defects_v[:,3], zdir='z', c='g',s=50)
            
  
class Writer:
	def __init__(self,nematic=False,connected=False):
		self.nematic=nematic
		self.connected=connected
		
	def writeConfigurationVTK(conf,outfile):
		# Data which goes into file: positions, directors, velocities
		# radii
		r = conf.radius
		# positions
		x = conf.rval[:,0]
		y = conf.rval[:,1]
		z = conf.rval[:,2]
		# directors
		nx = conf.nval[:,0]
		ny = conf.nval[:,1]
		nz = conf.nval[:,2]
		# velocities
		vx = conf.vval[:,0]
		vy = conf.vval[:,1]
		vz = conf.vval[:,2]
		
		# Preparking the vtk structure
		Points = vtk.vtkPoints()
		
		Radii = vtk.vtkDoubleArray()
		Radii.SetNumberOfComponents(1)
		Radii.SetName('Radius')

		Velocities = vtk.vtkDoubleArray()
		Velocities.SetNumberOfComponents(3)
		Velocities.SetName("Velocity")

		Directors = vtk.vtkDoubleArray()
		Directors.SetNumberOfComponents(3)
		Directors.SetName("Directors")
		
		if self.nematic:
			NDirectors = vtk.vtkDoubleArray()
			NDirectors.SetNumberOfComponents(3)
			NDirectors.SetName("NDirectors")
		
		# Adding the data to the vtk structures
		for (xx,yy,zz,rr) in zip(x,y,z,r):
			Points.InsertNextPoint(xx,yy,zz)
			Radii.InsertNextValue(rr)
		for (vvx,vvy,vvz) in zip(vx,vy,vz):
			Velocities.InsertNextTuple3(vvx,vvy,vvz)
		for (nnx,nny,nnz) in zip(nx,ny,nz):	
			if self.nematic:
				Directors.InsertNextTuple3(0.5*nnx,0.5*nny,0.5*nnz)
				NDirectors.InsertNextTuple3(-0.5*nnx,-0.5*nny,-0.5*nnz)
			else:
				Directors.InsertNextTuple3(nnx,nny,nnz)
		# Connected, using convex hull (? ask Rastko ...?)
		if self.connected:
			Lines = vtk.vtkCellArray()
			Line = vtk.vtkLine()
			points = np.column_stack((x,y,z)) 
			hull = ConvexHull(points)
			edges = []
			for h in hull.simplices:
			i, j, k = h
			if not sorted([i,j]) in edges: edges.append(sorted([i,j]))
			if not sorted([i,k]) in edges: edges.append(sorted([i,k]))
			if not sorted([j,k]) in edges: edges.append(sorted([j,k]))
			for (i,j) in edges:
			Line.GetPointIds().SetId(0,i)
			Line.GetPointIds().SetId(1,j)
			Lines.InsertNextCell(Line)
			
		# Putting the results into a polydata structure
		polydata = vtk.vtkPolyData()
		polydata.SetPoints(Points)
		if self.connected:
			polydata.SetLines(Lines)
		polydata.GetPointData().AddArray(Radii)
		polydata.GetPointData().AddArray(Velocities)
		polydata.GetPointData().AddArray(Directors)
		if self.nematic:
			polydata.GetPointData().AddArray(NDirectors)
		polydata.Modified()
		
		# Finally, output via binary writer
		writer = vtk.vtkXMLPolyDataWriter()
		#outname = '.'.join(f.split('.')[:-1])
		writer.SetFileName(outfile)
		if vtk.VTK_MAJOR_VERSION <= 5:
			writer.SetInput(polydata)
		else:
			writer.SetInputData(polydata)
		writer.SetDataModeToAscii()
		writer.Write()
		
	def writeDefects(self,defects_n, defects_v,numdefect_n,numdefect_v,outfile):
		# Preparing the vtp output
		# Create point structure in vtk
		Points = vtk.vtkPoints()
		print "Created Points"
		# Create (something) associated to the points, with different values for each
		Number = vtk.vtkDoubleArray()
		Number.SetNumberOfComponents(1)
		Number.SetName('Number')
		Size = vtk.vtkDoubleArray()
		Size.SetNumberOfComponents(1)
		Size.SetName('Size')
		print "Created Number"
		# Put one point at the centre, and the ndefect ones around it
		Points.InsertNextPoint(0,0,0)
		Number.InsertNextValue(0)
		Size.InsertNextValue(0)
		for u in range(numdefect_n):
			Points.InsertNextPoint(defects_n[u,1],defects_n[u,2],defects_n[u,3])
			Number.InsertNextValue(1)
			Size.InsertNextValue(1.0)
		for u in range(numdefect_v):
			Points.InsertNextPoint(defects_v[u,1],defects_v[u,2],defects_v[u,3])
			Number.InsertNextValue(2)
			Size.InsertNextValue(1.0)
		print "Added Particles and Numbers"
		
		lines = vtk.vtkCellArray()
		line = vtk.vtkLine()
		for i in range(numdefect_n):
			line = vtk.vtkLine()
			line.GetPointIds().SetId(0,0)
			line.GetPointIds().SetId(1,i+1)
			lines.InsertNextCell(line)
		for i in range(numdefect_v):
			line = vtk.vtkLine()
			line.GetPointIds().SetId(0,0)
			line.GetPointIds().SetId(1,numdefect_n+i+1)
			lines.InsertNextCell(line)
		print "Added lines"
		
		polydata = vtk.vtkPolyData()
		polydata.SetPoints(Points)
		polydata.SetLines(lines)
		polydata.GetPointData().AddArray(Number)
		polydata.GetPointData().AddArray(Size)
		print "Finished Polydata"
		polydata.Modified()
		writer = vtk.vtkXMLPolyDataWriter()
		writer.SetFileName(outfile)
		# Python 2.7 vs. 3 incompatibility?
		if vtk.VTK_MAJOR_VERSION <= 5:
			writer.SetInput(polydata)
		else:
			writer.SetInputData(polydata)
		writer.SetDataModeToAscii()
		writer.Write()
		print "Wrote File"
		
	def writePatches(self,tess,outname):
		print outname
		points = vtk.vtkPoints()
		polygons = vtk.vtkCellArray()
		v=0
		polygon = vtk.vtkPolygon()
		havePoly=[]
		for k in range(len(tess.ParList)):
			nedge=len(tess.ParList[k])
			if nedge<2:
			print nedge
			print k
			print tess.ParList[k]
			else:
			havePoly.append(k)
			#for k in range(300):
			# Create the points of the polygon: the loop centers
			polygon = vtk.vtkPolygon()
			for l in tess.ParList[k]:
				points.InsertNextPoint(tess.LoopCen[l][0],tess.LoopCen[l][1],tess.LoopCen[l][2])
			polygon.GetPointIds().SetNumberOfIds(nedge)
			for l in range(nedge):
				#print l
				polygon.GetPointIds().SetId(l,v+l)
			
			polygons.InsertNextCell(polygon)
			v+=nedge
		# Create the matching polydata 
		polygonPolyData = vtk.vtkPolyData()
		polygonPolyData.SetPoints(points)
		polygonPolyData.SetPolys(polygons)
		# Add stresses ...
		eng, press,stress = tess.conf.compute_energy_and_pressure()
		pressure = vtk.vtkDoubleArray()
		pressure.SetNumberOfComponents(1)
		pressure.SetName('Pressure')
		for k in havePoly:
			pressure.InsertNextValue(press[k])
		polygonPolyData.GetCellData().AddArray(pressure)
			
		writer = vtk.vtkXMLPolyDataWriter()
		writer.SetFileName(outname)
		if vtk.VTK_MAJOR_VERSION <= 5:
			writer.SetInput(polygonPolyData)
		else:
			writer.SetInputData(polygonPolyData)
		writer.SetDataModeToAscii()
		writer.Write()	
  
    

      


  

  
# Scripting version: Only execute if this is called as a script. Otherwise, it attempts to go through here when loading as a module 
# and throws errors because some arguments aren't defined
if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--input", type=str, help="Input file with particle velocity field")
  parser.add_argument("-o", "--output", type=str, default="defects", help="Output file (text file)")
  parser.add_argument("-k", "--k", type=float, default=1.0, help="soft potential strength")
  parser.add_argument("-R", "--sphere_r", type=float, default=28.2094791, help="radius of sphere for spherical system")
  parser.add_argument("-r", "--particle_r", type=float, default=1.0, help="radius of particle ")
  args = parser.parse_args()

  print
  print "\tActive Particles on Curved Spaces (APCS)"
  print "\tPolar and nematic defect finding algoritm"
  print 
  print "\tSilke Henkes"
  print "\tUniversity of Aberdeen"
  print "\t(c) 2014"
  print "\t----------------------------------------------"
  print 
  print "\tInput : ", args.input
  print "\tOutput : ", args.output
  print "\tSpring constant : ", args.k
  print "\tRadius of the sphere : ", args.sphere_r
  print "\tRadius of the particle : ", args.particle_r
  print 

  outname = '.'.join((args.input).split('.')[:-1]) + '_data.vtp'
  print outname
  
  defects_n, defects_v,numdefect_n,numdefect_v=getDefects(args.input,args.sphere_r,args.particle_r,outname,'polar',True,True)
  
  outname = '.'.join((args.input).split('.')[:-1]) + '_defects.vtp'
  print outname
  #writer.SetFileName(args.output+'/'+outname+'.vtp')
  #writer.SetFileName(args.output+'.vtp')
  
  writeDefects(defects_n, defects_v,numdefect_n,numdefect_v,outname)
  
  
  

  plt.show()