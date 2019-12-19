# ***************************************************************************
# *
# *  Copyright (C) 2013-2016 University of Dundee
# *  All rights reserved. 
# *
# *  This file is part of SAMoS (Soft Active Matter on Surfaces) program.
# *
# *  SAMoS is free software; you can redistribute it and/or modify
# *  it under the terms of the GNU General Public License as published by
# *  the Free Software Foundation; either version 2 of the License, or
# *  (at your option) any later version.
# *
# *  SAMoS is distributed in the hope that it will be useful,
# *  but WITHOUT ANY WARRANTY; without even the implied warranty of
# *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# *  GNU General Public License for more details.
# *
# *  You should have received a copy of the GNU General Public License
# *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
# *
# *****************************************************************************

from numpy import linalg as LA
from Geometry import *
from Configuration import *
import glob
try:
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib.cm as cm
	HAS_MATPLOTLIB = True
except:
	HAS_MATPLOTLIB = False
	pass

# Profiles for the cornea case. 
# Warning: Implicitly spherical geometry. Results in other geometries will be garbage (should not crash though)
class Cornea:
	def __init__(self,directory,conffile,skip,howmany,ignore=True,maxtype=3):
		self.ignore=ignore
		self.param = Param(directory+conffile)
		files = sorted(glob.glob(directory + self.param.dumpname+'*.dat'))[skip:(skip+howmany)]
		if len(files) == 0:
  			files = sorted(glob.glob(directory + self.param.dumpname+'*.dat.gz'))[skip:(skip+howmany)]
		# Read the local data
		geometries={'sphere':GeometrySphere,'plane':GeometryPlane,'plane_periodic':GeometryPeriodicPlane,'none':Geometry,'tube':GeometryTube,'peanut':GeometryPeanut,'hourglass':GeometryHourglass}
		# Create the right geometry environment (TBC):
		self.geom=geometries[self.param.constraint](self.param)
		print self.geom
		# Number of files that we are dealing with
		self.Nsnap=len(files)
		print self.Nsnap
		# Find out if this is a simulation with a variable number of particles:
		# For now, otherwise we get a mess ...
		self.Nvariable=True
		# Deal with the radii in an appropriate manner:
		# First check that the potential actually uses radii:
		# If yes, read them in from the initial file (to be overwritten if there are separate ones in each file)
		if self.ignore:
			self.monodisperse=True
		else:
			if self.param.pot_params['use_particle_radii']==True:
				self.monodisperse=False
			else:
				self.monodisperse=True
		# Unfortunately, need a first read-through to gauge what kind of data size we need
		# I want to use numpy arrays for the analysis due to speed
		self.Nval=np.zeros((self.Nsnap,)).astype(int)
		u=0
		for f in files:
			#print "Pre - Processing file : ", f
			data = ReadData(f)
			tp=np.array(data.data[data.keys['type']])
			self.Nval[u]=len(np.where(tp<maxtype)[0])
			#tracers = [index for index,value in enumerate(tp) if value==2]
			u+=1
		print self.Nval
		self.N=int(np.amax(self.Nval))
		print "Handling a total of maximum " + str(self.N) + " particles!"
		# Produce empty arrays for initialization
		self.rval=np.zeros((self.Nsnap,self.N,3))
		self.vval=np.zeros((self.Nsnap,self.N,3))
		self.nval=np.zeros((self.Nsnap,self.N,3))
		self.flag=np.zeros((self.Nsnap,self.N)).astype(int)
		self.radius=np.zeros((self.Nsnap,self.N))
		u=0
		for f in files:
			print "Processing file : ", f
			data = ReadData(f)
			tp=np.array(data.data[data.keys['type']])
			# use only TA cells for analysis
			useparts = np.where(tp<maxtype)[0]
			# First read in everything ....
			x, y, z = np.array(data.data[data.keys['x']]), np.array(data.data[data.keys['y']]), np.array(data.data[data.keys['z']])
			vx, vy, vz = np.array(data.data[data.keys['vx']]), np.array(data.data[data.keys['vy']]), np.array(data.data[data.keys['vz']])
			nx, ny, nz = np.array(data.data[data.keys['nx']]), np.array(data.data[data.keys['ny']]), np.array(data.data[data.keys['nz']])	
			fl = np.array(data.data[data.keys['flag']])
			rad=np.array(data.data[data.keys['radius']])
			# and then keep what we would like to analyse
			self.radius[u,:self.Nval[u]]=rad[useparts]
			self.flag[u,:self.Nval[u]]=fl[useparts]
			
			rval_u = np.column_stack((x,y,z))
			vval_u = np.column_stack((vx,vy,vz))
			nval_u = np.column_stack((nx,ny,nz))
			self.rval[u,:self.Nval[u],:]=rval_u[useparts,:]
			self.vval[u,:self.Nval[u],:]=vval_u[useparts,:]
			self.nval[u,:self.Nval[u],:]=nval_u[useparts,:]
			
			u+=1
			
	def getFlowField(self,inisnap,dsnap=1,debug=True):
		# attempt to compute the flow field between two snapshots, based on the uniquely labeled particles present in both
		flag1=list(self.flag[inisnap,:self.Nval[inisnap]])
		flag2=list(self.flag[inisnap+dsnap,:self.Nval[inisnap+dsnap]])
		# intersection: do this manually, possibly easier ...
		# labels in order for both
		label1=[]
		label2=[]
		index=[]
		hasdied=[]
		for k1 in range(len(flag1)):
			try:
				k2=flag2.index(flag1[k1])
				label1.append(k1)
				label2.append(k2)
				index.append(flag1[k1])
			except:
				print "particle " + str(flag1[k1]) + " died."
				hasdied.append(k1)
		# now compute the flow field from the difference in position
		# time that has passed
		self.param.dt=0.001 # eff it for now, it's always that ...
		deltat = self.param.dt*self.param.dump['freq']*dsnap
		flow_field = (self.rval[inisnap+dsnap,label2,:]-self.rval[inisnap,label1,:])/deltat
		## just eff it, close the holes in there with an average...
		#flow_field = np.zeros((self.Nval[inisnap],3))
		#flow_field[label1,:]=flow_field0
		## missing buggers
		#for k in hasdied:
			#neighbours=self.conf.getNeighbours(i,mult,dmax)[0]
		
		if debug:
			plt.figure(figsize=(8,8))
			velangle = np.arctan2(flow_field[:,1],flow_field[:,0])
			#arctan2 gives results between -pi to pi, normalise ...
			colors = (velangle+np.pi)/(2*np.pi)
			colormap = cm.hsv
			plt.quiver(self.rval[inisnap,label1,0],self.rval[inisnap,label1,1],flow_field[:,0],flow_field[:,1],color=0.8*colormap(colors))
			plt.gca().set_aspect('equal')
			plt.title('Skipping ' + str(dsnap))
			
			plt.figure(figsize=(8,8))
			velangle = np.arctan2(self.vval[inisnap,:,1],self.vval[inisnap,:,0])
			#arctan2 gives results between -pi to pi, normalise ...
			colors = (velangle+np.pi)/(2*np.pi)
			colormap = cm.hsv
			plt.quiver(self.rval[inisnap,label1,0],self.rval[inisnap,label1,1],self.vval[inisnap,:,0],self.vval[inisnap,:,1],color=0.8*colormap(colors))
			plt.gca().set_aspect('equal')
			plt.title('Instantaneous velocity field')
			
		return flow_field, label1
	
		
	#def __init__(self,conf,geom,nematic=False,debug=False):
		#self.conf=conf
		#self.geom=geom
		##self.rval=self.conf.rval
		#ez = np.array([0,0,1])  # lab frame z-axis
		## Simply get the axis as the mean crossproduct or r and v; assuming alignment. This should also not flip.
		#if not nematic:
			#self.direction=np.sum(np.cross(conf.rval,conf.vval),axis=0)
		## Otherwise we can't do this since it's going to be close to 0
		#else:
			#print "doing nematic case!"
			## Be a bit better here. The moment of inertia tensor should have *some* kind of signature
			## Take this as an initial guess of the axis. Compare the two later.
			#Itensor=np.einsum('ik,ij->kj',conf.rval,conf.rval)
			##print Itensor
			## Compute its eigenvalues and eigenvectors
			#eigval,eigvec=LA.eig(Itensor)
			##print eigval
			## Careful, eigenvectors are defined this way round
			##print eigvec[:,0]
			##print eigvec[:,1]
			##print eigvec[:,2]
			## The eigenvector we want is the one with the *smallest* eigenvalue, as our shape
			## is a squashed thing with the z axis in the squattest direction
			#idx=np.argmin(eigval)
			#inertialz=eigvec[:,idx]
			#print inertialz 
			
			#directions=np.cross(conf.rval,conf.vval)
			## Those should then now be mostly either aligned or antialigned
			## Hope for the best and rectify them ...
			## Use the guess from the inertia tensor now
			#normdir=np.sqrt(directions[:,0]**2+directions[:,1]**2+directions[:,2]**2)
			#dirnorm=((directions).transpose()/(normdir).transpose()).transpose()
			#orient=np.round(dirnorm[:,0]*inertialz[0]+dirnorm[:,1]*inertialz[1]+dirnorm[:,2]*inertialz[2])
			#print orient
			#print sum(orient)
			#print sum(orient**2)/len(conf.rval)
			#self.direction=np.einsum('ij,i->j',directions,orient)
			##self.direction=np.empty((len(conf.vval),3))
			##self.direction[:,0]=directions[:,0]*orient
			##self.direction[:,1]=directions[:,1]*orient
			##self.direction[:,2]=directions[:,2]*orient
		#self.orderpar=self.direction/len(conf.rval)
		#print self.orderpar
		##self.direction = self.direction/np.linalg.norm(self.direction)
		#self.direction=inertialz
		#print self.direction
		## to plot, make this a line ...
		
		#axis = np.cross(self.direction,ez)
		#axis = axis/np.linalg.norm(axis)
		#rot_angle = np.arccos(np.dot(self.direction,ez))
		#print rot_angle
		#axis0 = np.empty(np.shape(self.conf.rval))
		#axis0[:,0] = axis[0]
		#axis0[:,1] = axis[1]
		#axis0[:,2] = axis[2]
		#if debug:
			## Debugging output
			#if HAS_MATPLOTLIB:
				#fig = plt.figure()
				#ax = fig.add_subplot(111, projection='3d')
				#ax.scatter(self.conf.rval[:,0], self.conf.rval[:,1], self.conf.rval[:,2], zdir='z', c='b')
				#ax.plot([-self.geom.R*self.direction[0],self.geom.R*self.direction[0]],[-self.geom.R*self.direction[1],self.geom.R*self.direction[1]],[-self.geom.R*self.direction[2],self.geom.R*self.direction[2]],'o-r')
				#ax.plot([0,0],[0,0],[-self.geom.R,self.geom.R],'o-g')
				#ax.plot([-self.geom.R*axis[0],self.geom.R*axis[0]],[-self.geom.R*axis[1],self.geom.R*axis[1]],[-self.geom.R*axis[2],self.geom.R*axis[2]],'o-k')
			#else:
				#print 'Error: Matplotlib does not exist on this machine, cannot plot system'
		#self.conf.rotateFrame(axis0,rot_angle)
		## Need to redo the cell list after the rotation
		#self.conf.redoCellList()
		
		##rhat=((rval).transpose()/(np.sqrt(np.sum(rval**2,axis=1))).transpose()).transpose()
		#vel = np.sqrt(self.conf.vval[:,0]**2 + self.conf.vval[:,1]**2 + self.conf.vval[:,2]**2)
		#velnorm=((self.conf.vval).transpose()/(vel).transpose()).transpose()
  
		#self.theta,self.phi,self.etheta,self.ephi=geom.TangentBundle(self.conf.rval)
		## Alpha, the angle between the local polarity and the equator; here represented by ephi
		#self.alpha=-np.arcsin(np.sum(self.conf.nval*self.etheta, axis=1))
		## Same thing for the velocity
		## No - add pi/2 to get something that does not add up to zero 
		##alpha_v=np.arccos(np.sum(velnorm*etheta, axis=1))
  
	def getProfiles(self,nbin,debug=False):
		eng, press,ncon,stress = self.conf.compute_energy_and_pressure()
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
		normz=2*np.pi*self.geom.R*abs(np.cos(theta_out))
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
			if HAS_MATPLOTLIB:
				fig = plt.figure()
				ax = fig.add_subplot(111, projection='3d')
				ax.scatter(self.conf.rval[:,0], self.conf.rval[:,1], self.conf.rval[:,2], zdir='z', c='b')
			else:
				print 'Error: Matplotlib does not exist on this machine, cannot plot system'
			
		return [theta_out,rho_profile,vel_profile,eng_profile,press_profile,s_tt_profile,s_tp_profile,s_pt_profile,s_pp_profile,alpha_profile,alpha_v_profile,self.direction,self.orderpar]
        
