# * *************************************************************
# *  
# *   Soft Active Mater on Surfaces (SAMoS)
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

import sys
import argparse
import pickle

from Writer import *
from Geometry import *
from Hessian import *
from read_param import *
from read_data import *
try:
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
	HAS_MATPLOTLIB=True
	
	import matplotlib
	matplotlib.rcParams['text.usetex'] = 'false'
	matplotlib.rcParams['lines.linewidth'] = 2
	matplotlib.rcParams['axes.linewidth'] = 2
	matplotlib.rcParams['xtick.major.size'] = 8
	matplotlib.rcParams['ytick.major.size'] = 8
	matplotlib.rcParams['font.size']=16.0
	matplotlib.rcParams['legend.fontsize']=14.0

	cdict = {'red':   [(0.0,  0.0, 0.5),
					  (0.35,  1.0, 0.75),
					  (0.45,  0.75, 0.0),
					  (1.0,  0.0, 0.0)],

			'green': [(0.0,  0.0, 0.0),
					  (0.35,  0.0, 0.5),
					  (0.5, 1.0, 1.0),
					  (0.8,  0.5, 0.0),
					  (1.0,  0.0, 0.0)],

			'blue':  [(0.0,  0.0, 0.0),
					  (0.5,  0.0, 0.0),
					  (0.7, 0.5, 1.0),
					  (1.0,  0.25, 0.0)]}
except:
	HAS_MATPLOTLIB=False
	pass

class SimRun:
	def __init__(self,directory,conffile,inputfile,radiusfile,skip,ignore=False,tracer=False,debug=False):
		self.debug=debug
		self.tracer=tracer
		self.ignore=ignore
		self.param = Param(directory+conffile)
		files = sorted(glob(directory + self.param.dumpname+'*.dat'))[skip:]
		if len(files) == 0:
  			files = sorted(glob(directory + self.param.dumpname+'*.dat.gz'))[skip:]
		# Read the local data
		geometries={'sphere':GeometrySphere,'plane':GeometryPeriodicPlane,'none':Geometry,'tube':GeometryTube,'peanut':GeometryPeanut,'hourglass':GeometryHourglass}
		# Create the right geometry environment (TBC):
		self.geom=geometries[self.param.constraint](self.param)
		print self.geom
		# Number of files that we are dealing with
		self.Nsnap=len(files)
		# Find out if this is a simulation with a variable number of particles:
		if self.param.npopulation>0:
			self.Nvariable=True
		else:
			self.Nvariable=False
		# Deal with the radii in an appropriate manner:
		# First check that the potential actually uses radii:
		# If yes, read them in from the initial file (to be overwritten if there are separate ones in each file)
		if self.ignore:
			self.monodisperse=True
		else:
			if self.param.pot_params['use_particle_radii']==True:
				print "Reading radii from initial file!"
				data_ini=ReadData(directory+radiusfile)
				self.monodisperse=False
				self.radius=data_ini.data[data_ini.keys['radius']]
				if self.Nvariable==False:
					self.N=len(self.radius)
					print "Constant number of " + str(self.N) + " particles!"
			else:
				self.monodisperse=True
		# Unfortunately, need a first read-through to gauge what kind of data size we need
		# I want to use numpy arrays for the analysis due to speed
		if self.Nvariable:
			self.Nval=np.zeros((self.Nsnap,))
			u=0
			for f in files:
				 #print "Pre - Processing file : ", f
				 data = ReadData(f)
				 x= np.array(data.data[data.keys['x']])
				 self.Nval[u]=len(x)
				 #print self.Nval[u]	
				 u+=1
			self.N=int(np.amax(self.Nval))
			print "Handling a total of maximum " + str(self.N) + " particles!"
		elif self.monodisperse:
			data = ReadData(files[0])
			x= np.array(data.data[data.keys['x']])
			self.N=len(x)
			print "Constant number of " + str(self.N) + " particles!"
		if tracer:
			data = ReadData(files[0])
			tp=np.array(data.data[data.keys['type']])
			tracers = [index for index,value in enumerate(tp) if value==2]
			print tracers
			self.Ntracer=len(tracers)
			print "Tracking " + str(self.Ntracer) + " tracer particles!"
		# Produce empty arrays for initialization
		self.rval=np.zeros((self.Nsnap,self.N,3))
		self.vval=np.zeros((self.Nsnap,self.N,3))
		#self.nval=np.empty((self.Nsnap,self.N,3))
		#self.vhat=np.empty((self.Nsnap,self.N,3))
		if self.Nvariable:
			self.flag=np.zeros((self.Nsnap,self.N))
			if not self.monodisperse:
				 self.radius=np.zeros((self.Nsnap,self.N))
		if tracer:
			self.tracers=np.zeros((self.Nsnap,self.Ntracer))
			self.rtracers=np.zeros((self.Nsnap,self.Ntracer,3))
			self.vtracers=np.zeros((self.Nsnap,self.Ntracer,3))
		u=0
		for f in files:
			print "Processing file : ", f
			data = ReadData(f)
			x, y, z = np.array(data.data[data.keys['x']]), np.array(data.data[data.keys['y']]), np.array(data.data[data.keys['z']])
			vx, vy, vz = np.array(data.data[data.keys['vx']]), np.array(data.data[data.keys['vy']]), np.array(data.data[data.keys['vz']])
			#if data.keys.has_key('nx'):
			#	nx, ny, nz = np.array(data.data[data.keys['nx']]), np.array(data.data[data.keys['ny']]), np.array(data.data[data.keys['nz']])	
			if data.keys.has_key('flag'):
				fl = data.data[data.keys['flag']]
				if self.Nvariable:
					self.flag[u,:self.Nval[u]]=fl
				else:
					self.flag[u,:]=fl
			if data.keys.has_key('type'):
				tp = data.data[data.keys['type']]
				if tracer:
					tracers = [index for index,value in enumerate(tp) if value==2]
					self.tracers[u,:]=tracers
			if data.keys.has_key('radius'):
				rad=data.data[data.keys['radius']]
				if self.Nvariable:
					self.radius[u,:self.Nval[u]]=rad
				#else:
					#self.radius[u,:]=rad
			#if u>0 and takeDrift:
				#rval_u0=rval_u
				#print u
			rval_u = np.column_stack((x,y,z))
			vval_u = np.column_stack((vx,vy,vz))
			#nval_u = np.column_stack((nx,ny,nz))
			#vel = np.sqrt(vval_u[:,0]**2 + vval_u[:,1]**2 + vval_u[:,2]**2)
			#vhat_u=((vval_u).transpose()/(vel).transpose()).transpose()
			if self.Nvariable:
				self.rval[u,:self.Nval[u],:]=rval_u
				self.vval[u,:self.Nval[u],:]=vval_u
				if tracer:
					#self.rtracers[u,:]=self.rval[u,self.tracers[u,:]]
					#self.vtracers[u,:]=self.rval[u,self.tracers[u,:]]
					self.rtracers[u,:,:]=self.rval[u,tracers,:]
					self.vtracers[u,:,:]=self.rval[u,tracers,:]
			else:
				self.rval[u,:,:]=rval_u
				self.vval[u,:,:]=vval_u
			u+=1
			
	#def takeDrift():
		## taking off the drift
		#if u>0:
			## Bring back people who crossed the line
			#diff_real =self.geom.ApplyPeriodic12(rval_u0,rval_u)
			##print np.max(diff_real)
			##print diff_real
			#rval_ud=rval_u0+diff_real
			#drift=np.sum(diff_real,axis=0)/self.N
			##print drift
			#rval_ud-=drift
			#self.rval[u,:,:]=rval_ud
		#else:
			#self.rval[u,:,:]=rval_u
		
	def getMSD(self):
		self.msd=np.empty((self.Nsnap,))
		# in case of tracers, simply use the tracer rval isolated above
		for u in range(self.Nsnap):
			smax=self.Nsnap-u
			if self.Nvariable:
				if self.tracer:
					if self.geom.periodic:
						self.msd[u]=np.sum(np.sum(np.sum((self.geom.ApplyPeriodic33(self.rtracers[:smax,:,:],self.rtracers[u:,:,:]))**2,axis=2),axis=1),axis=0)/(self.Ntracer*smax)
					else:
						self.msd[u]=np.sum(np.sum(np.sum((self.rtracers[u:,:,:]-self.rtracers[:smax,:,:])**2,axis=2),axis=1),axis=0)/(self.Ntracer*smax)
				else:
					print "Sorry: MSD for dividing particles is ambiguous and currently not implemented!"
			else:
				if self.geom.periodic:
					self.msd[u]=np.sum(np.sum(np.sum((self.geom.ApplyPeriodic33(self.rval[:smax,:,:],self.rval[u:,:,:]))**2,axis=2),axis=1),axis=0)/(self.N*smax)
				else:
					self.msd[u]=np.sum(np.sum(np.sum((self.rval[u:,:,:]-self.rval[:smax,:,:])**2,axis=2),axis=1),axis=0)/(self.N*smax)
		xval=np.linspace(0,self.Nsnap*self.param.dt*self.param.dump['freq'],num=self.Nsnap)
		if self.debug:
			fig=plt.figure()
			plt.loglog(xval,self.msd,'r.-',lw=2)
			plt.loglog(xval,self.msd[1]/(1.0*xval[1])*xval,'-',lw=2,color=[0.5,0.5,0.5])
			plt.xlabel('time')
			plt.ylabel('MSD')
			
			plt.show()
		return xval, self.msd
			
	def getVelcorr(self,dx):
		# start with the isotropic one - since there is no obvious polar region
		# and n is not the relevant variable, and v varies too much
		rangebin=0.5*np.sqrt(self.param.lx**2+self.param.ly**2)
		npts=int(round(rangebin/dx))
		print npts
		bins=np.linspace(0,rangebin,npts)
		velcorr=np.zeros((npts,))
		velav=np.zeros((self.Nsnap,3))
		for u in range(self.Nsnap):
			velcount=np.zeros((npts,))
			velcorr0=np.zeros((npts,))
			velav[u,:]=np.sum(self.vval[u,:,:],axis=0)/self.N
			print "Average velocity: ", velav[u,:]
			for k in range(self.N):
				vdot=np.sum(self.vval[u,k,:]*self.vval[u,:,:],axis=1)
				dr=self.geom.GeodesicDistance12(self.rval[u,k,:],self.rval[u,:,:])
				drbin=(np.round(dr/dx)).astype(int)
				for l in range(npts):
					pts=np.nonzero(drbin==l)[0]
					velcorr0[l]+=sum(vdot[pts])
					velcount[l]+=len(pts)
			print velcount
			isdata=[index for index, value in enumerate(velcount) if value>0]
			print isdata
			velcorr0[isdata]=velcorr0[isdata]/velcount[isdata] - np.sum(velav[u,:]*velav[u,:])
			print velcorr0
			velcorr[isdata]+=velcorr0[isdata]/velcorr0[0]
		velcorr/=self.Nsnap
		if self.debug:
			fig=plt.figure()
			isdata=[index for index, value in enumerate(velcount) if value>0]
			plt.plot(bins[isdata],velcorr[isdata],'.-')
			#plt.show()
			plt.xlabel("r-r'")
			plt.ylabel('Correlation')
		return bins,velcorr,fig
	  
	# Project our displacements or any stuff like that onto the eigenmodes of a hessian matrix, which has been calculated separately
	# we will need self.eigval and self.eigvec
	# I assume that the global skip has already taken care of any of the transient stuff
	# I am *not* removing any dreaded rattlers, because they should be part of the whole thing. 
	def projectModes(self,Hessian):
		if self.Nvariable:
			print "Hessians and dividing particles don't mix! Stopping here!"
			self.proj=0
			self.projv=0
		else:
			# self.rval and self.vval is where the fun is, self.rval=np.zeros((self.Nsnap,self.N,3))
			self.proj=np.zeros((3*Hessian.N,self.Nsnap))
			self.projv=np.zeros((3*Hessian.N,self.Nsnap))
			#proj2=np.zeros((3*Hessian.N,self.Nsnap))
			for u in range(self.Nsnap):
				dr=self.geom.ApplyPeriodic2d(self.rval[u,:,:]-Hessian.rval)
				# aah. modulo periodic boundary conditions
				dv=self.vval[u,:,:]
				# serious WTF
				#if self.debug:
					#if u==100:
						#plt.figure()
						#plt.quiver(self.rval[u,:,0],self.rval[u,:,1],dr[:,0],dr[:,1])
						#plt.title('Displacements')
						
						#plt.figure()
						#plt.quiver(self.rval[u,:,0],self.rval[u,:,1],dv[:,0],dv[:,1])
						#plt.title('Velocities')
				# now project onto the modes
				# This is the organisation of our matrix
				#plt.quiver(self.rval[:,0],self.rval[:,1],self.eigvec[0:3*self.N:3,u],self.eigvec[1:3*self.N:3,u])
				self.proj[:,u]=np.einsum('i,ij->j',dr[:,0],Hessian.eigvec[0:3*Hessian.N:3,:]) + np.einsum('i,ij->j',dr[:,1],Hessian.eigvec[1:3*Hessian.N:3,:])
				self.projv[:,u]=np.einsum('i,ij->j',dv[:,0],Hessian.eigvec[0:3*Hessian.N:3,:]) + np.einsum('i,ij->j',dv[:,1],Hessian.eigvec[1:3*Hessian.N:3,:])
				#for v in range(3*Hessian.N):
					#proj2[v,u]=np.sum(dr[:,0]*Hessian.eigvec[0:3*Hessian.N:3,u]) + np.einsum('i,ij->j',dr[:,1],Hessian.eigvec[1:3*Hessian.N:3,:])
			# projection normalization
			self.proj/=self.Nsnap
			self.projv/=self.Nsnap
			self.proj2av=np.sum(self.proj**2,axis=1)
			self.projv2av=np.sum(self.projv**2,axis=1)
		return self.proj,self.projv,self.proj2av,self.projv2av
			
	def plotProjections(self,Hessian,nmodes=5):
	  
		multmap=LinearSegmentedColormap('test',cdict,N=nmodes) 
		tval=np.linspace(0,self.Nsnap-1,self.Nsnap)
		plt.figure()
		for u in range(nmodes):
			plt.plot(tval,self.proj[u,:],'.-',color=multmap(u),label='mode '+str(u))
		plt.legend()
		plt.xlabel('time')
		plt.ylabel('projection')
		plt.title('Displacements')
		
		plt.figure()
		for u in range(nmodes):
			plt.plot(tval,self.projv[u,:],'.-',color=multmap(u),label='mode '+str(u))
		plt.legend()
		plt.xlabel('time')
		plt.ylabel('projection')
		plt.title('Velocities')
		  
		
		tau=1.0/float(self.param.nu)
		v=float(self.param.v0)
		
		plt.figure()
		plt.loglog(Hessian.eigval,0.5*Hessian.eigval*self.proj2av,'.-r',label='projection')
		# directly add the prediction here
		plt.loglog(Hessian.eigval,v**2*tau/(4.0*(1.0+Hessian.eigval*tau)),'-k',label='prediction')
		plt.xlabel(r'$\lambda$')
		plt.ylabel('Energy')
		plt.xlim(0,12)
		#plt.ylim(1e-10,1e-5)
		plt.title('Energy from displacement projections, v0=' + str(self.param.v0)+ ' tau=' + str(tau))
		
		plt.figure()
		plt.loglog(Hessian.eigval,self.proj2av,'o-r',label='displacements')
		plt.loglog(Hessian.eigval,self.projv2av,'o-g',label='velocities')
		plt.xlabel(r'$\lambda$')
		plt.ylabel('projection square')
		#plt.xlim(-0.01,12)
		plt.legend()
		plt.title('Square projections, v0=' + str(self.param.v0)+ ' tau=' + str(tau))
		
		
		plt.figure()
		plt.semilogy(Hessian.eigval,Hessian.eigval**2*self.proj2av,'o-r',label='displacements')
		plt.semilogy(Hessian.eigval,self.projv2av,'o-g',label='velocitites')
		plt.xlabel(r'$\lambda$')
		plt.ylabel('velocity projections (expected)')
		plt.xlim(-0.01,12)
		plt.legend()
		plt.title('Square velocity projections, v0=' + str(self.param.v0)+ ' tau=' + str(tau))
		
		
			
			
	


