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
	def __init__(self,directory,conffile,inputfile,boundaryfile,skip,tracer=False,ignore=False,takeDrift=False):
		self.tracer=tracer
		self.ignore=ignore
		self.takeDrift=takeDrift
		self.param = Param(directory+conffile)
		files = sorted(glob(directory + self.param.dumpname+'*.dat'))[skip:]
		if len(files) == 0:
  			files = sorted(glob(directory + self.param.dumpname+'*.dat.gz'))[skip:]
		# Read the local data
		# Outdated list of geometries
		#geometries={'sphere':GeometrySphere,'plane':GeometryPeriodicPlane,'none':Geometry,'tube':GeometryTube,'peanut':GeometryPeanut,'hourglass':GeometryHourglass}
		geometries={'sphere':GeometrySphere,'plane':GeometryPlane,'plane_periodic':GeometryPeriodicPlane,'none':Geometry,'tube':GeometryTube,'peanut':GeometryPeanut,'hourglass':GeometryHourglass}
		# Create the right geometry environment (TBC):
		self.geom=geometries[self.param.constraint](self.param)
		print self.geom
		# Number of files that we are dealing with
		self.Nsnap=len(files)
		# Find out if this is a simulation with a variable number of particles:
		# For now, otherwise we get a mess ...
		self.Nvariable=True
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
		self.flag=np.zeros((self.Nsnap,self.N))
		if self.Nvariable:
			#self.flag=np.zeros((self.Nsnap,self.N))
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
			if self.Nvariable:
				self.rval[u,:self.Nval[u],:]=rval_u
				self.vval[u,:self.Nval[u],:]=vval_u
				if tracer:
					self.rtracers[u,:,:]=self.rval[u,tracers,:]
					self.vtracers[u,:,:]=self.rval[u,tracers,:]
			else:
				self.rval[u,:,:]=rval_u
				self.vval[u,:,:]=vval_u
			u+=1
		# Take the drift off as a post-processing step if desired
		# Rather: Just compute it, and take it off in the MSD etc
		if self.takeDrift:
			self.drift=np.zeros((self.Nsnap,3))
			if self.Nvariable:
				print "Variable N: Taking off the drift is meaningless. Doing nothing."
			else:	
				for u in range(1,self.Nsnap):
					 dr=self.geom.ApplyPeriodic2d(self.rval[u,:,:]-self.rval[u-1,:,:])
					 drift0=np.sum(dr,axis=0)/self.N
					 self.drift[u,:]=self.drift[u-1,:]+drift0
					 print self.drift[u,:]
			  
		
	def getMSD(self,verbose=True):
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
				# Drift is only meaningful here
				#print "Doing the non-tracer MSD: this may take some time"
				if self.takeDrift:
					hmm=(self.drift[:smax,:]-self.drift[u:,:])
					takeoff=np.einsum('j,ik->ijk',np.ones((self.N,)),hmm)
					if self.geom.periodic:
						dr=self.geom.ApplyPeriodic3d(self.rval[:smax,:,:]-self.rval[u:,:,:])-takeoff
					else:
						dr=self.rval[:smax,:,:]-self.rval[u:,:,:]-takeoff
					self.msd[u]=np.sum(np.sum(np.sum(dr**2,axis=2),axis=1),axis=0)/(self.N*smax)
				else:
					if self.geom.periodic:
						self.msd[u]=np.sum(np.sum(np.sum((self.geom.ApplyPeriodic33(self.rval[:smax,:,:],self.rval[u:,:,:]))**2,axis=2),axis=1),axis=0)/(self.N*smax)
					else:
						self.msd[u]=np.sum(np.sum(np.sum((self.rval[u:,:,:]-self.rval[:smax,:,:])**2,axis=2),axis=1),axis=0)/(self.N*smax)
		xval=np.linspace(0,self.Nsnap*self.param.dt*self.param.dump['freq'],num=self.Nsnap)
		if verbose:
			fig=plt.figure()
			plt.loglog(xval,self.msd,'r.-',lw=2)
			plt.loglog(xval,self.msd[1]/(1.0*xval[1])*xval,'-',lw=2,color=[0.5,0.5,0.5])
			plt.xlabel('time')
			plt.ylabel('MSD')
			
			#plt.show()
		return xval, self.msd

	
	# Definition of the self-intermediate scattering function (Flenner + Szamel)
	# 1/N <\sum_n exp(iq[r_n(t)-r_n(0)]>_t,n
	def SelfIntermediate(self,qval,verbose=True):
		# This is single particle, single q, shifted time step. Equivalent to the MSD, really
		SelfInt=np.empty((self.Nsnap,),dtype=complex)
		for u in range(self.Nsnap):
			smax=self.Nsnap-u
			if self.Nvariable:
				if self.tracer:
					if self.geom.periodic:
						SelfInt[u]=np.sum(np.sum(np.exp(1.0j*qval[0]*(self.geom.ApplyPeriodicX(-self.rtracers[:smax,:,0]+self.rtracers[u:,:,0]))+1.0j*qval[1]*(self.geom.ApplyPeriodicY(-self.rtracers[:smax,:,1]+self.rtracers[u:,:,1]))+1.0j*qval[2]*(self.geom.ApplyPeriodicZ(-self.rtracers[:smax,:,2]+self.rtracers[u:,:,2]))),axis=1),axis=0)/(self.Ntracer*smax)
					else:
						SelfInt[u]=np.sum(np.sum(np.exp(1.0j*qval[0]*(-self.rtracers[:smax,:,0]+self.rtracers[u:,:,0])+1.0j*qval[1]*(-self.rtracers[:smax,:,1]+self.rtracers[u:,:,1])+1.0j*qval[2]*(-self.rtracers[:smax,:,2]+self.rtracers[u:,:,2])),axis=1),axis=0)/(self.Ntracer*smax)
				else:
					print "Sorry: Self-intermediate scattering function for dividing particles is ambiguous and currently not implemented!"
			else:
				if self.geom.periodic:
					SelfInt[u]=np.sum(np.sum(np.exp(1.0j*qval[0]*(self.geom.ApplyPeriodicX(-self.rval[:smax,:,0]+self.rval[u:,:,0]))+1.0j*qval[1]*(self.geom.ApplyPeriodicY(-self.rval[:smax,:,1]+self.rval[u:,:,1]))+1.0j*qval[2]*(self.geom.ApplyPeriodicZ(-self.rval[:smax,:,2]+self.rval[u:,:,2]))),axis=1),axis=0)/(self.N*smax)
				else:
					SelfInt[u]=np.sum(np.sum(np.exp(1.0j*qval[0]*(-self.rval[:smax,:,0]+self.rval[u:,:,0])+1.0j*qval[1]*(-self.rval[:smax,:,1]+self.rval[u:,:,1])+1.0j*qval[2]*(-self.rval[:smax,:,2]+self.rval[u:,:,2])),axis=1),axis=0)/(self.N*smax)
		tval=np.linspace(0,self.Nsnap*self.param.dt*self.param.dump['freq'],num=self.Nsnap)
		# Looking at the absolute value of it here
		SelfInt2=(np.real(SelfInt)**2 + np.imag(SelfInt)**2)**0.5
		qnorm=np.sqrt(qval[0]**2+qval[1]**2+qval[2]**2)
		if verbose:
			fig=plt.figure()
			plt.semilogx(tval,SelfInt2,'.-r',lw=2)
			plt.xlabel('time')
			plt.ylabel('F_s(k,t)')
			plt.title('Self-intermediate, k = ' + str(qnorm))
			#plt.show()
		return tval, SelfInt2
		
					
	def FourierTrans(self,whichframe,qmax=0.3,verbose=True):
		# Note to self: only low q values will be interesting in any case. 
		# The stepping is in multiples of the inverse box size. Assuming a square box.
		print "Fourier transforming positions"
		dq=1.0/self.geom.Lx
		nq=int(qmax/dq)
		print "Stepping Fourier transform with step " + str(dq)+ ", resulting in " + str(nq)+ " steps."
		qx, qy, qrad, ptsx, ptsy=self.makeQrad(dq,qmax,nq)
		fourierval=np.zeros((nq,nq),dtype=complex)
		for kx in range(nq):
			for ky in range(nq):
				# And, alas, no FFT since we are most definitely off grid. And averaging is going to kill everything.
				fourierval[kx,ky]=np.sum(np.exp(1j*(qx[kx]*self.rval[whichframe,:,0]+qy[ky]*self.rval[whichframe,:,1])))
		plotval=np.real(fourierval)**2+np.imag(fourierval)**2
		# Produce a radial averaging to see if anything interesting happens
		nq2=int(2**0.5*nq)
		valrad=np.zeros((nq2,))
		for l in range(nq2):
			valrad[l]=np.mean(plotval[ptsx[l],ptsy[l]])
		
		if verbose:
			plt.figure()
			plt.pcolor(qx,qy,plotval)
			plt.colorbar()
			plt.title('Positions')
		return qrad,valrad
	  
	def makeQrad(self,dq,qmax,nq):
		nq2=int(2**0.5*nq)
		qmax2=2**0.5*qmax
		qx=np.linspace(0,qmax,nq)
		qy=np.linspace(0,qmax,nq)
		qrad=np.linspace(0,qmax2,nq2)
		# do this silly counting once and for all
		binval=np.empty((nq,nq))
		for kx in range(nq):
			for ky in range(nq):
				qval=np.sqrt(qx[kx]**2+qy[ky]**2)
				binval[kx,ky]=round(qval/dq)
		#print binval
		ptsx=[]
		ptsy=[]
		# do the indexing arrays
		for l in range(nq2):
			pts0x=[]
			pts0y=[]
			for kx in range(nq):
				hmm=np.nonzero(binval[kx,:]==l)[0]
				for v in range(len(hmm)):
					pts0y.append(hmm[v])
					pts0x.append(kx)
			ptsx.append(pts0x)
			ptsy.append(pts0y)
		return qx, qy, qrad, ptsx, ptsy
	  
	# Well, that seems to do fuck all
	def getDynStruct(self,qmax,omegamax,verbose=True,nmax=50):
		# Following the template of Wysocki, Winkler, Gompper
		dq=1.0/self.geom.Lx
		nq=int(qmax/dq)
		if nq>nmax:
			print "Coarsening q interval to reduce computational load"
			nq=nmax
			dq=qmax/nq
		nq2=int(2**0.5*nq)
		print "Stepping space Fourier transform with step " + str(dq)+ ", resulting in " + str(nq)+ " steps."
		dom=1.0/(self.Nsnap*self.param.dt*self.param.dump['freq'])
		nom1=int(omegamax/dom)
		nom=2*int(omegamax/dom)+1
		print "Stepping time Fourier transform with step " + str(dom)+ ", resulting in " + str(nom)+ " steps."
		# Formally: S(q,omega) = 1/N int dt \rho_q(t) \rho*_q(0) e^i\omega t, where \rho_q(t) = \int dr \rho(r,t) e^iq r
		# The second part is what we already had for the positional static structure factor
		# For simplicity reasons, do the radial averaging before taking the time transform
		qx, qy, qrad, ptsx, ptsy=self.makeQrad(dq,qmax,nq)
		rhorad=np.zeros((self.Nsnap,nq2),dtype=complex)
		for u in range(self.Nsnap):
			if (u%10==0):
				print u
			fourierval=np.empty((nq,nq),dtype=complex)
			for kx in range(nq):
				for ky in range(nq):
					# And, alas, no FFT since we are most definitely off grid. And averaging is going to kill everything.
					fourierval[kx,ky]=np.sum(np.exp(1j*(qx[kx]*self.rval[u,:,0]+qy[ky]*self.rval[u,:,1])))
			for l in range(nq2):
				rhorad[u,l]=np.mean(fourierval[ptsx[l],ptsy[l]])
		# Do our little shifted averaging procedure at constant q now
		rhocorr=np.zeros((self.Nsnap,nq2),dtype=complex)
		for u in range(self.Nsnap):
			smax=self.Nsnap-u
			rhocorr[u,:]=np.sum(rhorad[u:,:]*rhorad[:smax,:],axis=0)/smax
		# Cute. Now do the time tranform:
		DynStruct=np.zeros((nom,nq2),dtype=complex)
		tval=np.linspace(0,self.Nsnap*self.param.dt*self.param.dump['freq'],num=self.Nsnap)
		omega=np.empty((nom,))
		for no in range(0,nom):
			omega[no]=(-nom1+no)*dom
			DynStruct[no,:]=np.einsum('ij,i', rhocorr, np.exp(1j*omega[no]*tval))
		print omega
		# OK, what have we got? Take the absolute value and look
		PlotDynStruct=np.real(DynStruct)**2+np.imag(DynStruct)**2
		if verbose:
			plt.figure()
			plt.pcolor(qrad,omega,np.log10(PlotDynStruct))
			plt.colorbar()
			plt.title('Dynamical structure factor')
			
			plt.figure()
			plt.pcolor(qrad,tval,np.log10(np.real(rhocorr)))
			plt.colorbar()
			plt.title('Density correlation function')
		if verbose:
			plt.figure()
			plt.plot(omega,np.log10(PlotDynStruct[:,0]),'.-k')
			plt.plot(omega,np.log10(PlotDynStruct[:,1]),'.-r')
			plt.plot(omega,np.log10(PlotDynStruct[:,5]),'.-g')
			plt.plot(omega,np.log10(PlotDynStruct[:,10]),'.-b')
			plt.xlabel('omega')
			plt.ylabel('structure factor')
			plt.title('Dynamical structure factor')

			
			plt.figure()
			plt.plot(tval,np.log10(rhocorr[:,0]),'.-k')
			plt.plot(tval,np.log10(rhocorr[:,1]),'.-g')
			plt.plot(tval,np.log10(rhocorr[:,5]),'.-r')
			plt.plot(tval,np.log10(rhocorr[:,10]),'.-b')
			plt.title('Density correlation function')
		return omega,qrad,DynStruct
		
	# Four point structure factor
	def FourPoint(self,a,qmax=3.14,verbose=True,nmax=20):
		# As written, this thing only works with tracer particles, since I need to track them through the whole simulation
		# Following the template of Wysocki, Winkler, Gompper
		dq=1.0/self.geom.Lx
		nq=int(qmax/dq)
		if nq>nmax:
			print "Coarsening q interval to reduce computational load"
			nq=nmax
			dq=qmax/nq
		nq2=int(2**0.5*nq)
		print "Stepping space Fourier transform with step " + str(dq)+ ", resulting in " + str(nq)+ " steps."
		qx, qy, qrad, ptsx, ptsy=self.makeQrad(dq,qmax,nq)
		FourPoint=np.zeros((nq2,self.Nsnap))
		for u in range(self.Nsnap):
			if (u%10==0):
				print u
			smax=self.Nsnap-u
			if self.Nvariable:
				if self.tracer:
					# First filter out the particles we are dealing with: only those that have moved less than distance a
					if self.geom.periodic:
						dists=np.sqrt(np.sum(self.geom.ApplyPeriodic3d(self.rtracers[u:,:,:]-self.rtracers[:smax,:,:])**2,axis=2))
					else:
						dists=np.sqrt(np.sum((self.rtracers[u:,:,:]-self.rtracers[:smax,:,:])**2,axis=2))
					# So now replace those by 1s or 0s depending on how far they have gone
					# this is 0 or 1 as we want, except for negative ones where it's gone far
					hmm=1.0-np.round(dists/a)
					# remove the negative ones
					weights=0.5*(hmm+abs(hmm))
					# Properly: S_4 = 1/N < \sum_n \sum_m w_n(t) w_m(t) e^{iq.(r_n(0)-r_m(0)}>_m,n,t
					# For us this triple sum over m, n and time shift needs to be done here. 
					# Then do the radial averaging in the last step.
					# we can decompose the fourier transform into
					# S_4 = 1/N < \sum_n w_n(t) e^{iq.r_n(0)} \sum_m w_m(t) e^{-iq r_m(0)}>_t
					# So do these first, then radially average, finally take the t average 
					fourierval=np.zeros((nq,nq,smax),dtype=complex)
					fourierrad=np.zeros((nq2,smax),dtype=complex)
					for kx in range(nq):
						for ky in range(nq):
							fourierval[kx,ky,:]=np.sum(weights*np.exp(1j*qx[kx]*self.rtracers[u:,:,0]+1j*qy[ky]*self.rtracers[u:,:,1]),axis=1)
					for l in range(nq2):
						fourierrad[l,:]=np.mean(fourierval[ptsx[l],ptsy[l],:])
					# So now finally multiply and do the shifted average over time. PBC should have been sorted out right in dists? Or not?
					# Should be real at that point
					FourPoint[:,u]=np.real(np.sum(fourierrad*np.conjugate(fourierrad),axis=1))/(self.Ntracer*smax)
				else:
					print "Sorry: Four point function for dividing particles is ambiguous and currently not implemented!"
			else:
				# First filter out the particles we are dealing with: only those that have moved less than distance a
				#print "before distances"
				if self.geom.periodic:
					dists=np.sqrt(np.sum(self.geom.ApplyPeriodic3d(self.rval[u:,:,:]-self.rval[:smax,:,:])**2,axis=2))
				else:
					dists=np.sqrt(np.sum((self.rval[u:,:,:]-self.rval[:smax,:,:])**2,axis=2))
				# So now replace those by 1s or 0s depending on how far they have gone
				# this is 0 or 1 as we want, except for negative ones where it's gone far
				#print "before weights"
				hmm=1.0-np.round(dists/a)
				# remove the negative ones
				weights=0.5*(hmm+abs(hmm))
				# Properly: S_4 = 1/N < \sum_n \sum_m w_n(t) w_m(t) e^{iq.(r_n(0)-r_m(0)}>_m,n,t
				# For us this triple sum over m, n and time shift needs to be done here. 
				# Then do the radial averaging in the last step.
				# we can decompose the fourier transform into
				# S_4 = 1/N < \sum_n w_n(t) e^{iq.r_n(0)} \sum_m w_m(t) e^{-iq r_m(0)}>_t
				# So do these first, then radially average, finally take the t average 
				#print "before Fourier"
				fourierval=np.zeros((nq,nq,smax),dtype=complex)
				fourierrad=np.zeros((nq2,smax),dtype=complex)
				for kx in range(nq):
					for ky in range(nq):
						fourierval[kx,ky,:]=np.sum(weights[:,:]*np.exp(1j*(qx[kx]*self.rval[u:,:,0])+qy[ky]*self.rval[u:,:,1]),axis=1)
				#print "before radial average"
				for l in range(nq2):
					fourierrad[l,:]=np.mean(fourierval[ptsx[l],ptsy[l],:])
				# So now finally multiply and do the shifted average over time. PBC should have been sorted out right in dists? Or not?
				# Should be real at that point
				#print "before fourpoint"
				FourPoint[:,u]=np.real(np.sum(fourierrad*np.conjugate(fourierrad),axis=1))/(self.N*smax)
		tval=np.linspace(0,self.Nsnap*self.param.dt*self.param.dump['freq'],num=self.Nsnap)
		if verbose:
			plt.figure()
			vmap=LinearSegmentedColormap('test',cdict,N=nq2) 
			for q in range(0,nq2):
				plt.loglog(tval,FourPoint[q,:],'.-',color=vmap(q))
			plt.xlabel('t')
			plt.ylabel('FourPoint')			
		return tval, FourPoint
	
	def FourierTransVel(self,whichframe,qmax=0.3,verbose=True):
		# Note to self: only low q values will be interesting in any case. 
		# The stepping is in multiples of the inverse box size. Assuming a square box.
		print "Fourier transforming velocities"
		dq=1.0/self.geom.Lx
		nq=int(qmax/dq)
		print "Stepping Fourier transform with step " + str(dq)+ ", resulting in " + str(nq)+ " steps."
		qx, qy, qrad, ptsx, ptsy=self.makeQrad(dq,qmax,nq)
		fourierval=np.zeros((nq,nq,2),dtype=complex)
		for kx in range(nq):
			for ky in range(nq):
				# And, alas, no FFT since we are most definitely off grid. And averaging is going to kill everything.
				fourierval[kx,ky,0]=np.sum(np.exp(1j*(qx[kx]*self.rval[whichframe,:,0]+qy[ky]*self.rval[whichframe,:,1]))*self.vval[whichframe,:,0])
				fourierval[kx,ky,1]=np.sum(np.exp(1j*(qx[kx]*self.rval[whichframe,:,0]+qy[ky]*self.rval[whichframe,:,1]))*self.vval[whichframe,:,1])
		# Sq = \vec{v_q}.\vec{v_-q}, assuming real and symmetric
		# = \vec{v_q}.\vec{v_q*} = v
		Sq=np.real(fourierval[:,:,0])**2+np.imag(fourierval[:,:,0])**2+np.real(fourierval[:,:,1])**2+np.imag(fourierval[:,:,1])**2
		plotval_x=np.sqrt(np.real(fourierval[:,:,0])**2+np.imag(fourierval[:,:,0])**2)
		plotval_y=np.sqrt(np.real(fourierval[:,:,1])**2+np.imag(fourierval[:,:,1])**2)
		# Produce a radial averaging to see if anything interesting happens
		nq2=int(2**0.5*nq)
		valrad=np.zeros((nq2,2))
		Sqrad=np.zeros((nq2,))
		for l in range(nq2):
			valrad[l,0]=np.mean(plotval_x[ptsx[l],ptsy[l]])
			valrad[l,1]=np.mean(plotval_y[ptsx[l],ptsy[l]])
			Sqrad[l]=np.mean(Sq[ptsx[l],ptsy[l]])
		
		if verbose:
			plt.figure()
			plt.pcolor(qx,qy,plotval_x)
			plt.colorbar()
			plt.title('Velocities - x')
			plt.figure()
			plt.pcolor(qx,qy,plotval_y)
			plt.colorbar()
			plt.title('Velocities - y')
		return qrad,valrad,Sqrad
	  
	def getVelcorrSingle(self,whichframe,dx,xmax,verbose=True):
		# start with the isotropic one - since there is no obvious polar region
		# and n is not the relevant variable, and v varies too much
		print "Velocity correlation function for frame " + str(whichframe)
		npts=int(round(xmax/dx))
		bins=np.linspace(0,xmax,npts)
		velcorr=np.zeros((npts,))
		velcount=np.zeros((npts,))
		velav=np.sum(self.vval[whichframe,:,:],axis=0)/self.N
		for k in range(self.N):
			vdot=np.sum(self.vval[whichframe,k,:]*self.vval[whichframe,:,:],axis=1)
			dr=self.geom.GeodesicDistance12(self.rval[whichframe,k,:],self.rval[whichframe,:,:])
			drbin=(np.round(dr/dx)).astype(int)
			for l in range(npts):
				pts=np.nonzero(drbin==l)[0]
				velcorr[l]+=sum(vdot[pts])
				velcount[l]+=len(pts)
		isdata=[index for index, value in enumerate(velcount) if value>0]
		velcorr[isdata]=velcorr[isdata]/velcount[isdata] - np.sum(velav*velav)
		if verbose:
			fig=plt.figure()
			isdata=[index for index, value in enumerate(velcount) if value>0]
			plt.plot(bins[isdata],velcorr[isdata],'.-')
			#plt.show()
			plt.xlabel("r-r'")
			plt.ylabel('Correlation')
		return bins,velcorr
	  
	# This one is highly impractical
	def getVelcorr(self,dx,verbose=True):
		# start with the isotropic one - since there is no obvious polar region
		# and n is not the relevant variable, and v varies too much
		rangebin=0.5*np.sqrt(self.param.lx**2+self.param.ly**2)
		npts=int(round(rangebin/dx))
		print npts
		bins=np.linspace(0,rangebin,npts)
		velcorr=np.zeros((npts,))
		velav=np.zeros((self.Nsnap,3))
		for u in range(self.Nsnap):
			print "snapshot " + str(u)
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
			isdata=[index for index, value in enumerate(velcount) if value>0]
			velcorr0[isdata]=velcorr0[isdata]/velcount[isdata] - np.sum(velav[u,:]*velav[u,:])
			velcorr[isdata]+=velcorr0[isdata]/velcorr0[0]
		velcorr/=self.Nsnap
		if verbose:
			fig=plt.figure()
			isdata=[index for index, value in enumerate(velcount) if value>0]
			plt.plot(bins[isdata],velcorr[isdata],'.-')
			#plt.show()
			plt.xlabel("r-r'")
			plt.ylabel('Correlation')
		return bins,velcorr
	  
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
		
		
			
			
	


