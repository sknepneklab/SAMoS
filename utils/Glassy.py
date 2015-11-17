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
from read_param import *
from read_data import *
try:
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
	HAS_MATPLOTLIB=True
except:
	HAS_MATPLOTLIB=False
	pass

class SimRun:
	def __init__(self,takeDrift,directory,conffile,inputfile,skip,debug=False):
		self.debug=debug
		self.param = Param(directory+conffile)
		files = sorted(glob(directory + self.param.dumpname+'*.dat'))[skip:]
		if len(files) == 0:
  			files = sorted(glob(directory + self.param.dumpname+'*.dat.gz'))[skip:]
		# Read the local data
		geometries={'sphere':GeometrySphere,'plane':GeometryPeriodicPlane,'none':Geometry,'tube':GeometryTube,'peanut':GeometryPeanut,'hourglass':GeometryHourglass}
		# Create the right geometry environment (TBC):
		self.geom=geometries[self.param.constraint](self.param)
		print self.geom
		
		
		print "Reading radii!"
		data_ini=ReadData(directory+inputfile)
		self.radius=data_ini.data[data_ini.keys['radius']]	
		self.monodisperse=False
		self.N=len(self.radius)
		self.Nsnap=len(files)
		self.rval=np.empty((len(files),self.N,3))
		self.vval=np.empty((len(files),self.N,3))
		self.nval=np.empty((len(files),self.N,3))
		self.vhat=np.empty((len(files),self.N,3))
		u=0
		for f in files:
			print "Processing file : ", f
			data = ReadData(f)
			x, y, z = np.array(data.data[data.keys['x']]), np.array(data.data[data.keys['y']]), np.array(data.data[data.keys['z']])
			vx, vy, vz = np.array(data.data[data.keys['vx']]), np.array(data.data[data.keys['vy']]), np.array(data.data[data.keys['vz']])
			nx, ny, nz = np.array(data.data[data.keys['nx']]), np.array(data.data[data.keys['ny']]), np.array(data.data[data.keys['nz']])	
			if data.keys.has_key('flag'):
				self.flag = data.data[data.keys['flag']]
			if u>0 and takeDrift:
				rval_u0=rval_u
				print u
			rval_u = np.column_stack((x,y,z))
			vval_u = np.column_stack((vx,vy,vz))
			nval_u = np.column_stack((nx,ny,nz))
			vel = np.sqrt(vval_u[:,0]**2 + vval_u[:,1]**2 + vval_u[:,2]**2)
			vhat_u=((vval_u).transpose()/(vel).transpose()).transpose()
			# taking off the drift
			if takeDrift:
				if u>0:
					# Bring back people who crossed the line
					diff_real =self.geom.ApplyPeriodic12(rval_u0,rval_u)
					#print np.max(diff_real)
					#print diff_real
					rval_ud=rval_u0+diff_real
					drift=np.sum(diff_real,axis=0)/self.N
					#print drift
					rval_ud-=drift
					self.rval[u,:,:]=rval_ud
				else:
					self.rval[u,:,:]=rval_u
			else:
				self.rval[u,:,:]=rval_u
			self.vval[u,:,:]=vval_u
			self.nval[u,:,:]=nval_u
			self.vhat[u,:,:]=vhat_u
			u+=1
		
	def getMSD(self):
		self.msd=np.empty((self.Nsnap,))
		for u in range(self.Nsnap):
			smax=self.Nsnap-u
			#self.geom.periodic=True
			if self.geom.periodic:
				self.msd[u]=np.sum(np.sum(np.sum((self.geom.ApplyPeriodic33(self.rval[:smax,:,:],self.rval[u:,:,:]))**2,axis=2),axis=1),axis=0)/(self.N*smax)
			else:
				self.msd[u]=np.sum(np.sum(np.sum((self.rval[u:,:,:]-self.rval[:smax,:,:])**2,axis=2),axis=1),axis=0)/(self.N*smax)
		if self.debug:
			fig=plt.figure()
			xval=np.linspace(0,self.Nsnap*self.param.dt*self.param.dump['freq'],num=self.Nsnap)
			plt.loglog(xval,self.msd,'.-')
			plt.show()
			
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
			print velav[u,:]
			for k in range(self.N):
				vdot=np.sum(self.vval[u,k,:]*self.vval[u,:,:],axis=1)
				dr=self.geom.GeodesicDistance12(self.rval[u,k,:],self.rval[u,:,:])
				binplace = np.digitize(dr, bins)
				for l in range(npts):
					velcorr0[binplace[l]-1]+=sum(vdot[np.where(binplace == l)])
					velcount[binplace[l]-1]+=len(np.where(binplace == l))
			print velcount
			velcorr0=velcorr0/velcount - np.sum(velav[u,:]*velav[u,:])
			print velcorr0
			velcorr+=velcorr0
		velcorr/=self.Nsnap
		if self.debug:
			fig=plt.figure()
			plt.plot(bins,velcorr/velcorr[0],'.-')
			plt.show()
		
	# Mostly for dividing systems: simple statistics, including density, pressure, velocity
	def getStatistics(self):
		 rhoval=np.zeros((self.Nsnap,))
		 pval=np.zeros((self.Nsnap,))
		 v2val=np.zeros((self.Nsnap,))
		 for u in range(self.Nsnap):
			rhoval[u]=np.sum(
	
	#def getSelfInt(self):
				
		
