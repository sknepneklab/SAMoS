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

from Configuration import *
#from CellList import *
from numpy import linalg as LA

try:
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
	from matplotlib.colors import LinearSegmentedColormap
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
	HAS_MATPLOTLIB=True
except:
	HAS_MATPLOTLIB=False
	pass


# For plane 2d Hessans only
class Hessian2d:
    
	def __init__(self,conf,rattlers=[],debug=False):
		self.conf=conf
		self.rval=self.conf.rval
		self.rattlers=rattlers
		self.N=self.conf.N
		self.Nrigid=self.N-len(rattlers)
		self.geom=self.conf.geom
		self.inter=self.conf.inter
		self.debug=debug
                
	def makeMatrix(self):
		# This matrix is in principle 3N by 3N. We will have to be careful later on in throwing out extra off-surface modes
		print "Hessian: Info - allocating the " + str(2*self.N) + " by " + str(2*self.N) + " 2d Hessian matrix."
		self.Hessian=np.zeros((2*self.Nrigid,2*self.N))
		# Construct it particle by particle (easiest way: everything here is always going to be N**2, no matter what, due to diagonalization algorithm)
		# Follow the formula derived in my notes
		# The unit normal for everybody will certainly be used
		Normal=self.geom.UnitNormal(self.rval)
		# So far, we really only have sphere or plane
		if (self.geom.manifold=='plane'):
			print "Hessian: Calculating Hessian on a plane!"
		else:
			print "Hessian: Error: 2d Hessian has not yet been implemented on " + self.geom.manifold + " manifolds!"  
		fsum=0.0
		fav=0.0
		for i in range(self.N):
			if i not in self.rattlers:
				if (i%200==0):
					print i
				# get some of the constants that are necessary here:
				neighbours, drvec, dr=self.conf.getNeighbours(i,self.inter.getMult(),self.inter.getDmax())
				# contact normal vectors
				nij=np.transpose(np.transpose(drvec)/np.transpose(dr))
				# Forces
				fvec=self.inter.getForce(i,neighbours,drvec,dr)
				fsum+=sum(fvec)
				# Projected onto the contact normal
				fval=np.sum(fvec*nij,axis=1)
				fav+=sum(fval)
				# Stiffnesses
				kij=self.inter.getStiffness(i,neighbours,drvec,dr)
				# equilibrium distances are given by dr already
				# Alright: elements are labeled as follows: Contact ij has sub-square 3i, 3i+1, 3i+2 and 3j, 3j+1, 3j+2
				diagsquare=np.zeros((2,2))
				for j in range(len(neighbours)):
					n=nij[j,:]
					N=Normal[neighbours[j],:]
					subsquare=np.zeros((2,2))
					# xx, xy and xz
					subsquare[0,0]=-fval[j]/dr[j]*(1-n[0]*n[0]-N[0]*N[0])+kij[j]*n[0]*n[0]
					subsquare[0,1]=-fval[j]/dr[j]*(0-n[0]*n[1]-N[0]*N[1])+kij[j]*n[0]*n[1]
					# yx, yy and yz
					subsquare[1,0]=-fval[j]/dr[j]*(0-n[1]*n[0]-N[1]*N[0])+kij[j]*n[1]*n[0]
					subsquare[1,1]=-fval[j]/dr[j]*(1-n[1]*n[1]-N[1]*N[1])+kij[j]*n[1]*n[1]
					# Stick into the big matrix
					label=neighbours[j]
					self.Hessian[2*i:(2*i+2),2*label:(2*label+2)]=-subsquare
					# Add the required bits to the diagonal part of the matrix
					# xx, xy and xz
					diagsquare[0,0]+=fval[j]/dr[j]*(1-n[0]*n[0]-N[0]*N[0])-kij[j]*n[0]*n[0]
					diagsquare[0,1]+=fval[j]/dr[j]*(0-n[0]*n[1]-N[0]*N[1])-kij[j]*n[0]*n[1]
					# yx, yy and yz
					diagsquare[1,0]+=fval[j]/dr[j]*(0-n[1]*n[0]-N[1]*N[0])-kij[j]*n[1]*n[0]
					diagsquare[1,1]+=fval[j]/dr[j]*(1-n[1]*n[1]-N[1]*N[1])-kij[j]*n[1]*n[1]
				#print diagsquare
				self.Hessian[2*i:(2*i+2),2*i:(2*i+2)]=-diagsquare
		fav/=self.N
		print "Hessian: Estimating distance from mechanical equilibrium of initial configuration "
		print "Scaled force sum is " + str(fsum/fav)
			
	def getModes(self):
		# Let's have a look if what we get is in any way reasonable
		# Eigenvalues and eigenvectors
		# Only symmetrise to calculate - for clarity and debugging above
		HessianSym=0.5*(self.Hessian+np.transpose(self.Hessian))
		if self.debug:
			plt.figure()
			plt.pcolor(HessianSym)
		# Use routines for hermitian eigenvector decomposition
		# Default is ascending order, which suits us
		print "Starting Diagonalisation!"
		self.eigval, self.eigvec = LA.eigh(HessianSym)
		if self.debug:
			# start with some debugging output
			plt.figure()
			eigrank=np.linspace(0,2*self.N,2*self.N)
			plt.plot(eigrank,self.eigval,'.-')
                wx=np.zeros((self.N,))
                wy=np.zeros((self.N,))
                for u in range(self.N):
                        # dimensional contributions
                        wx[u]=np.sum(self.eigvec[0:2*self.N:2,u]**2)
                        wy[u]=np.sum(self.eigvec[1:2*self.N:2,u]**2)
		print "The smallest eigenvalue is: " + str(np.amin(self.eigval))
		print self.eigval
		print wx
		print wy
		#print wz
	
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
		
	# project the modes into Fourier space to see how it's scaling
	def ModesFourierLongTrans(self,whichmode,qmax=0.3,verbose=True):
		eps=0.001
		print "Fourier transforming mode" + str(whichmode)
		dq=2.0*np.pi/self.geom.Lx
		nq=int(qmax/dq)
		print "Stepping Fourier transform with step " + str(dq)+ ", resulting in " + str(nq)+ " steps."
		qx, qy, qrad, ptsx, ptsy=self.makeQrad(dq,qmax,nq)
		fourierlong0=np.zeros((nq,nq),dtype=complex)
		fouriertrans0=np.zeros((nq,nq),dtype=complex)
		# #self.eigvec[0:2*self.N:2,u]
		eigx=self.eigvec[0:2*self.N:2,whichmode]
		eigy=self.eigvec[1:2*self.N:2,whichmode]
		for kx in range(nq):
			for ky in range(nq):
				# we need to be doing longitudinal and transverse here
				# Both have the same FT, but the local bits are q . e, and q X e
				fourierlong0[kx,ky]=np.sum(np.exp(1j*(qx[kx]*self.rval[:,0]+qy[ky]*self.rval[:,1]))*(eigx*qx[kx]+eigy*qy[ky])/np.sqrt(qx[kx]**2+qy[ky]**2+eps))/len(self.rval[:,0])
				fouriertrans0[kx,ky]=np.sum(np.exp(1j*(qx[kx]*self.rval[:,0]+qy[ky]*self.rval[:,1]))*(eigx*qy[ky]-eigy*qx[kx])/np.sqrt(qx[kx]**2+qy[ky]**2+eps))/len(self.rval[:,0])
		sqlong=np.real(fourierlong0**2+np.imag(fourierlong0)**2)
		sqtrans=np.real(fouriertrans0**2+np.imag(fouriertrans0)**2)
		# Produce a radial averaging to see if anything interesting happens
		nq2=int(2**0.5*nq)
		valrad=np.zeros((nq2,2))
		fourierlong=np.zeros((nq2,))
		fouriertrans=np.zeros((nq2,))
		for l in range(nq2):
			fourierlong[l]=np.mean(sqlong[ptsx[l],ptsy[l]])
			fouriertrans[l]=np.mean(sqtrans[ptsx[l],ptsy[l]])
		if verbose:
			plt.figure()
			plt.plot(qrad,fourierlong,'.-k')
			plt.plot(qrad,fouriertrans,'.-r')
			plt.xlabel('q')
			plt.ylabel('|xi_q|^2')
			plt.title('Mode ' + str(whichmode))
		return qrad,fourierlong,fouriertrans
		
	# project the modes into Fourier space to see how it's scaling
	def ModesFourier(self,whichmode,qmax=0.3,verbose=True):
		eps=0.001
		print "Fourier transforming mode" + str(whichmode)
		dq=2.0*np.pi/self.geom.Lx
		nq=int(qmax/dq)
		print "Stepping Fourier transform with step " + str(dq)+ ", resulting in " + str(nq)+ " steps."
		qx, qy, qrad, ptsx, ptsy=self.makeQrad(dq,qmax,nq)
		fouriertrans=np.zeros((nq,nq,2),dtype=complex)
		eigx=self.eigvec[0:2*self.N:2,whichmode]
		eigy=self.eigvec[1:2*self.N:2,whichmode]
		for kx in range(nq):
			for ky in range(nq):
				# we need to be doing longitudinal and transverse here
				# Both have the same FT, but the local bits are q . e, and q X e
				fouriertrans[kx,ky,0]=np.sum(np.exp(1j*(qx[kx]*self.rval[:,0]+qy[ky]*self.rval[:,1]))*eigx)/self.N #len(self.rval[:,0])
				fouriertrans[kx,ky,1]=np.sum(np.exp(1j*(qx[kx]*self.rval[:,0]+qy[ky]*self.rval[:,1]))*eigy)/self.N # len(self.rval[:,0])
                Sq=np.real(fouriertrans[:,:,0])**2+np.imag(fouriertrans[:,:,0])**2+np.real(fouriertrans[:,:,1])**2+np.imag(fouriertrans[:,:,1])**2
		# Produce a radial averaging to see if anything interesting happens
		nq2=int(2**0.5*nq)
		Sqrad=np.zeros((nq2,))
		for l in range(nq2):
                        Sqrad[l]=np.mean(Sq[ptsx[l],ptsy[l]])
		return qrad,Sqrad
	
	# get the elastic moduli cleanly once and for all. We won't diagonalise the dynamical matrix, but Fourier transform it instead
	# We have the dynamical matrix in real space as self.Hessian
	def getModuli(self,qmax=1.5,verbose=True):
                print "Fourier transforming Hessian"
		dq=2.0*np.pi/self.geom.Lx
		nq=int(qmax/dq)
		print "Stepping Fourier transform with step " + str(dq)+ ", resulting in " + str(nq)+ " steps."
                qx, qy, qrad, ptsx, ptsy=self.makeQrad(dq,qmax,nq)
		print "After qrad"
		longitudinal0=np.zeros((nq,nq))
		transverse0=np.zeros((nq,nq))
		for k in range(nq):
                        kx=qx[k]
			for l in range(nq):
                            ky=qy[l]
                            if verbose:
                                print kx
                                print ky
                            # In Fourier space, for a given k (vector), we define the 2x2 k hessian as
                            khessian=np.zeros((2,2),dtype=complex)
                            # Hessian is defined with an upfront minus sign for some reason ...
                            #khessian[0,0]=-np.sum(np.exp(1j*kx*self.rval[:,0])*np.sum(self.Hessian[0:2*self.N:2,0:2*self.N:2]*np.exp(1j*kx*self.rval[:,0]),axis=1),axis=0)
                            #khessian[0,1]=-np.sum(np.exp(1j*kx*self.rval[:,0])*np.sum(self.Hessian[0:2*self.N:2,1:2*self.N:2]*np.exp(1j*ky*self.rval[:,1]),axis=1),axis=0)
                            #khessian[1,0]=-np.sum(np.exp(1j*ky*self.rval[:,1])*np.sum(self.Hessian[1:2*self.N:2,0:2*self.N:2]*np.exp(1j*kx*self.rval[:,0]),axis=1),axis=0)
                            #khessian[1,1]=-np.sum(np.exp(1j*ky*self.rval[:,1])*np.sum(self.Hessian[1:2*self.N:2,1:2*self.N:2]*np.exp(1j*ky*self.rval[:,1]),axis=1),axis=0)
                            khessian[0,0]=np.dot(np.exp(1j*kx*self.rval[:,0]),np.dot(self.Hessian[0:2*self.N:2,0:2*self.N:2],np.exp(-1j*kx*self.rval[:,0])))/self.N
                            khessian[0,1]=np.dot(np.exp(1j*kx*self.rval[:,0]),np.dot(self.Hessian[0:2*self.N:2,1:2*self.N:2],np.exp(-1j*ky*self.rval[:,1])))/self.N
                            khessian[1,0]=np.dot(np.exp(1j*ky*self.rval[:,1]),np.dot(self.Hessian[1:2*self.N:2,0:2*self.N:2],np.exp(-1j*kx*self.rval[:,0])))/self.N
                            khessian[1,1]=np.dot(np.exp(1j*ky*self.rval[:,1]),np.dot(self.Hessian[1:2*self.N:2,1:2*self.N:2],np.exp(-1j*ky*self.rval[:,1])))/self.N
                            # Its eigenvalues are then (B+\mu) k^2 and k^2, in principle
                            eigk, eigveck = LA.eig(khessian)
                            # We won't get proper, pure longitudinal and transverse eigenvectors
                            # project them onto k, and take the one more along k as the longitudinal one
                            proj1= kx*eigveck[0,0]+ky*eigveck[0,1]
                            proj2= kx*eigveck[1,0]+ky*eigveck[1,1]
                            if abs(proj2)>abs(proj1):
                                longitudinal0[k,l]=eigk[1]
                                transverse0[k,l]=eigk[0]
                            else:
                                longitudinal0[k,l]=eigk[0]
                                transverse0[k,l]=eigk[1]
                            if verbose:
                                print "Found eigenvalues long " + str(longitudinal0[k,l]) + " and transverse " + str(transverse0[k,l])
                nq2=int(2**0.5*nq)
		longitudinal=np.zeros((nq2,))
		transverse=np.zeros((nq2,))
		for l in range(nq2):
			longitudinal[l]=np.mean(longitudinal0[ptsx[l],ptsy[l]])
			transverse[l]=np.mean(transverse0[ptsx[l],ptsy[l]])
		
		if verbose:
			plt.figure()
			plt.plot(qrad**2,longitudinal,'ok')
			plt.plot(qrad**2,transverse,'or')
                return qrad, longitudinal, transverse
                
                
                
		
		