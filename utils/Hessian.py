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

#try:
	#import matplotlib.pyplot as plt
	#from mpl_toolkits.mplot3d import Axes3D
	#from matplotlib.colors import LinearSegmentedColormap
	#matplotlib.rcParams['text.usetex'] = 'false'
	#matplotlib.rcParams['lines.linewidth'] = 2
	#matplotlib.rcParams['axes.linewidth'] = 2
	#matplotlib.rcParams['xtick.major.size'] = 8
	#matplotlib.rcParams['ytick.major.size'] = 8
	#matplotlib.rcParams['font.size']=16.0
	#matplotlib.rcParams['legend.fontsize']=14.0

	#cdict = {'red':   [(0.0,  0.0, 0.5),
					  #(0.35,  1.0, 0.75),
					  #(0.45,  0.75, 0.0),
					  #(1.0,  0.0, 0.0)],

			#'green': [(0.0,  0.0, 0.0),
					  #(0.35,  0.0, 0.5),
					  #(0.5, 1.0, 1.0),
					  #(0.8,  0.5, 0.0),
					  #(1.0,  0.0, 0.0)],

			#'blue':  [(0.0,  0.0, 0.0),
					  #(0.5,  0.0, 0.0),
					  #(0.7, 0.5, 1.0),
					  #(1.0,  0.25, 0.0)]}
	#HAS_MATPLOTLIB=True
#except:
	#HAS_MATPLOTLIB=False
	#pass



class Hessian:
    
	def __init__(self,conf,rattlers=[],debug=False):
		self.conf=conf
		self.rval=self.conf.rval
		self.rattlers=rattlers
		self.N=self.conf.N
		self.Nrigid=self.N-len(rattlers)
		self.geom=self.conf.geom
		self.inter=self.conf.inter
		self.debug=debug
                
	def makeMatrix(self,addRestoring=True,ksurf=10.0):
		# This matrix is in principle 3N by 3N. We will have to be careful later on in throwing out extra off-surface modes
		print "Hessian: Info - allocating the " + str(3*self.N) + " by " + str(3*self.N) + " Hessian matrix."
		self.Hessian=np.zeros((3*self.Nrigid,3*self.N))
		# Construct it particle by particle (easiest way: everything here is always going to be N**2, no matter what, due to diagonalization algorithm)
		# Follow the formula derived in my notes
		# The unit normal for everybody will certainly be used
		Normal=self.geom.UnitNormal(self.rval)
		# So far, we really only have sphere or plane
		addCurvature=False
		if (self.geom.manifold=='sphere'):
			addCurvature=True
			Rval=self.geom.R
			print "Hessian: Calculating Hessian on a sphere!"
		elif (self.geom.manifold=='plane'):
			print "Hessian: Calculating Hessian on a plane!"
		else:
			print "Hessian: Error: Hessian has not yet been implemented on " + self.geom.manifold + " manifolds!"  
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
				diagsquare=np.zeros((3,3))
				for j in range(len(neighbours)):
					n=nij[j,:]
					N=Normal[neighbours[j],:]
					subsquare=np.zeros((3,3))
					# xx, xy and xz
					subsquare[0,0]=-fval[j]/dr[j]*(1-n[0]*n[0]-N[0]*N[0])+kij[j]*n[0]*n[0]
					subsquare[0,1]=-fval[j]/dr[j]*(0-n[0]*n[1]-N[0]*N[1])+kij[j]*n[0]*n[1]
					subsquare[0,2]=-fval[j]/dr[j]*(0-n[0]*n[2]-N[0]*N[2])+kij[j]*n[0]*n[2]
					# yx, yy and yz
					subsquare[1,0]=-fval[j]/dr[j]*(0-n[1]*n[0]-N[1]*N[0])+kij[j]*n[1]*n[0]
					subsquare[1,1]=-fval[j]/dr[j]*(1-n[1]*n[1]-N[1]*N[1])+kij[j]*n[1]*n[1]
					subsquare[1,2]=-fval[j]/dr[j]*(0-n[1]*n[2]-N[1]*N[2])+kij[j]*n[1]*n[2]
					# zx, zy and zz
					subsquare[2,0]=-fval[j]/dr[j]*(0-n[2]*n[0]-N[2]*N[0])+kij[j]*n[2]*n[0]
					subsquare[2,1]=-fval[j]/dr[j]*(0-n[2]*n[1]-N[2]*N[1])+kij[j]*n[2]*n[1]
					subsquare[2,2]=-fval[j]/dr[j]*(1-n[2]*n[2]-N[2]*N[2])+kij[j]*n[2]*n[2]
					# Stick into the big matrix
					label=neighbours[j]
					self.Hessian[3*i:(3*i+3),3*label:(3*label+3)]=-subsquare
					# Add the required bits to the diagonal part of the matrix
					# xx, xy and xz
					diagsquare[0,0]+=fval[j]/dr[j]*(1-n[0]*n[0]-N[0]*N[0])-kij[j]*n[0]*n[0]
					diagsquare[0,1]+=fval[j]/dr[j]*(0-n[0]*n[1]-N[0]*N[1])-kij[j]*n[0]*n[1]
					diagsquare[0,2]+=fval[j]/dr[j]*(0-n[0]*n[2]-N[0]*N[2])-kij[j]*n[0]*n[2]
					# yx, yy and yz
					diagsquare[1,0]+=fval[j]/dr[j]*(0-n[1]*n[0]-N[1]*N[0])-kij[j]*n[1]*n[0]
					diagsquare[1,1]+=fval[j]/dr[j]*(1-n[1]*n[1]-N[1]*N[1])-kij[j]*n[1]*n[1]
					diagsquare[1,2]+=fval[j]/dr[j]*(0-n[1]*n[2]-N[1]*N[2])-kij[j]*n[1]*n[2]
					# zx, zy and zz
					diagsquare[2,0]+=fval[j]/dr[j]*(0-n[2]*n[0]-N[2]*N[0])-kij[j]*n[2]*n[0]
					diagsquare[2,1]+=fval[j]/dr[j]*(0-n[2]*n[1]-N[2]*N[1])-kij[j]*n[2]*n[1]
					diagsquare[2,2]+=fval[j]/dr[j]*(1-n[2]*n[2]-N[2]*N[2])-kij[j]*n[2]*n[2]
					# Add the curvature term if required
					# Sooooo. The derivation says there is just a -fval/Rval n N term here, due to the tilt of the normal with parallel transport
					# However, this term very explicitly punishes out-from-surface deviations. 
					# Making this larger gives much cleaner results. Soo - huh?
					if (addCurvature):
						# xx, xy and xz
						diagsquare[0,0]+=-fval[j]/Rval*n[0]*N[0]
						diagsquare[0,1]+=-fval[j]/Rval*n[1]*N[0]
						diagsquare[0,2]+=-fval[j]/Rval*n[2]*N[0]
						# yx, yy and yz
						diagsquare[1,0]+=-fval[j]/Rval*n[0]*N[1]
						diagsquare[1,1]+=-fval[j]/Rval*n[1]*N[1]
						diagsquare[1,2]+=-fval[j]/Rval*n[2]*N[1]
						# zx, zy and zz
						diagsquare[2,0]+=-fval[j]/Rval*n[0]*N[2]
						diagsquare[2,1]+=-fval[j]/Rval*n[1]*N[2]
						diagsquare[2,2]+=-fval[j]/Rval*n[2]*N[2]
					
					if (addRestoring):
						# Manual restoring force along the normal
						#ksurf=10
						diagsquare[0,0]+=-ksurf*N[0]*N[0]
						diagsquare[0,1]+=-ksurf*N[1]*N[0]
						diagsquare[0,2]+=-ksurf*N[2]*N[0]
						# yx, yy and yz
						diagsquare[1,0]+=-ksurf*N[0]*N[1]
						diagsquare[1,1]+=-ksurf*N[1]*N[1]
						diagsquare[1,2]+=-ksurf*N[2]*N[1]
						# zx, zy and zz
						diagsquare[2,0]+=-ksurf*N[0]*N[2]
						diagsquare[2,1]+=-ksurf*N[1]*N[2]
						diagsquare[2,2]+=-ksurf*N[2]*N[2]
				#print diagsquare
				self.Hessian[3*i:(3*i+3),3*i:(3*i+3)]=-diagsquare
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
		#HessianASym=0.5*(self.Hessian-np.transpose(self.Hessian))
		#print HessianASym
		# Use routines for hermitian eigenvector decomposition
		# Default is ascending order, which suits us
		print "Starting Diagonalisation!"
		self.eigval, self.eigvec = LA.eigh(HessianSym)
		if self.debug:
			# start with some debugging output
			plt.figure()
			eigrank=np.linspace(0,3*self.N,3*self.N)
			plt.plot(eigrank,self.eigval,'.-')
                wx=np.zeros((self.N,))
                wy=np.zeros((self.N,))
                wz=np.zeros((self.N,))
                for u in range(self.N):
                        # dimensional contributions
                        wx[u]=np.sum(self.eigvec[0:3*self.N:3,u]**2)
                        wy[u]=np.sum(self.eigvec[1:3*self.N:3,u]**2)
                        wz[u]=np.sum(self.eigvec[2:3*self.N:3,u]**2)
		print "The smallest eigenvalue is: " + str(np.amin(self.eigval))
		print self.eigval
		print wx
		print wy
		print wz
		
	def plotModes(self,omegamax=3.0,npts=100):
		# Straight here: The projection ratios on the sphere/plane
		projrat=np.zeros((3*self.N,))
		# Get the tangent bundle and coordinates
		# Note: this works for anything! In the plane, these are x, y, ex, ey
		theta,phi,etheta,ephi=self.geom.TangentBundle(self.rval)
		for u in range(3*self.N):
			thproj=etheta[:,0]*self.eigvec[0:3*self.N:3,u]+etheta[:,1]*self.eigvec[1:3*self.N:3,u]+etheta[:,2]*self.eigvec[2:3*self.N:3,u]
			phiproj=ephi[:,0]*self.eigvec[0:3*self.N:3,u]+ephi[:,1]*self.eigvec[1:3*self.N:3,u]+ephi[:,2]*self.eigvec[2:3*self.N:3,u]
			projrat[u]=np.sum(thproj**2)+np.sum(phiproj**2)	
		# For this histogram and only this histogram: replace negative eigenvalues by zeros
		eigplot=self.eigval
		badlist=list(np.where(eigplot<0.0)[0])
		eigplot[badlist]=0.0
		omega=np.real(np.sqrt(eigplot))
		ombin=np.linspace(0,omegamax,npts)
		dom=ombin[1]-ombin[0]
		omhist, bin_edges = np.histogram(omega,bins=ombin)
		omlabel=(np.round(omega/dom)).astype(int)
		projhist=np.zeros((npts-1,))
		projcount=np.zeros((npts-1,))
		for l in range(npts-1):
			pts=np.nonzero(omlabel==l)[0]
			projhist[l]+=sum(projrat[pts])
			projcount[l]+=len(pts)
		isdata=[index for index, value in enumerate(projcount) if value>0]
		projhist[isdata]/=projcount[isdata]
					
		plt.figure()
		plt.plot(ombin[1:]-dom/2,omhist,'o-k',label='DOS')
		plt.xlim(0,omegamax)
		plt.xlabel('omega')
		plt.ylabel('D(omega)')
		plt.title('Density of states')
		
		plt.figure()
		plt.plot(ombin[isdata]+dom/2,projhist[isdata],'o-r',label='Surface projection')
		plt.ylim(0,1.1)
		plt.xlim(0,omegamax)
		plt.xlabel('omega')
		plt.ylabel('projection')
		plt.title('Surface projection value')
		
		# Plotting eigenvectors (#0, #N and #2N)
		if self.geom.manifold=='plane':
			usepts=[0,1,2,3,self.N,2*self.N-1,3*self.N-1]
			for u in usepts:
				plt.figure()
				plt.quiver(self.rval[:,0],self.rval[:,1],self.eigvec[0:3*self.N:3,u],self.eigvec[1:3*self.N:3,u])
				# dimensional contributions
				wx=np.sum(self.eigvec[0:3*self.N:3,u]**2)
				wy=np.sum(self.eigvec[1:3*self.N:3,u]**2)
				wz=np.sum(self.eigvec[2:3*self.N:3,u]**2)
				plt.title('Eigenvector #' + str(u) + ' ' + str(wx+wy) + ' on plane')
				
		if self.geom.manifold=='sphere':
			# Get the tangent bundle and coordinates
			theta,phi,etheta,ephi=self.geom.TangentBundle(self.rval)
			usepts=[0,1,2,3,self.N,2*self.N-1,3*self.N-1]
			for u in usepts:
				fig = plt.figure()
				ax = fig.add_subplot(111, projection='3d')
				ax.quiver(self.rval[:,0], self.rval[:,1], self.rval[:,2], self.eigvec[0:3*self.N:3,u], self.eigvec[1:3*self.N:3,u], self.eigvec[2:3*self.N:3,u])
				# See how much of the mode is along the sphere here
				thproj=etheta[:,0]*self.eigvec[0:3*self.N:3,u]+etheta[:,1]*self.eigvec[1:3*self.N:3,u]+etheta[:,2]*self.eigvec[2:3*self.N:3,u]
				phiproj=ephi[:,0]*self.eigvec[0:3*self.N:3,u]+ephi[:,1]*self.eigvec[1:3*self.N:3,u]+ephi[:,2]*self.eigvec[2:3*self.N:3,u]
				frac=np.sum(thproj**2)+np.sum(phiproj**2)
				plt.title('Eigenvector #' + str(u) + ' ' + str(frac) +' on sphere') 
		
	
                        
	
		
		