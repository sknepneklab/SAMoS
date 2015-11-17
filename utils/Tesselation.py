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

from Configuration import *
from CellList import *

MMAX=3.0
class Tesselation:
    
	def __init__(self,conf,debug=False):
		self.conf=conf
		self.rval=self.conf.rval
		self.geom=self.conf.geom
		self.debug=debug
                self.ordered_patches = False
		
	def findLoop(self,closeHoles=False):
		neighList=[]
		self.Ival=[]
		self.Jval=[]
		if closeHoles:
			Inei=[[] for k in range(len(self.rval))]
		else:
			Inei=[]
		count=0
		# Identify all neighbours and add them to a list. Keep i->j and j->i separate
		# The label is in neighList, the particle numbers are in Ival and Jval
		#if self.conf.monodisperse:
		if self.conf.param.potential=='soft':
			dmax=2*self.conf.sigma
			print dmax
		elif self.conf.param.potential=='morse':
			#dmax=16*self.conf.sigma**2
			dmax=2*self.conf.sigma
		else:
			dmax=2*self.conf.sigma
			print "Warning: unimplemented potential, defaulting to maximum contact distance 2"
		print dmax	
		if self.conf.param.potential=='morse':
			re=self.conf.param.pot_params['re']
		#cl = CellList(2.0*np.sqrt(dmax),self.conf.param.box)
		#for i in range(len(self.rval)):
			#cl.add_vertex(self.rval[i,:],i)
		for i in range(len(self.rval)):
			neighbours=[]
			dist=self.geom.GeodesicDistance12(self.rval[i,:],self.rval)
			if closeHoles:
				mult=1.0
				while len(neighbours)<4 and mult<MMAX:
					if self.conf.monodisperse:
						neighbours=[index for index,value in enumerate(dist) if value <mult*dmax]
					else:
						if self.conf.param.potential=='soft':
							neighbours=[index for index,value in enumerate(dist) if value < mult*(self.conf.radius[i]+self.conf.radius[index])]
						elif self.conf.param.potential=='morse':
							neighbours=[index for index,value in enumerate(dist) if value < mult*(re*self.conf.radius[i]+re*self.conf.radius[index])]
						else:
							neighbours=[index for index,value in enumerate(dist) if value <0.8*mult*(self.conf.radius[i]+self.conf.radius[index])]	
					neighbours.remove(i)
					mult=1.1*mult
                                #print i, ' --> ', neighbours
				# find the new neighbours (because of our asymmetric contact algorithm there are neighbours which haven't been found
				# for a particle if it is far away from the others, and the others are in a high density region)
				jexisting=[]
				for n in Inei[i]:
					jexisting.append(self.Jval[n])
				neighs_new = [a for a in neighbours if a not in jexisting]
				neighList.extend([u for u in range(count,count+len(neighs_new))])
				self.Ival.extend([i for k in range(len(neighs_new))])
				self.Jval.extend(neighs_new)
				Inei[i].extend([u for u in range(count,count+len(neighs_new))])
				count+=len(neighs_new)
				# and the reverse contacts
				neighList.extend([u for u in range(count,count+len(neighs_new))])
				self.Jval.extend([i for k in range(len(neighs_new))])
				self.Ival.extend(neighs_new)
				# Need to find to which particles the other contacts of I (those <i) have been assigned
				for a in range(len(neighs_new)):
					Inei[neighs_new[a]].append(count+a)
				count+=len(neighs_new)
			else:
				if self.conf.monodisperse:
					neighbours=[index for index,value in enumerate(dist) if value <1.4*dmax]
				else:
					if self.conf.param.potential=='soft':
						neighbours=[index for index,value in enumerate(dist) if value <(self.conf.radius[i]+self.conf.radius[index])]		
					elif self.conf.param.potential=='morse':
						neighbours=[index for index,value in enumerate(dist) if value <(re*self.conf.radius[i]+re*self.conf.radius[index])]
					else:
						neighbours=[index for index,value in enumerate(dist) if value <0.8*(self.conf.radius[i]+self.conf.radius[index])]	
				neighbours.remove(i)
				print len(neighbours)
				neighList.extend([u for u in range(count,count+len(neighbours))])
				self.Ival.extend([i for k in range(len(neighbours))])
				#self.Jval.extend(neighs[neighbours])
				self.Jval.extend(neighbours)
				Inei.append([u for u in range(count,count+len(neighbours))])
				count+=len(neighbours)
		# Identify loops based on the neighbour list. Kick out any (one-way) contacts that have occured so far
		Jarray=np.array(self.Jval)
		self.LoopList=[]
		# The dual: which loops belong to which particle
		self.ParList=[[] for k in range(len(self.rval))]
		self.LoopCen=[]
		self.l=0
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
				if self.conf.geom.periodic:
					dr0hat=self.geom.ApplyPeriodic11(self.rval[self.Ival[idx],:],self.rval[self.Jval[idx],:])
				else:
					dr0hat=self.rval[self.Jval[idx],:]-self.rval[self.Ival[idx],:]
				dr0hat/=np.sqrt(np.sum(dr0hat**2))
				jnei0=Inei[self.Jval[idx]]
				jnei=list(Jarray[jnei0])  
				if self.conf.geom.periodic:
					drvec=self.geom.ApplyPeriodic12(self.rval[self.Jval[idx],:],self.rval[jnei,:])
				else:
					drvec=self.rval[jnei,:]-self.rval[self.Jval[idx],:]
				drhat=((drvec).transpose()/(np.sqrt(np.sum(drvec**2,axis=1))).transpose()).transpose()
				cbeta=np.einsum('kj,j->k',drhat,self.conf.e2[self.Jval[idx],:])
				sbeta=np.einsum('kj,j->k',drhat,self.conf.e1[self.Jval[idx],:])
				cbeta0=np.dot(dr0hat,self.conf.e2[self.Jval[idx],:])
				sbeta0=np.dot(dr0hat,self.conf.e1[self.Jval[idx],:])
			
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
					self.ParList[Jarray[idx]].append(self.l)
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
			looppos=self.rval[llist]
			# Apply periodic boundary conditions here too ...
			# Arbitrary decision - go with the side of the first one
			if self.geom.periodic:
				dl=self.geom.ApplyPeriodic12(looppos[0,:],looppos)
				looppos=looppos[0,:]+dl
			lcen=[np.mean(looppos[:,0]), np.mean(looppos[:,1]),np.mean(looppos[:,2])]
			if self.conf.geom.manifold=='sphere':
				llen=np.sqrt(lcen[0]**2+lcen[1]**2+lcen[2]**2)
				lcen/=llen
				lcen*=self.conf.geom.R
			self.LoopCen.append(lcen)
			self.LoopList.append(llist)
			self.l+=1
		return self.LoopList,self.Ival,self.Jval
      
	# Much prettier: a loop that is too big (as measured by the mean square distance of the distances to the particles)
	# Deconstruct it into lots of little loops (virtual ones), with defined centers
	def makeEdges(self,rmax):   
		for l0 in range(len(self.LoopList)):
			llist=self.LoopList[l0]
			looppos=self.rval[llist]
			if self.geom.periodic:
				dlvec=self.geom.ApplyPeriodic12(self.LoopCen[l0],looppos)
			else:
				dlvec=looppos-self.LoopCen[l0]
			isLong=np.sqrt(np.sum(np.sum(dlvec**2,axis=1)))/len(llist)
			#if len(llist)>5:
				#print llist
				#print isLong
			#if isLong>rmax:
			if len(llist)>20:
				print "Loop " + str(l0) + " with particles " + str(llist) + " is too big! "
				for k in range(len(llist)):
					kside=k-1
					if kside<0:
						kside=len(llist)-1
					# Attempting to catch the inward pointing loops: the have to be global boundary ~sqrt(N)
					if len(llist)<0.5*np.sqrt(len(self.rval)):
						newcen=0.5*(self.rval[llist[k]]+self.rval[llist[kside]])-self.conf.sigma*dlvec[k,:]/np.sqrt(np.sum(dlvec[k,:]**2))
					else:
						newcen=0.5*(self.rval[llist[k]]+self.rval[llist[kside]])+self.conf.sigma*dlvec[k,:]/np.sqrt(np.sum(dlvec[k,:]**2))
					self.LoopCen.append(newcen)
					try:
						self.ParList[llist[k]].remove(l0)
					except ValueError:
						pass
					self.ParList[llist[k]].append(self.l)
					try:
						self.ParList[llist[kside]].remove(l0)
					except ValueError:
						pass
					self.ParList[llist[kside]].append(self.l)
					self.l+=1
        
	# While we are at it, we can construct the dual tesselation here.
	# All that's missing is to order the patches for the particles counterclockwise
	def OrderPatches(self):
		LoopCen1=np.array(self.LoopCen)
		for i in range(len(self.rval)):
			parray=np.array(self.ParList[i])
			#print parray
			if self.geom.periodic:
				drvec=self.geom.ApplyPeriodic12(self.rval[i,:],LoopCen1[self.ParList[i]])
			else:
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
			self.ParList[i]=parray[lorder] 
			# Use the new ParList structure where loops belong to particles are stored
		self.ordered_patches = True
                return self.LoopList,self.LoopCen,self.ParList,self.Ival,self.Jval


	def ComputePatchArea(self):
		if not self.ordered_patches: 
			raise Exception('Patches have to be ordered in order to cumpute their area.')
		self.area = []
		for k in xrange(len(self.ParList)):
			if len(self.ParList[k])>0:
				xc, yc, zc = 0.0, 0.0, 0.0
				for l in self.ParList[k]:
					xc += self.LoopCen[l][0]
					yc += self.LoopCen[l][1]
					zc += self.LoopCen[l][2]
				xc /= len(self.ParList[k])
				yc /= len(self.ParList[k])
				zc /= len(self.ParList[k])
				N = len(self.ParList[k])
				area = 0.0
				for i in xrange(N):
					#print self.ParList[i]
					x1 = self.LoopCen[self.ParList[k][i]][0] - xc
					y1 = self.LoopCen[self.ParList[k][i]][1] - yc
					z1 = self.LoopCen[self.ParList[k][i]][2] - zc
					x2 = self.LoopCen[self.ParList[k][(i+1)%N]][0] - xc
					y2 = self.LoopCen[self.ParList[k][(i+1)%N]][1] - yc 
					z2 = self.LoopCen[self.ParList[k][(i+1)%N]][2] - zc
					n = np.cross([x1,y1,z1],[x2,y2,z2])
					area += 0.5*np.sqrt(np.dot(n,n))
			else:
				area=0.0
			self.area.append(area)
			
	def computeContractile(self,alpha):
		self.ComputePatchArea()
		N = len(self.rval)
		stress = np.zeros((N,3,3))
		pressure = np.zeros((N,))
		for i in range(N):
			# simple contractile part n x n
			#stress[i,:,:]=alpha*np.einsum('ij,kl->self.nval,self.nval)
			pressure[i]=alpha*self.conf.param.pot_params['k']*np.sum(self.conf.nval[i,:]**2)*self.area[i]
		return pressure

