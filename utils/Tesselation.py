# ################################################################
#  
#    Active Particles on Curved Spaces (APCS)
#    
#    Author: Silke Henkes
#
#    ICSMB, Department of Physics
#    University of Aberdeen 
#  
#    (c) 2014, 2015
#    
#    This program cannot be used, copied, or modified without
#    explicit permission of the authors.
#  
# ################################################################

from Configuration import *
from CellList import *

class Tesselation:
    
	def __init__(self,conf):
		self.conf=conf
		self.rval=self.conf.rval
		
	def findLoop(self):
		neighList=[]
		self.Ival=[]
		self.Jval=[]
		Inei=[]
		count=0
		# Identify all neighbours and add them to a list. Keep i->j and j->i separate
		# The label is in neighList, the particle numbers are in Ival and Jval
		if self.conf.monodisperse:
			if self.conf.param.potential=='soft':
				dmax=4*self.conf.sigma**2
			elif self.conf.param.potential=='morse':
				dmax=16*self.conf.sigma**2
			else:
				dmax=4*self.conf.sigma**2
				print "Warning: unimplemented potential, defaulting to maximum contact distance 2"
		if self.conf.param.potential=='morse':
			re=self.conf.param.pot_param['re']
		cl = CellList(2.0*np.sqrt(dmax),self.conf.param.box)
		for i in range(len(self.rval)):
			cl.add_vertex(self.rval[i,:],i)
		for i in range(len(self.rval)):
			dist=np.sum((self.rval-self.rval[i,:])**2,axis=1)
			#neighs = np.array(cl.get_neighbours(self.rval[i,:]))
			#dist=self.conf.geom.GeodesicDistance(self.rval[neighs,:],self.rval[i,:])
			#dist=self.conf.geom.GeodesicDistance(self.rval[neighs,:],self.rval[i,:])
			if self.conf.monodisperse:		
				neighbours=[index for index,value in enumerate(dist) if value <dmax]
			else:
				if self.conf.param.potential=='soft':
					neighbours=[index for index,value in enumerate(dist) if value < (self.conf.radius[i]+self.conf.radius[index])**2]
				elif self.conf.param.potential=='morse':	
					neighbours=[index for index,value in enumerate(dist) if value < (re*self.conf.radius[i]+re*self.conf.radius[index])**2]
				else:
					neighbours=[index for index,value in enumerate(dist) if value < (self.conf.radius[i]+self.conf.radius[index])**2]
					print "Warning: unimplemented potential, defaulting to maximum contact distance 2"
			neighbours.remove(i)
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
			looppos=self.rval[llist]
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
			print parray
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