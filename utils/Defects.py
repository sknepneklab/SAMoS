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
from Tesselation import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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
       
    # MISSING: getting the symmetry type directly from the parameter list, not as an input parameter
    # MISSING: an intelligent storage method for batch defects. Use same as cluster storage for friction project. 
    # However: necessites switch to binary pickle files, not text files for output
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
		