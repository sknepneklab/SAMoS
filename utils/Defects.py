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
from Tesselation import *
try:
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
	HAS_MATPLOTLIB=True
except:
	HAS_MATPLOTLIB=False
	pass

class Defects:
	def __init__(self,tess,conf,verbose=False):
		self.LoopList=tess.LoopList
		self.conf=conf
		self.normal=conf.geom.UnitNormal(conf.rval)
		self.theta,self.phi,self.etheta,self.ephi = conf.geom.TangentBundle(conf.rval)
		self.numdefect_n=0
		self.numdefect_v=0
		# Defect storage, up to 100
		# For n and velocity
		self.defects_n=[]
		self.defects_v=[]
		if self.conf.geom.manifold == 'sphere':
                        self.THETAMIN=2*self.conf.inter.sigma/self.conf.geom.R
                        if verbose:
                            print "Angle cutoff for polar defects: " + str(self.THETAMIN)
                        #fig = plt.figure()
                        #self.ax = fig.gca(projection='3d')
                        #self.ax.quiver(self.conf.rval[:,0], self.conf.rval[:,1], self.conf.rval[:,2], self.etheta[:,0], self.etheta[:,1], self.etheta[:,2])
                        #self.ax.quiver(self.conf.rval[:,0], self.conf.rval[:,1], self.conf.rval[:,2], self.ephi[:,0], self.ephi[:,1], self.ephi[:,2], color='r')
                        #plt.plot(self.conf.rval[:,0],self.phi,'.')
                        #plt.show()
		self.verbose=verbose
       
    # MISSING: getting the symmetry type directly from the parameter list, not as an input parameter
    # MISSING: an intelligent storage method for batch defects. Use same as cluster storage for friction project. 
    # However: necessites switch to binary pickle files, not text files for output
	def getDefects(self,symtype,field): 
		# Generalized algorithm for defects of any type
		# Count the defect charge. Times two, to use integers and easier if statements
		# Need to reinitialize defects in case multiple trackings are done on the same system
		if field == 'orientation':
                        print "Looking for orientation field defects!"
                        self.numdefect_n=0
                        self.defects_n=[]
                        char_n=0
                elif field == 'velocity':
                        print "Looking for velocity field defects!"
                        self.numdefect_v=0
                        self.defects_v=[]
                        char_v=0
                else:
                        print "Unknown symmetry type, aborting!"
                        return 1
		if symtype=='oldnematic':
			print "Tracking nematic defects with the Goldenfeld algorithm!"
		elif symtype=='polar':
			print "Tracking polar defects!"
		elif symtype=='nematic':
			print "Tracking nematic defects!"
		else:
			print "Unknown alignment symmetry type! Not tracking defects!"
			return 1
		for u in range(len(self.LoopList)):
			thisLoop=self.LoopList[u]
			if symtype=='oldnematic':
				charge=self.getDefectsGoldenfeld(thisLoop,field)
			else:
				charge=self.computeDefect(thisLoop,field,symtype)
			if abs(charge)>0:
				print "Found Defect in " + field + " field!"
				print charge
				if field == 'orientation':
                                        char_n+=charge
                                else:
                                        char_v+=charge
				# Construct the geometric centre of the defect
				r0s=self.conf.rval[thisLoop]
				if self.conf.geom.periodic:
					r0s[1:,:]=self.conf.geom.ApplyPeriodic12(r0s[0,:],r0s[1:,:])+r0s[0,:]
				rmval=np.mean(r0s,axis=0)
				if self.conf.geom.manifold=='sphere':
					rabs=np.sqrt(rmval[0]**2+rmval[1]**2+rmval[2]**2)
					rmval=rmval/rabs*self.conf.geom.R
				# Charge of the defect
				defbit=[charge]
				defbit.extend(list(rmval))
				if field == 'orientation':
                                        self.defects_n.append(defbit)
                                        self.numdefect_n+=1
                                else:
                                        self.defects_v.append(defbit)
                                        self.numdefect_v+=1
                if field == 'orientation':
                        print 'Number of orientation field defects: ' + str(self.numdefect_n)
                        print 'Total charge of orientation field defects: ' + str(char_n)
                        return self.defects_n, self.numdefect_n
                else:
                        print 'Number of velocity field defects: ' + str(self.numdefect_v)
                        print 'Total charge of velocity field defects: ' + str(char_v)
                        return self.defects_v,self.numdefect_v
                                        
			
        
	def computeDefect(self,thisLoop,field,symtype): 
		# Generalized algorithm for defects of any type
		# Count the defect charge. Times two, to use integers and easier if statements
		# nval
		thetatot=0
		t0=thisLoop[-1]
		#print "Starting loop " + str(thisLoop) + " in region " + str(self.conf.rval[t0,:])
		thetadum=[]
                for t in thisLoop[0:len(thisLoop)]:
                        # This is the old version based on small loops
                        if field == 'orientation':
                                ctheta=np.dot(self.conf.nval[t,:],self.conf.nval[t0,:])
                                stheta=np.dot(self.normal[t,:],np.cross(self.conf.nval[t,:],self.conf.nval[t0,:]))
                        else:
                                ctheta=np.dot(self.conf.vhat[t,:],self.conf.vhat[t0,:])
                                stheta=np.dot(self.normal[t,:],np.cross(self.conf.vhat[t,:],self.conf.vhat[t0,:]))
                        # Parallel transport does not seem to reliably work
                        #if self.conf.geom.manifold == 'sphere':
                        ##if tmp == 'plane':
                                #if field == 'orientation':
                                        #ctheta=np.dot(self.conf.nval[t,:],self.conf.nval[t0,:])
                                        #stheta=np.dot(self.normal[t,:],np.cross(self.conf.nval[t,:],self.conf.nval[t0,:]))
                                #else:
                                        #ctheta=np.dot(self.conf.vhat[t,:],self.conf.vhat[t0,:])
                                        #stheta=np.dot(self.normal[t,:],np.cross(self.conf.vhat[t,:],self.conf.vhat[t0,:]))
                        #else:
                                ## Need to take into account parallel transport. This is easiest done by just expressing the relevant vectors
                                ## in the proper local coordinate system (associated to the local normal)
                                ## This has been tested for spheres so far
                                #if field == 'orientation':
                                        #s0=np.dot(self.conf.nval[t0,:],self.etheta[t0,:])
                                        #s1=np.dot(self.conf.nval[t,:],self.etheta[t,:])
                                        #c0=np.dot(self.conf.nval[t0,:],self.ephi[t0,:])
                                        #c1=np.dot(self.conf.nval[t,:],self.ephi[t,:])
                                #else:
                                        #s0=np.dot(self.conf.vhat[t0,:],self.etheta[t0,:])
                                        #s1=np.dot(self.conf.vhat[t,:],self.etheta[t,:])
                                        #c0=np.dot(self.conf.vhat[t0,:],self.ephi[t0,:])
                                        #c1=np.dot(self.conf.vhat[t,:],self.ephi[t,:])
                                #ctheta=c1*c0+s1*s0
                                #stheta=s1*c0-c1*s0
                        if abs(stheta)>1:
                                stheta=np.sign(stheta)
                        if abs(ctheta)>1:
                                ctheta=np.sign(ctheta)
                        if np.isnan(stheta):
                                stheta=0
                        if np.isnan(ctheta):
                                ctheta=1
                        if symtype == 'nematic':
                                theta=np.arcsin(stheta*np.sign(ctheta))
                        elif symtype == 'polar':
                                theta=np.arccos(ctheta)*np.sign(stheta)
                        else:
                                print "Unknown symmetry type, doing nothing"
                                theta=0.0
                        #print theta
                        thetatot+=theta
                        #print thetatot
                        t0=t
                #if abs(thetatot)>0.5:
                        #print thetatot
                if symtype == 'nematic':
                        charge=0.5*int(round(thetatot/(np.pi)))
                else:
                        charge=int(round(thetatot/(2*np.pi)))
		if self.verbose:
                        if abs(charge)>0.1:
                                print "Found potential defect:"
                                print charge
                                print self.theta[t]
                                print self.phi[t]
                                print self.conf.rval[t,:]
                                #if charge>0:
                                        #self.ax.scatter(self.conf.rval[t,0], self.conf.rval[t,1], self.conf.rval[t,2],marker='o',s=50,color='r')
                                #else:
                                        #self.ax.scatter(self.conf.rval[t,0], self.conf.rval[t,1], self.conf.rval[t,2],marker='o',s=50,color='g')
		# except at the poles: that stuff is garbage?
                if self.conf.geom.manifold == 'sphere':
                    if abs(charge)>0.1:
                        if self.theta[t]<self.THETAMIN:
                            charge=0
                            if self.verbose:
                                    print "nuked this " + field + " field defect"
                        if self.theta[t]>(np.pi-self.THETAMIN):
                            charge=0
                            if self.verbose:
                                    print "nuked this " + field + " field defect"
                return charge
            
            
        def getDefectsGoldenfeld(self,thisLoop,field): 
		# Should already be ordered counterclockwise
		# Following a version of the Goldenfeld algorithm, with nx,ny,nz as is playing the role of the order parameter. The sphere is in cartesian space
		# The old nematic algorithm, based on the hemispheres            
		# The polarization vector nval
		ctheta=1
		coord=[]
		coord.append(nval[thisLoop[0],:])
		for t in range(1,len(thisLoop)):
                        if field == 'orientation':
                                ctheta=np.dot(self.conf.nval[thisLoop[t],:],np.sign(ctheta)*self.conf.nval[thisLoop[t-1],:])
                                # Nematic: append the order parameter, rotated through the *smaller* angle
                                coord.append(np.sign(ctheta)*self.conf.nval[thisLoop[t],:])
                                # Find out if the last point and the starting point are in the same hemisphere. 
                        elif field == 'velocity':
                                ctheta=np.dot(self.conf.vhat[thisLoop[t],:],np.sign(ctheta)*self.conf.vhat[thisLoop[t-1],:])
				# Nematic: append the order parameter, rotated through the *smaller* angle
				coord.append(np.sign(ctheta)*self.conf.vhat[thisLoop[t],:])
				# Find out if the last point and the starting point are in the same hemisphere. 
                        else:
                                print "Unknown field, doing nothing"
                                return 0
		cdefect=np.dot(coord[t],coord[0])
		if cdefect<0:
			charge=0.5
		else:
			charge=0.0
		return charge
            
	def PlotDefects(self):
		if HAS_MATPLOTLIB:
			fig = plt.figure()
			ax = fig.add_subplot(111, projection='3d')
			print self.defects_n
			ax.scatter(self.conf.rval[:,0], self.conf.rval[:,1], self.conf.rval[:,2], zdir='z', c='b',s=4)
			ax.scatter(self.defects_n[:][1], self.defects_n[:][2], self.defects_n[:][3], zdir='z', c='r',s=50)
			ax.scatter(self.defects_v[:][1], self.defects_v[:][2], self.defects_v[:][3], zdir='z', c='g',s=50)
		else:
			print 'Error: Matplotlib does not exist on this machine, cannot plot system and defects'
		