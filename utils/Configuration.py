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

from Geometry import *
from read_param import *
from read_data import *
from CellList import *
from Interaction import *


class Configuration:
	def __init__(self,param,filename,ignore=False,debug=False):
		self.param=param
		# Outdated list of geometries
		#geometries={'sphere':GeometrySphere,'plane':GeometryPeriodicPlane,'none':Geometry,'tube':GeometryTube,'peanut':GeometryPeanut,'hourglass':GeometryHourglass}
		geometries={'sphere':GeometrySphere,'plane':GeometryPlane,'plane_periodic':GeometryPeriodicPlane,'none':Geometry,'tube':GeometryTube,'peanut':GeometryPeanut,'hourglass':GeometryHourglass}
		print "Processing file : ", filename
		data = ReadData(filename)
		x, y, z = np.array(data.data[data.keys['x']]), np.array(data.data[data.keys['y']]), np.array(data.data[data.keys['z']])
		vx, vy, vz = np.array(data.data[data.keys['vx']]), np.array(data.data[data.keys['vy']]), np.array(data.data[data.keys['vz']])
		try:
			nx, ny, nz = np.array(data.data[data.keys['nx']]), np.array(data.data[data.keys['ny']]), np.array(data.data[data.keys['nz']])
		except KeyError:
			nx, ny, nz = np.zeros(np.shape(x)),np.zeros(np.shape(y)),np.zeros(np.shape(z))
		self.monodisperse=False
		self.N=len(x)
		if not data.keys.has_key('radius'): 
			# MISSING: read them in from the initial configuration
			self.radius = np.array([1.0 for i in range(self.N)])
			self.monodisperse=True
			#self.sigma=1.0
		else: 
			self.radius = np.array(data.data[data.keys['radius']])	
			#self.sigma = np.mean(self.radius)
		if data.keys.has_key('type'):
			self.ptype = data.data[data.keys['type']]
		else:
			self.ptype = np.ones((self.N,))
		if data.keys.has_key('flag'):
			self.flag = data.data[data.keys['flag']]
		self.rval = np.column_stack((x,y,z))
		self.vval = np.column_stack((vx,vy,vz))
		self.nval = np.column_stack((nx,ny,nz))
		# Create the right geometry environment (TBC):
		self.geom=geometries[param.constraint](param)
		print self.geom
		# Create the Interaction class
		self.inter=Interaction(self.param,self.radius,ignore)
		
		if self.geom.periodic:
			# Apply periodic geomtry conditions just in case (there seem to be some rounding errors floating around)
			self.rval=self.geom.ApplyPeriodic2d(self.rval)
			self.rval=self.geom.ApplyPeriodic12(np.array([0.0,0.0,0.0]),self.rval)
		# unit normal to the surface (only sphere so far)
		vel = np.sqrt(self.vval[:,0]**2 + self.vval[:,1]**2 + self.vval[:,2]**2)
		self.vhat=((self.vval).transpose()/(vel).transpose()).transpose()
		
		# Create the cell list
		cellsize=param.nlist_rcut
		if cellsize>5*self.inter.sigma:
			cellsize=5*self.inter.sigma
			print "Warning! Reduced the cell size to manageable proportions (5 times mean radius). Re-check if simulating very long objects!"
		self.clist=CellList(self.geom,cellsize)
		# Populate it with all the particles:
		for k in range(self.N):
			self.clist.add_particle(self.rval[k,:],k)
		#self.clist.printMe()
		
		if debug:
			fig = plt.figure()
			ax = fig.add_subplot(111, projection='3d')
			ax.scatter(self.rval[:,0], self.rval[:,1], self.rval[:,2], zdir='z', c='b')
		
	# Tangent bundle: Coordinates in an appropriate coordinate system defined on the manifold
	# and the coordinate vectors themselves, for all the particles
	def getTangentBundle(self):
		self.x1,self.x2,self.e1,self.e2=self.geom.TangentBundle(self.rval)
		return self.x1,self.x2,self.e1,self.e2
	
	def rotateFrame(self,axis,rot_angle):
		self.rval = self.geom.RotateVectorial(self.rval,axis,-rot_angle)
		self.vval = self.geom.RotateVectorial(self.vval,axis,-rot_angle)
		self.nval = self.geom.RotateVectorial(self.nval,axis,-rot_angle)
		self.nval=((self.nval).transpose()/(np.sqrt(np.sum(self.nval**2,axis=1))).transpose()).transpose()
		self.vel = np.sqrt(self.vval[:,0]**2 + self.vval[:,1]**2 + self.vval[:,2]**2)
		
	# Need to redo boxes after something like that
	def redoCellList(self):
		del self.clist
		# Create the cell list
		cellsize=self.param.nlist_rcut
		if cellsize>5*self.inter.sigma:
			cellsize=5*self.inter.sigma
			print "Warning! Reduced the cell size to manageable proportions (5 times mean radius). Re-check if simulating very long objects!"
		self.clist=CellList(self.geom,cellsize)
		# Populate it with all the particles:
		for k in range(self.N):
			self.clist.add_particle(self.rval[k,:],k)
		#self.clist.printMe()
        
	def getNeighbours(self,i,mult,dmax):
		# Find potential neighbours from neighbour list first
		cneighbours=self.clist.get_neighbours(self.rval[i,:])
		#print "Cell list neighbours: " + str(len(cneighbours))
		drvec0=self.geom.ApplyPeriodic2d(self.rval[cneighbours,:]-self.rval[i,:])
		dist=np.sqrt(drvec0[:,0]**2+drvec0[:,1]**2+drvec0[:,2]**2)
		#dist=self.geom.GeodesicDistance12(self.rval[cneighbours,:],self.rval[i,:])
		#print "Mean cutoff: " + str(mult*dmax)
		if self.monodisperse: 
			neighbours=[cneighbours[index] for index,value in enumerate(dist) if value <mult*dmax]
		else:
			neighbours=[cneighbours[index] for index,value in enumerate(dist) if value < mult*(self.radius[i]+self.radius[cneighbours[index]])]
		## Stupid one for debugging purposes:
		#dist=self.geom.GeodesicDistance12(self.rval,self.rval[i,:])
		#neighbours = [index for index, value in enumerate(dist) if value < mult*(self.radius[i]+self.radius[index])]
		neighbours.remove(i)
		#print "Contact neighbours: " + str(len(neighbours))
		#print neighbours
		drvec=self.geom.ApplyPeriodic2d(self.rval[neighbours,:]-self.rval[i,:])
		dr=np.sqrt(drvec[:,0]**2+drvec[:,1]**2+drvec[:,2]**2)
		return neighbours, drvec, dr
	      
	def compute_energy_and_pressure(self):
		eng = np.zeros(self.N)
		press = np.zeros(self.N)
		ncon = np.zeros(self.N)
		stress = np.zeros((self.N,3,3))
		for i in range(self.N):
			neighbours, drvec, dr=self.getNeighbours(i,self.inter.getMult(),self.inter.getDmax())
			ncon[i]=len(neighbours)
			eng[neighbours]+=self.inter.getEnergy(i,neighbours,drvec,dr)
			press_val,stress_val=self.inter.getStresses(i,neighbours,drvec,dr)
			stress[neighbours,:,:]+=stress_val
			press[neighbours]+=press_val
		return [eng, press, ncon,stress]
	
	
	# Flat case statistics (or other geometry statistics, if desired)
	def getStatsBasic(self,debug=False):
		vel2 = self.vval[:,0]**2 + self.vval[:,1]**2 + self.vval[:,2]**2
		vel2av=np.mean(vel2)
		phival=np.pi*np.sum(self.radius**2)/self.geom.area
		ndensity=self.N/self.geom.area
		eng, press,ncon,stress = self.compute_energy_and_pressure()
		pressure=np.sum(press)/self.geom.area
		fmoment=np.mean(press)
		energy=np.mean(eng)
		energytot=np.sum(eng)
		zav=np.mean(ncon)
		return vel2av, phival,ndensity, pressure,fmoment,energy,energytot,zav
	  
	def getStatsBand(self,debug=False):
		ez = np.array([0,0,1])  # lab frame z-axis
		# The order parameter with v_0 still in it. Normalize in final polish
		orderparV=np.sum(vval,axis=0)/len(vval)
		orderpar=np.sum(nval,axis=0)/len(nval)
		print orderpar
		print orderparV
		direction = orderpar/np.linalg.norm(orderpar)
		directionV = orderparV/np.linalg.norm(orderparV)
		axisorth= np.cross(direction,directionV)
		axisval=np.linalg.norm(axisorth)
		alpha=np.arcsin(axisval)
		axisorth=axisorth/axisval
		axisnorm=np.cross(ez,directionV)
		axisnorm/=np.linalg.norm(axisnorm)
		
		print directionV
		print axisorth
		
		vel = np.sqrt(self.vval[:,0]**2 + self.vval[:,1]**2 + self.vval[:,2]**2)
		velnorm=((self.vval).transpose()/(vel).transpose()).transpose()
		
		eng, press,stress = self.compute_energy_and_pressure()
		print np.shape(stress)
		# Project the stresses into the e,theta,phi components. The rr component hast to be 0, and the r cross components
		# belong to the projection. So they are not all that interesting. 
		# We want the theta theta, theta phi, phi theta ant phi phi components (implicitly testing symmetries ...)
		# I give up on the notation. Stress is (N,3,3), the axes are (N,3). We want e_i sigma_ij e_j
		s_tt=np.sum(axisnorm*np.einsum('kij,j->ki',stress,axisnorm),axis=1)
		s_tp=np.sum(axisnorm*np.einsum('kij,j->ki',stress,directionV),axis=1)
		s_pt=np.sum(directionV*np.einsum('kij,j->ki',stress,axisnorm),axis=1)
		s_pp=np.sum(directionV*np.einsum('kij,j->ki',stress,directionV),axis=1)
		print np.shape(s_tt)
		# Mean density really makes no sense? Determined by the initial conditions in periodic boundary conditions.
		# I do not wish to set up artificial bins in a translationally invariant system
		vel_av=np.mean(vel)
		eng_av=np.mean(eng)
		press_av=np.mean(press)
		s_tt_av=np.mean(s_tt)
		s_tp_av=np.mean(s_tp)
		s_pt_av=np.mean(s_pt)
		s_pp_av=np.mean(s_pp)
		
		# Debugging output
		if debug==True:
			if HAS_MATPLOTLIB:
				fig = plt.figure()
				ax = fig.add_subplot(111, projection='3d')
				ax.scatter(rval[:,0], rval[:,1], rval[:,2], zdir='z', c='b')
			else:
				print 'Error: Matplotlib does not exist on this machine, cannot plot system'
			
		return [vel_av,eng_av,press_av,s_tt_av,s_tp_av,s_pt_av,s_pp_av,alpha,direction,directionV,orderpar,orderparV]
	