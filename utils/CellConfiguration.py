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
#from read_param_json import *
from read_param import *
from read_data import *
from read_faces import *
from CellList import *
from Interaction import *

try:
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
	HAS_MATPLOTLIB=True
except:
	HAS_MATPLOTLIB=False
	pass

class CellConfiguration:
	def __init__(self,param,filename_cells,filename_faces,ignore=False,debug=False):
		self.param=param
		# Read the local data
		geometries={'sphere':GeometrySphere,'plane':GeometryPlane,'plane_periodic':GeometryPeriodicPlane,'none':Geometry,'tube':GeometryTube,'peanut':GeometryPeanut,'hourglass':GeometryHourglass}
		print "Processing file : ", filename_cells
		data = ReadData(filename_cells)
		if data.keys.has_key('x'):
			x, y, z = np.array(data.data[data.keys['x']]), np.array(data.data[data.keys['y']]), np.array(data.data[data.keys['z']])
		else:
			print "Error: did not find positions in data file!"
			return 1
		self.N=len(x)
		if data.keys.has_key('vx'):
			vx, vy, vz = np.array(data.data[data.keys['vx']]), np.array(data.data[data.keys['vy']]), np.array(data.data[data.keys['vz']])
		else:
			vx, vy, vz = np.zeros(np.shape(x)),np.zeros(np.shape(y)),np.zeros(np.shape(z))
		try:
			nx, ny, nz = np.array(data.data[data.keys['nx']]), np.array(data.data[data.keys['ny']]), np.array(data.data[data.keys['nz']])
		except KeyError:
			nx, ny, nz = np.zeros(np.shape(x)),np.zeros(np.shape(y)),np.zeros(np.shape(z))
		if data.keys.has_key('fx'):
			fx, fy, fz = np.array(data.data[data.keys['fx']]), np.array(data.data[data.keys['fy']]), np.array(data.data[data.keys['fz']])
		# Now some cell-specific things: 
		# area  cell_area  cell_perim  cont_num  boundary 
		if data.keys.has_key('area'):
			self.area_native=np.array(data.data[data.keys['area']])
		else:
			# Assume the default of pi
			self.area_native=3.141592*np.ones(np.shape(x))
		if data.keys.has_key('cell_area'):
			self.area=np.array(data.data[data.keys['cell_area']])
		if data.keys.has_key('cell_perim'):
			self.perim=np.array(data.data[data.keys['cell_perim']])
		if data.keys.has_key('cont_num'):
			self.ncon=np.array(data.data[data.keys['cont_num']])
		if data.keys.has_key('boundary'):
			self.boundary=np.array(data.data[data.keys['boundary']])
		if data.keys.has_key('type'):
			self.ptype = data.data[data.keys['type']]
		else:
			self.ptype = np.ones((self.N,))
		if data.keys.has_key('flag'):
			self.flag = data.data[data.keys['flag']]
		self.rval = np.column_stack((x,y,z))
		self.vval = np.column_stack((vx,vy,vz))
		self.nval = np.column_stack((nx,ny,nz))
		self.fval = np.column_stack((fx,fy,fz))
		# Create the right geometry environment (TBC):
		self.geom=geometries[param.constraint](param)
		print self.geom
		# Create the Interaction class
		# Not yet created for the active vertex model
		#self.inter=Interaction(self.param,self.radius,ignore)
		
		# Again, not implemented yet. Keep it however, since this should eventually happen
		if self.geom.periodic:
			# Apply periodic geomtry conditions just in case (there seem to be some rounding errors floating around)
			self.rval=self.geom.ApplyPeriodic2d(self.rval)
			self.rval=self.geom.ApplyPeriodic12(np.array([0.0,0.0,0.0]),self.rval)
		# unit normal to the surface (only sphere so far)
		vel = np.sqrt(self.vval[:,0]**2 + self.vval[:,1]**2 + self.vval[:,2]**2)
		self.vhat=((self.vval).transpose()/(vel).transpose()).transpose()
		
		# Now get the connectivity here from the faces files
		#print "Processing file : ", filename_faces
		data_faces = ReadFaces(filename_faces)
		self.Faces=data_faces.Faces
		self.NFaces=data_faces.Nfaces
		
		# Create the Interaction class
		# This is essentially dummy for now
		self.radius=1
		self.inter=Interaction(self.param,self.radius,True)
		
		# Create the cell list
		# This is provisional, again. Here a lot more info from the nlist / triangulation should come in
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
        
	# For now: just get a couple of vital statistics out of this
	
	def getStatsCells(self,areabin,ratbin,conbin,maskradius=25):
		rdist=np.sqrt(self.rval[:,0]**2+self.rval[:,1]**2+self.rval[:,2]**2)
		# First find the ones that are not actually boundary
		inside=[index for index,value in enumerate(self.boundary) if value==0]
		#inside=[index for index,value in enumerate(rdist) if value<maskradius]
		#print maskradius
		Ninside=len(inside)
		print Ninside
		vel2 = self.vval[inside,0]**2 + self.vval[inside,1]**2 + self.vval[inside,2]**2
		vel2av=np.mean(vel2)
		# mean square force
		f2 = self.fval[inside,0]**2 + self.fval[inside,1]**2 + self.fval[inside,2]**2
		f2av=np.mean(f2)
		# Actual area
		areav=np.mean(self.area[inside])
		areadist,bins=np.histogram(self.area[inside],areabin,density=True)
		# The infamous perimeter / area ratio
		pratio=self.perim[inside]/np.sqrt(self.area[inside])
		# mean and distribution
		pratav=np.mean(pratio)
		pratdist,bins=np.histogram(pratio,ratbin,density=True)
		# contact number distribution (be careful with binning)
		# Well, current data does't have that keyword ..
		zav=0.0
		zdist=np.zeros((len(conbin)-1,))
		#zav=np.mean(self.ncon[inside])
		#zdist,bins=np.histogram(self.ncon[inside],conbin,density=True)
		# Boundary length as a fraction of total number of particles (useless for fixed)
		# Check that read as an int
		boundary=[index for index,value in enumerate(self.boundary) if value==1]
		bfrac=len(boundary)/self.N
		# Now get the actual length of the boundary ...
		borderface=self.Faces[self.NFaces-2]
		# Calculate its length, the old fashioned way
		bshift=[borderface[k-1] for k in range(len(borderface))]
                blens=np.sqrt(np.sum((self.rval[bshift,:]-self.rval[borderface,:])**2,axis=1))
                borderlen=np.sum(blens)
		return vel2av, f2av, areav, areadist, pratav,pratdist,zav,zdist,borderlen,bfrac,Ninside
	
	
	  
	
	