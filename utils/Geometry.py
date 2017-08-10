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

import numpy as np

# Default geometry class: 3 dimensional unconstrained space
class Geometry(object):
	def __init__(self,manifold,box,periodic,area=1.0):
		self.manifold=manifold
		self.periodic=periodic
		self.area=area
		self.Lx=box[0]
		self.Ly=box[1]
		self.Lz=box[2]
		print "Geometry: Created new geometry " + manifold + " for which periodic = " + str(self.periodic) + " and with box size " + str(self.Lx) + " x " +str(self.Ly) + " x " + str(self.Lz)
		
	def RotateMatrixVectorial(self,axis,theta):
		axlen=np.sqrt(axis[:,0]**2+axis[:,1]**2+axis[:,2]**2)
		#print axlen
		axis[:,0]=axis[:,0]/axlen
		axis[:,1]=axis[:,1]/axlen
		axis[:,2]=axis[:,2]/axlen
		a=np.cos(theta/2)
		b=-axis[:,0]*np.sin(theta/2)
		c=-axis[:,1]*np.sin(theta/2)
		d=-axis[:,2]*np.sin(theta/2)
		rotmat=np.empty((len(axis[:,0]),3,3))
		rotmat[:,0,0]=a*a+b*b-c*c-d*d
		rotmat[:,0,1]=2*(b*c-a*d)
		rotmat[:,0,2]=2*(b*d+a*c)
		rotmat[:,1,0]=2*(b*c+a*d)
		rotmat[:,1,1]=a*a+c*c-b*b-d*d
		rotmat[:,1,2]=2*(c*d-a*b)
		rotmat[:,2,0]=2*(b*d-a*c)
		rotmat[:,2,1]=2*(c*d+a*b)
		rotmat[:,2,2]=a*a+d*d-b*b-c*c
		return rotmat

	def RotateVectorial(self,v,n,phi):
		vrot=np.empty(np.shape(v))
		np.shape(vrot)
		rotmat=self.RotateMatrixVectorial(n,phi)
		np.shape(rotmat)
		vrot[:,0]=rotmat[:,0,0]*v[:,0]+rotmat[:,0,1]*v[:,1]+rotmat[:,0,2]*v[:,2]
		vrot[:,1]=rotmat[:,1,0]*v[:,0]+rotmat[:,1,1]*v[:,1]+rotmat[:,1,2]*v[:,2]
		vrot[:,2]=rotmat[:,2,0]*v[:,0]+rotmat[:,2,1]*v[:,1]+rotmat[:,2,2]*v[:,2]
		return vrot

	# Default: parallel transport just transports parallel; i.e. it does nothing in a flat geometry    
	def ParallelTransport(self,r1,r2,a2):
		return a2
		
	def ParallelTransportSingle(self,r1,r2,a2):
		return a2
		
	def GeodesicDistance1d(self,r1,r2):
		return np.sqrt(np.sum((r2-r1)**2,axis=0))
	
	def GeodesicDistance2d(self,r1,r2):
		return np.sqrt(np.sum((r2-r1)**2,axis=1))
	
		
	      
	# Just the cartesian distance in the plane, modulo periodic boundary conditions
	# Problem true to type right now ...
	# assume the first is a vector, the second a scalar
	# NEED TO GENERALIZE
	def ApplyPeriodic11(self,r1,r2):
		dr=r2-r1
		dr[0]-=self.Lx*np.round(dr[0]/self.Lx)
		dr[1]-=self.Ly*np.round(dr[1]/self.Ly)
		dr[2]-=self.Lz*np.round(dr[2]/self.Lz)
		return dr
	
	def ApplyPeriodic12(self,r1,r2):
		dr=r2-r1
		dr[:,0]-=self.Lx*np.round(dr[:,0]/self.Lx)
		dr[:,1]-=self.Ly*np.round(dr[:,1]/self.Ly)
		dr[:,2]-=self.Lz*np.round(dr[:,2]/self.Lz)
		return dr
		
	def ApplyPeriodic33(self,r1,r2):
		dr=r2-r1
		dr[:,:,0]-=self.Lx*np.round(dr[:,:,0]/self.Lx)
		dr[:,:,1]-=self.Ly*np.round(dr[:,:,1]/self.Ly)
		dr[:,:,2]-=self.Lz*np.round(dr[:,:,2]/self.Lz)
		return dr
	 
	# rather use these ones, kill the other ones eventually
	def ApplyPeriodic1d(self,dr):
		dr[0]-=self.Lx*np.round(dr[0]/self.Lx)
		dr[1]-=self.Ly*np.round(dr[1]/self.Ly)
		dr[2]-=self.Lz*np.round(dr[2]/self.Lz)
		return dr
	
	def ApplyPeriodic2d(self,dr):
		dr[:,0]-=self.Lx*np.round(dr[:,0]/self.Lx)
		dr[:,1]-=self.Ly*np.round(dr[:,1]/self.Ly)
		dr[:,2]-=self.Lz*np.round(dr[:,2]/self.Lz)
		return dr
		
	def ApplyPeriodic3d(self,dr):
		dr[:,:,0]-=self.Lx*np.round(dr[:,:,0]/self.Lx)
		dr[:,:,1]-=self.Ly*np.round(dr[:,:,1]/self.Ly)
		dr[:,:,2]-=self.Lz*np.round(dr[:,:,2]/self.Lz)
		return dr
	  
	# Sigh, looks like we need some specific x, y and z ones as well ...
	def ApplyPeriodicX(self,dr):
		dr-=self.Lx*np.round(dr/self.Lx)
		return dr
	
	def ApplyPeriodicY(self,dr):
		dr-=self.Ly*np.round(dr/self.Ly)
		return dr
	  
	def ApplyPeriodicZ(self,dr):
		dr-=self.Lz*np.round(dr/self.Lz)
		return dr
		

        
# Spherical geometry   
class GeometrySphere(Geometry):
	def __init__(self,param):
		self.R=param.r
		self.periodic=False
		self.area=4.0*np.pi*self.R**2
		print "GeometrySphere: Created new geometry sphere with radius " + str(self.R)
		super(GeometrySphere,self).__init__('sphere',param.box,self.periodic,self.area)
		
	# Fully vectorial version of parallel transport
	# 1.determine the cross product of the origins
	# 2.compute the magnitude of all the origin and cross vectors
	# 3.Compute the dot product of the origins
	# 4.The rotation axis is the direction of the cross product
	# 5.The rotation angle is the angle between the origin vectors, extracted from the dot product
	def ParallelTransport(self,r1,r2,a2):
		r2_x_r1=np.cross(r2,r1)
		#len_r2_x_r1=np.sqrt(r2_x_r1[:,0]**2+r2_x_r1[:,1]**2+r2_x_r1[:,2]**2)
		len_r2_x_r1=np.sqrt(np.sum(r2_x_r1**2,axis=1)) 
		#lenr1=np.sqrt(r1[:,0]**2+r1[:,1]**2+r1[:,2]**2)
		lenr1=np.sqrt(np.sum(r1**2,axis=1))
		#lenr2=np.sqrt(r2[:,0]**2+r2[:,1]**2+r2[:,2]**2)
		lenr2=np.sqrt(np.sum(r2**2,axis=1))
		dot_r1r2=r1[:,0]*r2[:,0]+r1[:,1]*r2[:,1]+r1[:,2]*r2[:,2]
		n=np.empty(np.shape(r1))
		n = r2_x_r1/len_r2_x_r1
		#n[:,0] = r2_x_r1[:,0]/len_r2_x_r1
		#n[:,1] = r2_x_r1[:,1]/len_r2_x_r1
		#n[:,2] = r2_x_r1[:,2]/len_r2_x_r1
		phi = np.arccos(dot_r1r2/(lenr1*lenr2))
		a2trans=rotate_vectorial(a2,n,-phi)
		return a2trans
		
	# same thing for one vector and a set (i.e. a particle and its neigbours)
	def ParallelTransportSingle(self,r1,r2,a2):
		r2_x_r1=np.cross(r2,r1)
		len_r2_x_r1=np.sqrt(np.sum(r2_x_r1**2,axis=1)) 
		lenr1=np.sqrt(np.sum(r1**2,axis=1))
		lenr2=np.sqrt(np.sum(r2**2,axis=1))
		dot_r1r2=np.dot(r1,r2)
		n=np.empty(np.shape(r1))
		n = r2_x_r1/len_r2_x_r1
		phi = np.arccos(dot_r1r2/(lenr1*lenr2))
		a2trans=rotate_vectorial(a2,n,-phi)
		return a2trans
		
	def GeodesicDistance1d(self,r1,r2):
		# That is the arc length for a sphere. Take the dot product, that's the cosine.
		# Invert with arccos and multiply by 2\pi R
		cosval=np.sum(r1*r2)/(np.sqrt(np.sum(r1**2))*np.sqrt(np.sum(r2**2)))
		return np.abs(np.arccos(cosval))*self.R
		
	def GeodesicDistance2d(self,r1,r2):
		# That is the arc length for a sphere. Take the dot product, that's the cosine.
		# Invert with arccos and multiply by 2\pi R
		cosval=np.sum(r1*r2,axis=1)/(np.sqrt(np.sum(r1**2,axis=1))*np.sqrt(np.sum(r2**2,axis=1)))
		return np.abs(np.arccos(cosval))*self.R
		
	# Determine the Euler angles, essentially. Find theta and phi for each particle,
	def TangentBundle(self,rval):
		rhat=((rval).transpose()/(np.sqrt(np.sum(rval**2,axis=1))).transpose()).transpose()
		# Angle theta with the z axis. arccos is between 0 and pi, so that's ok already
		theta=np.arccos(rhat[:,2])
		# From the euler angles: rx = sin theta cos phi
		# Choosing correct quadrant through the sign of ry=sin theta sin phi
		phi=np.sign(rhat[:,1]/(np.sin(theta)))*np.arccos(rhat[:,0]/(np.sin(theta)))
		# The other two of our trio of local coordinate vectors
		etheta = np.empty(np.shape(rval))
		etheta[:,0]=np.cos(theta)*np.cos(phi)
		etheta[:,1]=np.cos(theta)*np.sin(phi)
		etheta[:,2]=-np.sin(theta)
		ephi=np.empty(np.shape(rval))
		ephi[:,0]=-np.sin(phi)
		ephi[:,1]=np.cos(phi)
		ephi[:,2]=0
		return theta,phi,etheta,ephi
	
	# Complementary: The unit normal
	def UnitNormal(self,rval):
		rs = np.sqrt(rval[:,0]**2 + rval[:,1]**2 + rval[:,2]**2)
		rhat=((rval).transpose()/(rs).transpose()).transpose()
		return rhat
            
        # Complementary: The unit normal
	def UnitNormal1d(self,rval):
		rs = np.sqrt(rval[0]**2 + rval[1]**2 + rval[2]**2)
		rhat=rval/rs
		return rhat
       
       
# Plane without periodic boundary conditions
class GeometryPlane(Geometry):
	def __init__(self,param):
		self.Lx=param.lx
		self.Ly=param.ly
		self.periodic=False
		self.area=self.Lx*self.Ly
		print "GeometryPlane: Created new geometry plane with Lx = " + str(self.Lx) + " and Ly = " +str(self.Ly)
		super(GeometryPlane,self).__init__('plane',param.box,self.periodic,self.area)
		
		
	def TangentBundle(self,rval):
		x=rval[:,0]
		y=rval[:,1]
		ex = np.empty(np.shape(rval))
		ex[:,0]=1.0*np.ones(len(rval))
		ex[:,1]=1.0*np.zeros(len(rval))
		ex[:,2]=1.0*np.zeros(len(rval))
		ey = np.empty(np.shape(rval))
		ey[:,0]=1.0*np.zeros(len(rval))
		ey[:,1]=1.0*np.ones(len(rval))
		ey[:,2]=1.0*np.zeros(len(rval))
		return x,y,ex,ey
	
	# Unit normal: just the z direction
	def UnitNormal(self,rval):
		nvec=np.zeros(np.shape(rval))
		nvec[:,2]=1.0
		return nvec
		
		
	def GeodesicDistance12(self,r1,r2):
		dr=r2-r1
		return np.sqrt(dr[:,0]**2+dr[:,1]**2)
	
	def GeodesicDistance11(self,r1,r2):
		dr=r2-r1
		return np.sqrt(dr[0]**2+dr[1]**2)
	      
# Plane with periodic boundary conditions. By default, the plane is along x and y
class GeometryPeriodicPlane(Geometry):
	def __init__(self,param):
		self.Lx=param.lx
		self.Ly=param.ly
		self.periodic=True
		self.area=self.Lx*self.Ly
		print "GeometryPeriodicPlane: Created new geometry periodic plane with Lx = " + str(self.Lx) + " and Ly = " +str(self.Ly)
		super(GeometryPeriodicPlane,self).__init__('plane',param.box,self.periodic,self.area)
		
		
	def TangentBundle(self,rval):
		x=rval[:,0]
		y=rval[:,1]
		ex = np.empty(np.shape(rval))
		ex[:,0]=1.0*np.ones(len(rval))
		ex[:,1]=1.0*np.zeros(len(rval))
		ex[:,2]=1.0*np.zeros(len(rval))
		ey = np.empty(np.shape(rval))
		ey[:,0]=1.0*np.zeros(len(rval))
		ey[:,1]=1.0*np.ones(len(rval))
		ey[:,2]=1.0*np.zeros(len(rval))
		return x,y,ex,ey
	
	# Unit normal: just the z direction
	def UnitNormal(self,rval):
		nvec=np.zeros(np.shape(rval))
		nvec[:,2]=1.0
		return nvec
		
		
	def GeodesicDistance12(self,r1,r2):
		dr=self.ApplyPeriodic12(r1,r2)
		return np.sqrt(dr[:,0]**2+dr[:,1]**2)
	
	def GeodesicDistance11(self,r1,r2):
		dr=self.ApplyPeriodic11(r1,r2)
		return np.sqrt(dr[0]**2+dr[1]**2)

class GeometryTube(Geometry):
	def __init__(self,param):
		super(GeometryTube,self).__init__('tube',param.box,self.periodic)
		self.R=param.const_params['R']
		try:
			self.A=param.const_params['A']
		except KeyError:
			self.A = 0
		print "GeometryTube: Created new geometry tube with radius = " + str(self.R) + " and oscillation amplitude " + str(self.A)
		
		
		
class GeometryPeanut(Geometry):
	def __init__(self,param):	
		self.periodic=False
		super(GeometryPeanut,self).__init__('peanut',param.box,self.periodic)
		self.a=param.const_params['a']
		self.b=param.const_params['b']
		print "GeometryPeanut: Created new geometry peanut with parameters a = " + str(self.a) + " and b = " + str(self.b)
		# Determine the Euler angles, essentially. Find theta and phi for each particle,
		
	def TangentBundle(self,rval):
		x = rval[:,0]
		y = rval[:,1]
		z = rval[:,2]
		nx = (x**2 + y**2 + z**2 - self.a**2)*x
		ny = (x**2 + y**2 + z**2 + self.a**2)*y
		nz = (x**2 + y**2 + z**2 + self.a**2)*z
		n=np.array([nx,ny,nz])
		nhat=n/(np.sqrt(nx**2 + ny**2 + nz**2))
		# Angle theta with the z axis. arccos is between 0 and pi, so that's ok already
		phi=np.arctan2(y,z)
		# The other two of our trio of local coordinate vectors
		ephi = np.empty(np.shape(rval))
		ephi[:,0]=0
		ephi[:,1]=np.cos(phi)
		ephi[:,2]=-np.sin(phi)
		#((rval).transpose()/(rs).transpose()).transpose()
		ephi=(ephi.transpose()/(np.sqrt(ephi[:,0]**2 + ephi[:,1]**2 + ephi[:,2]**2).transpose())).transpose()
		splus=(x**2 + y**2 + z**2 + self.a**2)
		sminus=(x**2 + y**2 + z**2 - self.a**2)
		et = np.empty(np.shape(rval))
		et[:,0]=-np.sin(phi)*y*splus - np.cos(phi)*z*splus
		et[:,1]=np.sin(phi)*sminus*x + np.cos(phi)*z*splus
		et[:,2]=np.cos(phi)*sminus*x
		et=(et.transpose()/(np.sqrt(et[:,0]**2 + et[:,1]**2 + et[:,2]**2).transpose())).transpose()
		return x,phi,et,ephi
	
	# Complementary: The unit normal
	def UnitNormal(self,rval):
		x = rval[:,0]
		y = rval[:,1]
		z = rval[:,2]
		nx = (x**2 + y**2 + z**2 - self.a**2)*x
		ny = (x**2 + y**2 + z**2 + self.a**2)*y
		nz = (x**2 + y**2 + z**2 + self.a**2)*z
		n = np.array([nx, ny, nz])
		nmag = np.sqrt(nx**2 + ny**2 + nz**2)
		nhat=n/nmag
		return nhat.transpose()
        
		
class GeometryHourglass(Geometry):
	def __init__(self,param):
		if param.boxtype == 'periodic':
			self.periodic=True
		else:
			self.periodic=False
		print self.periodic
		self.R=param.const_params['R']
		self.A=param.const_params['A']
		self.H=param.box[2]
		self.n=1 # number of nodes - 1 for now
		print "GeometryHourglass: Created new geometry Hourglass with radius = " + str(self.R) + " and amplitude " +str(self.A)
		super(GeometryHourglass,self).__init__('hourglass',param.box,self.periodic)
		
	def TangentBundle(self,rval):
		x = rval[:,0]
		y = rval[:,1]
		z = rval[:,2]
		# Angle theta in the x y coordinates. 
		theta=np.arctan2(y,x)
		etheta = np.empty(np.shape(rval))
		etheta[:,0]=-np.sin(theta)
		etheta[:,1]=np.cos(theta)
		etheta[:,2]=0
		ez = np.empty(np.shape(rval))
		ez[:,0]=0
		ez[:,1]=0
		ez[:,2]=1
		return z,theta,ez,etheta
	
	# Complementary: The unit normal
	def UnitNormal(self,rval):
		axis = np.empty(np.shape(rval))
		axis[:,0]=2.0*rval[:,0]
		axis[:,1]=2.0*rval[:,1]
		axis[:,2]=-4.0*self.A*self.n*np.pi*np.cos((2.0*self.n*np.pi*rval[:,2])/self.H)*(self.R + self.A*np.sin((2.0*self.n*np.pi*rval[:,2])/self.H))/self.H;
		nlen = np.sqrt(np.sum(axis**2,axis=1))
		nhat = (axis.transpose()/nlen.transpose()).transpose()
		return nhat
		
		
	def GeodesicDistance12(self,r1,r2):
		dr=self.ApplyPeriodic12(r1,r2)
		return np.sqrt(dr[:,0]**2+dr[:,1]**2+dr[:,2]**2)
	
	def GeodesicDistance11(self,r1,r2):
		dr=self.ApplyPeriodic11(r1,r2)
		return np.sqrt(dr[0]**2+dr[1]**2+dr[2]**2)
		