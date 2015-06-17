# ################################################################
#  
#    Active Particles on Curved Spaces (APCS)
#    
#   Author: Silke Henkes
#   
#    ICSMB, Department of Physics
#    University of Aberdeen
#
#    Author: Rastko Sknepnek
#   
#    Division of Physics
#    School of Engineering, Physics and Mathematics
#    University of Dundee
#    
#    (c) 2013, 2014, 2015
#    
#    This program cannot be used, copied, or modified without
#    explicit permission of the author.
#  
# ################################################################

import numpy as np

# Default geometry class: 3 dimensional unconstrained space
class Geometry(object):
	def __init__(self,manifold,periodic):
		try:
			self.manifold=manifold
			self.periodic=periodic
			print "Created new geometry " + manifold + "for which periodic = " + self.periodic
		except:
			pass
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
		rotmat=rotate_matrix_vectorial(n,phi)
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
		
	def GeodesicDistance11(self,r1,r2):
		return self.GeodesicDistance(r1,r2)
	
	def GeodesicDistance21(self,r1,r2):
		return self.GeodesicDistance(r1,r2)
	# Default: just the cartesian distance
	def GeodesicDistance(self,r1,r2):
		return np.sqrt(np.sum((r2-r1)**2,axis=1))
        
# Spherical geometry   
class GeometrySphere(Geometry):
	def __init__(self,param):
		self.R=param.r
		self.periodic=False
		print "Created new geometry sphere with radius " + str(self.R)
		super(GeometrySphere,self).__init__('sphere',self.periodic)
		
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
        
# Plane with periodic boundary conditions. By default, the plane is along x and y
class GeometryPeriodicPlane(Geometry):
	def __init__(self,param):
		self.Lx=param.lx
		self.Ly=param.ly
		self.periodic=True
		print "Created new geometry periodic plane with Lx = " + str(self.Lx) + " and Ly = " +str(self.Ly)
		super(GeometryPeriodicPlane,self).__init__('plane',self.periodic)
		
		
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
		
	# Just the cartesian distance in the plane, modulo periodic boundary conditions
	# Problem true to type right now ...
	# assume the first is a vector, the second a scalar
	# NEED TO GENERALIZE
	def ApplyPeriodic11(self,r1,r2):
		dr=r2-r1
		dr[0]-=self.Lx*np.round(dr[0]/self.Lx)
		dr[1]-=self.Ly*np.round(dr[1]/self.Ly)
		return dr
	
	def ApplyPeriodic12(self,r1,r2):
		dr=r2-r1
		dr[:,0]-=self.Lx*np.round(dr[:,0]/self.Lx)
		dr[:,1]-=self.Ly*np.round(dr[:,1]/self.Ly)
		return dr
		
	def GeodesicDistance12(self,r1,r2):
		dr=self.ApplyPeriodic12(r1,r2)
		return np.sqrt(dr[:,0]**2+dr[:,1]**2)
	
	def GeodesicDistance11(self,r1,r2):
		dr=self.ApplyPeriodic11(r1,r2)
		return np.sqrt(dr[0]**2+dr[1]**2)

class GeometryTube(Geometry):
	def __init__(self,param):
		self.R=param.const_param['r']
		#self.periodic=False
		print "Created new geometry tube with radius = " + str(self.R)
		super(GeometryTube,self).__init__('tube')
		print "ERROR: Tube geometry has not yet been implemented! Geometry will default to 3d space."
		
class GeometryPeanut(Geometry):
	def __init__(self,param):
		super(GeometryPeanut,self).__init__('peanut')
		print "ERROR: Peanut geometry has not yet been implemented! Geometry will default to 3d space."
		
class GeometryHourglass(Geometry):
	def __init__(self,param):
		self.R=param.const_param['R']
		self.A=param.const_param['A']
		print "Created new geometry tube with radius = " + str(self.R) + " and amplitude " +str(self.A)
		super(GeometryHourglass,self).__init__('hourglass')
		print "ERROR: Hourglass geometry has not yet been implemented! Geometry will default to 3d space."
		