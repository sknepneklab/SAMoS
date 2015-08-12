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

from Geometry import *
from Configuration import *
try:
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
	HAS_MATPLOTLIB=True
except:
	HAS_MATPLOTLIB=False
	pass

# Profiles for the spherical case. 
# Warning: Implicitly spherical geometry. Results in other geometries will be garbage (should not crash though)
class Profiles:
	def __init__(self,conf,geom):
		self.conf=conf
		self.geom=geom
		#self.rval=self.conf.rval
		ez = np.array([0,0,1])  # lab frame z-axis
		# Simply get the axis as the mean crossproduct or r and v; assuming alignment. This should also not flip.
		self.direction=np.sum(np.cross(conf.rval,conf.vval),axis=0)
		self.orderpar=direction/len(conf.rval)
		print orderpar
		self.direction = self.direction/np.linalg.norm(self.direction)
		axis = np.cross(self.direction,ez)
		axis = axis/np.linalg.norm(axis)
		rot_angle = acos(np.dot(self.direction,ez))
		axis0 = np.empty(np.shape(self.conf.rval))
		axis0[:,0] = axis[0]
		axis0[:,1] = axis[1]
		axis0[:,2] = axis[2]
		self.conf.rotateFrame(axis0,-rot_angle)
		
		#rhat=((rval).transpose()/(np.sqrt(np.sum(rval**2,axis=1))).transpose()).transpose()
		#self.vel = np.sqrt(self.conf.vval[:,0]**2 + self.conf.vval[:,1]**2 + self.conf.vval[:,2]**2)
		velnorm=((self.conf.vval).transpose()/(vel).transpose()).transpose()
  
		self.theta,self.phi,self.etheta,self.ephi=geom.TangentBundle(self,rval)
		# Alpha, the angle between the local polarity and the equator; here represented by ephi
		self.alpha=-np.arcsin(np.sum(nval*etheta, axis=1))
		# Same thing for the velocity
		# No - add pi/2 to get something that does not add up to zero 
		#alpha_v=np.arccos(np.sum(velnorm*etheta, axis=1))
  
	def getProfiles(self,debug):
		eng, press,stress = conf.compute_energy_and_pressure()
		# Project the stresses into the e,theta,phi components. The rr component hast to be 0, and the r cross components
		# belong to the projection. So they are not all that interesting. 
		# We want the theta theta, theta phi, phi theta ant phi phi components (implicitly testing symmetries ...)
		# I give up on the notation. Stress is (N,3,3), the axes are (N,3). We want e_i sigma_ij e_j
		s_tt=np.sum(self.etheta*np.einsum('kij,kj->ki',stress,self.etheta),axis=1)
		s_tp=np.sum(self.etheta*np.einsum('...ij,...j->...i',stress,self.ephi),axis=1)
		s_pt=np.sum(self.ephi*np.einsum('...ij,...j->...i',stress,self.etheta),axis=1)
		s_pp=np.sum(self.ephi*np.einsum('...ij,...j->...i',stress,self.ephi),axis=1)
		
		# Setting up the binning. I changed this to go from -pi/2 to pi/2 consistently. This maybe makes less pretty pictures,
		# but the edges are going to be a lot cleaner. Also only one bin to handle accross multiple v0/J.
		# Can always rebin to less resolution if necessary
		# Position angle with the z axis
		theta_bin=np.linspace(0,np.pi,nbin+1)
		dtheta=theta_bin[1]-theta_bin[0]
		theta_out=theta_bin[:nbin]+dtheta/2-np.pi/2
		
		rho_profile, bin_edges = np.histogram(self.theta, bins=theta_bin,density=True)
		isdata=[index for index,value in enumerate(rho_profile) if (value >0)]
		normz=2*np.pi*self.conf.radius*abs(np.cos(theta_out))
		rho_profile[isdata]=rho_profile[isdata]/normz[isdata]
		rho_profile/=np.mean(rho_profile)
		vel_profile=np.zeros(np.shape(rho_profile))
		eng_profile=np.zeros(np.shape(rho_profile))
		press_profile=np.zeros(np.shape(rho_profile))
		s_tt_profile=np.zeros(np.shape(rho_profile))
		s_tp_profile=np.zeros(np.shape(rho_profile))
		s_pt_profile=np.zeros(np.shape(rho_profile))
		s_pp_profile=np.zeros(np.shape(rho_profile))
		alpha_profile=np.zeros(np.shape(rho_profile))
		alpha_v_profile=np.zeros(np.shape(rho_profile))
		for idx in range(nbin):
			inbin=[index for index,value in enumerate(self.theta) if (value >= theta_bin[idx]  and value<=theta_bin[idx+1])]
			#print len(inbin)
			if len(inbin)>0:
				vel_profile[idx]=np.mean(self.conf.vel[inbin])
				eng_profile[idx]=np.mean(eng[inbin])
				press_profile[idx]=np.mean(press[inbin])
				s_tt_profile[idx]=np.mean(s_tt[inbin])
				s_tp_profile[idx]=np.mean(s_tp[inbin])
				s_pt_profile[idx]=np.mean(s_pt[inbin])
				s_pp_profile[idx]=np.mean(s_pp[inbin])
				alpha_profile[idx]=np.mean(self.alpha[inbin])
				#alpha_v_profile[idx]=np.mean(alpha_v[inbin])
		
		# Debugging output
		if debug==True:
			if HAS_MATPLOTLIB:
				fig = plt.figure()
				ax = fig.add_subplot(111, projection='3d')
				ax.scatter(self.conf.rval[:,0], self.conf.rval[:,1], self.conf.rval[:,2], zdir='z', c='b')
			else:
				print 'Error: Matplotlib does not exist on this machine, cannot plot system'
			
		return [theta_out,rho_profile,vel_profile,eng_profile,press_profile,s_tt_profile,s_tp_profile,s_pt_profile,s_pp_profile,alpha_profile,alpha_v_profile,self.direction,self.orderpar]
        
