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

from read_param import *
import numpy as np

# Given importance of radii et al, we have to pass Interaction the information about
# radii (and positions)? Alternatively, pass it just the unchanging particle labels
class Interaction:

	def __init__(self,param,radius,ignore=False,debug=False):
		self.param=param
		self.radius=radius
		self.sigma= np.mean(self.radius)
		self.ignore=ignore
		self.debug=debug
		# First step: unpack parameters to find out what kind of types we have, and what kind of potentials
		if ((self.param.ntypes==1) or (self.ignore==True)):
			if self.param.potential=='soft':
				self.k=self.param.pot_params['k']
				self.dmax=2*self.sigma
				self.mult=1
			elif self.param.potential=='soft_attractive':
				self.k=self.param.pot_params['k']
				self.fact=self.param.pot_params['re_fact']-1.0
				self.rmax=1+2.0*self.fact
				self.dmax=2*self.sigma
				self.mult=self.rmax
			elif self.param.potential=='morse':
				self.D=self.param.pot_params['D']
				self.re=self.param.pot_params['re']
				self.a=self.param.pot.params['a']
				self.dmax=4*self.sigma
				self.mult=self.re
			elif self.param.potential=='gaussian':
				# give it at least a default neighbour radius ...
				self.dmax=2*self.sigma
				print "Warning! Gaussian interaction has not yet been implemented!" 
			elif self.param.potential=='rod':
				# CHANGE: Assuming here rods are at default aspect ratio 5
				self.dmax=10*self.sigma
				self.mult=1.0
				print "Warning! Rod interaction has not yet been implemented!"
			else:
				# give it at least a default neighbour radius ...
				self.dmax=2*self.sigma
				self.mult=1.0
				#print "Warning! Unknown interaction type! "
		else:
			# give it at least a default neighbour radius ...
			self.dmax=2*self.sigma
			self.mult=1.0
			print "Warning! Multiple types of particles interacting have not yet been implemented!"
		
	def getDmax(self):
		return self.dmax
	
	def getMult(self):
		return self.mult
		
	# also is -gradient potential	
	def getForce(self,i,neighbours,drvec,dr):
                Fvec=0.0*(np.array(neighbours).transpose()*(drvec).transpose()).transpose()
		if ((self.param.ntypes==1) or (self.ignore==True)):
			if self.param.potential=='soft':	
				diff=self.radius[i]+self.radius[neighbours]-dr
				fact = 0.5*self.k*diff
				Fvec=self.k*((diff/dr).transpose()*(drvec).transpose()).transpose()
				return Fvec
			elif self.param.potential=='soft_attractive':
				scale=self.radius[i]+self.radius[neighbours]
				diff=scale-dr
				dscaled=diff/scale
				rep = [index for index, value in enumerate(dscaled) if value > -self.fact]
				#print "Before upshot ones: " + str(len(rep))
				att = [index for index, value in enumerate(dscaled) if value <= -self.fact]
				#print "Attractive after upshot ones: " + str(len(att))
				factor=np.empty((len(neighbours),))
				# repulsive ones
				factor[rep] = self.k*diff[rep]
				# attractive ones
				factor[att]=-self.k*(self.rmax*scale[att]-dr[att])
				Fvec=((factor/dr).transpose()*(drvec).transpose()).transpose()
				return Fvec
			elif self.param.potential=='morse':
				fnorm=-2*self.a*self.D*np.exp(-self.a*(dr-self.re))*(1-np.exp(-self.a*(dr-self.re)))
				Fvec=((fnorm/dr).transpose()*(drvec).transpose()).transpose()
				return Fvec
			elif self.param.potential=='gaussian':
				print "Warning! Gaussian interaction has not yet been implemented! Returning zero force"
				return Fvec
			elif self.param.potential=='rod':
				print "Warning! Rod interaction has not yet been implemented! Returning zero force"
				return Fvec
			else:
				#print "Warning! Unknown interaction type! Returning zero force"
				return Fvec
		else:
			# Do the Morse right now only ... will serve as a template
			print "Warning! Multiple types of particles interacting have not yet been implemented! Returning zero force"
			return Fvec
				
	def getStresses(self,i,neighbours,drvec,dr):
		# Do these all the standard way, including the pressure as trace of the matrix
		stress=np.zeros((len(neighbours),3,3))
		# First get the forces:
		Fvec=self.getForce(i,neighbours,drvec,dr)
		for u in range(3):
			for v in range(3):
				stress[:,u,v]+=0.5*drvec[:,u]*Fvec[:,v]
		# Then get the pressure as the trace over the last two components
		press=np.trace(stress,axis1=1,axis2=2)
		return press,stress
	  
	def getEnergy(self,i,neighbours,drvec,dr):
		# Note: There is a 0.5 before the energy return statements because as written, every contact is counted twice in the calculation
		if ((self.param.ntypes==1) or (self.ignore==True)):
			if self.param.potential=='soft':	
				diff=self.radius[i]+self.radius[neighbours]-dr
				fact = 0.5*self.k*diff
				eng_val = fact*diff
				return 0.5*eng_val
			elif self.param.potential=='soft_attractive':
				scale=self.radius[i]+self.radius[neighbours]
				diff=scale-dr
				dscaled=diff/scale
				rep = [index for index, value in enumerate(dscaled) if value > -self.fact]
				#print "Before upshot ones: " + str(len(rep))
				att = [index for index, value in enumerate(dscaled) if value <= -self.fact]
				#print "Attractive after upshot ones: " + str(len(att))
				factor=np.empty((len(neighbours),))
				# repulsive ones
				factor[rep] = self.k*diff[rep]
				# attractive ones
				factor[att]=-self.k*(self.rmax*scale[att]-dr[att])
				eng_val=np.empty((len(neighbours),))
				eng_val[rep] = 0.5*factor[rep]*diff[rep]
				eng0=0.5*self.k*(self.fact*scale[att])**2
				eng_val[att] = eng0+(eng0-(factor[att]**2)/self.k)
				return 0.5*eng_val
			elif self.param.potential=='morse':
				eng_val=self.D*(1-np.exp(-self.a*(dr-self.re)))**2
				return 0.5*eng_val
			elif self.param.potential=='gaussian':
				print "Warning! Gaussian interaction has not yet been implemented! Returning zero energy"
				return 0
			elif self.param.potential=='rod':
				print "Warning! Rod interaction has not yet been implemented! Returning zero energy"
				return 0
			else:
				#print "Warning! Unknown interaction type! Returning zero energy"
				return 0
		else:
			# Do the Morse right now only ... will serve as a template
			print "Warning! Multiple types of particles interacting have not yet been implemented! Returning zero energy"
			return 0
	  
	def getStiffness(self,i,neighbours,drvec,dr):
		if ((self.param.ntypes==1) or (self.ignore==True)):
			if self.param.potential=='soft':
				return self.k*np.ones((len(neighbours),))
			elif self.param.potential=='soft_attractive':
				scale=self.radius[i]+self.radius[neighbours]
				diff=scale-dr
				dscaled=diff/scale
				rep = [index for index, value in enumerate(dscaled) if value > -self.fact]
				att = [index for index, value in enumerate(dscaled) if value <= -self.fact]
				stiff=np.zeros((len(neighbours),))
				stiff[rep]=self.k 
				stiff[att]=-self.k
				return stiff
			elif self.param.potential=='morse':
				stiff=2.0*self.a**2*self.D*np.exp(-self.a*(dr-self.re))*(2.0*np.exp(-self.a*(dr-self.re))-1)
				return stiff
			elif self.param.potential=='gaussian':
				print "Warning! Gaussian interaction has not yet been implemented! Returning zero stiffness"
				return np.zeros((len(neighbours),))
			elif self.param.potential=='rod':
				print "Warning! Rod interaction has not yet been implemented! Returning zero stiffness"
				return np.zeros((len(neighbours),))
			else:
				#print "Warning! Unknown interaction type! Returning zero stiffness"
				return np.zeros((len(neighbours),))
		else:
			# Do the Morse right now only ... will serve as a template
			print "Warning! Multiple types of particles interacting have not yet been implemented! Returning zero stiffness"