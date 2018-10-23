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

import os
from read_conf import *

# Reads in parameter files, and places them in an intuitively readable class; using read_conf.py
# Reading in the different parameter values
# Careful, the code parser cares about CaPItalIzAtIoN
# See read_conf for an explanation
# Yes, this is hamfisted ...
# May 2015: Added infrastructure for the additional things which have been added to configuration files since the last edit 
# (and stuff which should have been in there in the first place)
# Configuration file key line examples (M for mandatory, O for optional):

# M messages
	# messages /where/the/message/file/is.msg

# M boxes
	# box fixed { lx = 100;  ly = 100;  lz = 100 }
	# box periodic { lx = 30.0;  ly = 30.0;  lz = 30.0 }

# M input
	# input /where/the/input/file/is.txt

# M neighbour list
	# nlist { rcut = 2.4; pad = 0.5 }

# O groups
	# group g1 { type = 1 }
	# group g2 { type = 2 }
	# group g3 { type = 3 }

# O outputs
	# Obscure formats
		# dump rod_test { type=velocity; start=0; freq=10000; multi; header }
		# dump rod_test { type=xyz; start=0; freq=10000; multi; header }
		# dump rod_test { type=dcd; start=0; freq=1000  }
		# dump rod_test_dir { type=xyzv; start=0; freq=1000; scale = 0.25; header  }
		# dump rod_test_vel { type=xyzv; start=0; freq=1000; scale = 0.1; velocity; header  }
		# dump rod_test { type=full; start=0; freq=1000; multi; coordinate; velocity; director; header  }
	# Common text format
		# dump sphere_v0_0.02_j_0.1 { type=full; start=0; freq=10000; multi; coordinate; velocity; director; header  }
		# dump plane_circle_test { type=full; start=0; freq=2000; multi; id; flag; tp; radius; coordinate; velocity; director; header  }
		# dump plane_test { type=full; start=0; freq=5000; multi; id; flag; tp; radius; coordinate; velocity; director; header  }
	# Contacts (underused)
		# dump contacts_v0_0.2_j_0.01 { type=contact; start=0; freq=10000; rcut = 2.4; multi  }
	
# M (technically O) constraints
	# sphere
		# constraint sphere { r = 28.209479177387816 } (for a sphere)
	# plane
		# constraint plane { lx = 100; ly = 100 } (for a plane)
		# constraint plane {  }
	# hourglass
		# We constrain all particles to move only on the surface of a hourglass if radius R = 10 with amplitude A = 2.0
		# It is assumed that we have only one waist (node)
		# constraint hourglass { R = 10.0; A = 2.0  }
	# cylinder
		# We constrain all particles to move only on the xy cylinder with lx=100, ly = 100 (actually those are read from the box)
		# constraint cylinder { r = 5 }
	# MISSING: catenoid, wavy cylinder
	
# M pair potentials (used for group=all by default)
	# soft, with and without radii
		# pair_potential soft { k = 1 }
		# pair_potential soft { k = 1.0; use_particle_radii}
	# gaussian
		# pair_potential gaussian { A = 2.0; B = -1.0 }
		# pair_potential gaussian { A = 3.0; B = -1.0; alpha = 3.0; beta = 3.0; rA = 0.0; rB = 1.0 }
	# morse
		# pair_potential morse { D = 1.5; a = 1.0; re = 1 }
		# pair_potential morse { D = 1.0; a = 3.0; re = 2.0; use_particle_radii }
	# rods
		# pair_potential rod { k = 1.0 }
	# MISSING: hertzian rods, hertzian?
	
# O pair alignment (used for group=all by default. Careful! J has to be capitalized ...)
	# pair_align polar { j = 0.1 }
	# pair_align vicsek { rcut = 1.0 }
	# pair_align nematic { J = 0.1 }
	
# O external alignment (active jamming, technically. Any external field fits as well)
	# external_align aj { tau = 0.01 }
	
# O logging (note that a log file with time step is generated no matter what)
	# log sphere_J_0.1_v0_0.1.log { freq = 1000; velocity; soft_energy; nematic_align; vec_velocity }
	# log sphere_test.log { freq = 1000; velocity; soft_energy; polar_align; vec_velocity }
	# log plane_test.log { freq = 1000; velocity; soft_energy; polar_align }
	# log plane_test.log { freq = 1000; velocity; rod_energy }
	# log hourglass_test.log { freq = 1000; velocity; soft_energy; polar_align; vec_velocity }
	# log plane_test.log { freq = 1000; velocity; soft_energy; polar_align }
	# log plane_test.log { freq = 1000; velocity; morse_energy }
	
# O NVE integrator
	# integrator nve { dt=0.001; limit=0.0001  }, together with
	# disable nve { group=all } (or even possibly disable nve)
	# integrator nve { dt=0.001; limit=0.0001; group = g1 }, together with
	# disable nve { group=g1 }
	
# O NVE integrator run time
	# run 10000 (for the relaxation step)
	
# O group-wise pair interaction parameters. Seems to have been used consistently after the NVE stage
	# pair_param morse { type_1 = 1; type_2 = 1; D = 1.5; a = 1.0; re = 1 }
	# pair_param morse { type_1 = 1; type_2 = 2; D = 3.0; a = 1.0; re = 1 }
	# pair_param morse { type_1 = 2; type_2 = 2; D = 1.5; a = 1.0; re = 1 }
	# pair_param soft { type_1 = 1; type_2 = 1; k=1.0 }
	# pair_param soft { type_1 = 1; type_2 = 2; k=10.0 }
	# pair_param soft { type_1 = 1; type_2 = 3; k=10.0 }
	# pair_param soft { type_1 = 2; type_2 = 2; k=1.0 }
	# pair_param soft { type_1 = 2; type_2 = 3; k=1.0 }
	# pair_param soft { type_1 = 3; type_2 = 3; k=1.0 }
	# pair_param gaussian { type_1 = 1; type_2 = 1; A = 3.0; B = 0.0; alpha = 3.0; beta = 3.0; rA = 0.0; rB = 1.0 }
	# pair_param gaussian { type_1 = 1; type_2 = 2; A = 3.0; B = -1.0; alpha = 3.0; beta = 3.0; rA = 0.0; rB = 1.0 }
	# pair_param gaussian { type_1 = 2; type_2 = 2; A = 3.0; B = 0.0; alpha = 3.0; beta = 3.0; rA = 0.0; rB = 1.0 }
	
# O group-wise pair alignment parameters. Seems to have been used consistently after the NVE stage
	# align_param polar { type_1 = 1; type_2 = 1; J = 0.25 }
	# align_param polar { type_1 = 1; type_2 = 2; J = 0.25 }
	# align_param polar { type_1 = 2; type_2 = 2; J = 0.0 }
	# align_param polar { type_1 = 1; type_2 = 3; J = 0.0 }
	# align_param polar { type_1 = 2; type_2 = 3; J = 0.0 }
	# align_param polar { type_1 = 3; type_2 = 3; J = 0.0 }

# M main integrator
	# Polar sphere brownian integrators
		# integrator brownian { dt= 0.001; seed = 7;  nu = 0.002; mu = 1.0;  v0 = 0.02 }
		# integrator brownian { dt=0.001; seed = 1;  nu = 0.002; mu = 1.0;  v0 = 0.05; group = g1 }
		# integrator brownian { dt=0.001; seed = 1;  nu = 0.00; mu = 1.0;  v0 = 1.0; group = all }
	# Nematic sphere brownian integrators, tau is flip time
		# integrator brownian { dt=0.001; seed = 22960;  nu = 0.0; mu = 1.0;  v0 = 0.1; nematic; tau = 1.0 }
	# Nematic rod brownian integrator
		# nu sets the width of the distribution for random changes of velocity
		# mu is particle mobility
		# mur is rotational rod mobility
		# v0 is the intensity of the self-propelling velocity
		# integrator brownian { dt=0.001; seed = 1;  nu = 0.00; mu = 1.0; mur = 1.0; v0 = 1.0; group = all; tau = 1.0; nematic }
	# Thermal brownian integrator (need more info!)
		# and another one for passive particles, but this time at temperature T = 0.3
		# temperature_control tells integrator to use constnat temeprature set by paramter
		# min_val. in this case max_val and steps are ignored. if we choose temperature_control=linear
		# then temeprature in linearly interpolated between min_val and max_val 
		# integrator brownian { group = passive; dt=0.001; seed = 4;  nu = 0.00; mu = 1.0;  v0 = 0.0; temperature_control=constant; min_val=0.3; max_val=0.3; steps = 1000 }
	# Vicsek integrator
		# integrator vicsek { dt=0.01; seed = 37;  eta = 1.0; mu = 1.0;  v0 = 0.5 }
		
# O Particle dynamics (division and death)
	# Stochastic divide and death
		# population random { group = g1; division_rate = 1000.0; death_rate = 1000.0; freq = 1000 }
	# Density-controlled divide and death; rho_max is in contact number units
		# Simple everyone
			# population density { group = all; division_rate = 0.0003; split_distance=0.1,rho_max = 6.0,death_rate = 0.00006; freq = 1000  }
		# While switching groups
			# population density { group = g1; division_rate = 0.0; death_rate = 0.00025; freq = 1000; change_prob_1 = 0.0; change_prob_2 = 0.0 , old_group = g1; new_group = g1; old_type = 1; new_type = 1 }
			# population density { group = g2; division_rate = 0.1; death_rate = 0.0; freq = 1000; change_prob_1 = 0.0; split_distance = 0.05; change_prob_2 = 1.0; old_group = g2; new_group = g1; old_type = 2; new_type = 1  }
			# population density { group = g2; division_rate = 0.025; death_rate = 0.0; freq = 1000; split_distance = 0.0;rho_max = 1000.0, change_prob_2 = 1.0; old_group = g2; new_group = g1; old_type = 2; new_type = 1 } 

# M main running time
	# run 10000000 (note that's the second occurence; the one refering *NOT* to the relaxation step

class Param:
	def __init__(self,filename):
		#for file in os.listdir(folder):
			#if file.endswith(".conf"):
				#self.filename=folder + file
				#print self.filename
				
		self.filename=filename
		print filename
		conf = ReadConf(self.filename)
		
		# Message file
		self.messagefile = conf.key_words['messages'][0].name
		print "Message file: " + self.messagefile
		# Boxes
		self.boxtype = conf.key_words['box'][0].name 
		print "Box type: " + self.boxtype
		self.box=[]
		self.box.append(float(conf.key_words['box'][0].attributes[0].val))
		self.box.append(float(conf.key_words['box'][0].attributes[1].val))
		self.box.append(float(conf.key_words['box'][0].attributes[2].val))
		print "Box dimensions: " 
		print self.box
		# Neighbour list
		self.nlist_rcut = float(conf.key_words['nlist'][0].attributes[0].val)
		self.nlist_pad = float(conf.key_words['nlist'][0].attributes[1].val)
		print "Neighbour list rcut " + str(self.nlist_rcut) + " and padding " + str(self.nlist_pad)
		# Input file
		self.inputfile = conf.key_words['input'][0].name
		print "Input file: " + self.inputfile
		# Dump parameters
		self.dumpname=conf.key_words['dump'][0].name
		self.dump={}
		for l in range(len(conf.key_words['dump'][0].attributes)):
			try:
				self.dump[str.strip(conf.key_words['dump'][0].attributes[l].name)]=float(conf.key_words['dump'][0].attributes[l].val)
			except:
				try:
					self.dump[str.strip(conf.key_words['dump'][0].attributes[l].name)]=str.strip(conf.key_words['dump'][0].attributes[l].val)
				except: # no constraints
					pass
		
		# Groups
		try:
			self.ngroups=len(conf.key_words['group'])
			self.groupnames=[]
			self.grouptypes=[]
			for l in range(self.ngroups):
				self.groupnames.append(conf.key_words['group'][l].name) 
				self.grouptypes.append(int(conf.key_words['group'][l].attributes[0].val))
			print "Group names: " 
			print self.groupnames
			print "Group types: " 
			print self.grouptypes
		except KeyError: 
			self.ngroups=1
		print "Number of groups: " + str(self.ngroups)
				
		# Constraints
		try:
			self.constraint = conf.key_words['constraint'][0].name 
			print "Constraint: " + self.constraint
			self.const_params={}
			try:
				for l in range(len(conf.key_words['constraint'][0].attributes)):
					try:
						self.const_params[str.strip(conf.key_words['constraint'][0].attributes[l].name)]=float(conf.key_words['constraint'][0].attributes[l].val)
					except:
						try:
							self.const_params[str.strip(conf.key_words['constraint'][0].attributes[l].name)]=str.strip(conf.key_words['constraint'][0].attributes[l].val)
						except: # no constraints
							pass
				print "Constraint parameters " 
				print self.const_params
			except KeyError:
				pass
			# Legacy for the two most common constraints 
			if self.constraint == 'sphere':
				self.r=float(conf.key_words['constraint'][0].attributes[0].val)
				print 'Radius'
				print self.r
			elif self.constraint == 'plane':
				try: # try from the constraint
					self.lx=float(conf.key_words['constraint'][0].attributes[0].val)
					self.ly=float(conf.key_words['constraint'][0].attributes[1].val)
				except:# else use the box
					self.lx=self.box[0]
					self.ly=self.box[1]
				print 'Lx and Ly'
				print self.lx
				print self.ly  
				if self.boxtype=='periodic':
					self.constraint='plane_periodic'
		except KeyError:
			self.constraint='none'
			self.const_params={}
							
		# Default pair potentials and aligners
		try:
			self.potential=conf.key_words['pair_potential'][0].name 
			self.pot_params={}
			for l in range(len(conf.key_words['pair_potential'][0].attributes)):
				try:
					self.pot_params[str.strip(conf.key_words['pair_potential'][0].attributes[l].name)]=float(conf.key_words['pair_potential'][0].attributes[l].val)
				except:
					try:
						self.pot_params[str.strip(conf.key_words['pair_potential'][0].attributes[l].name)]=str.strip(conf.key_words['pair_potential'][0].attributes[l].val)
					except: # use_particle_radii: different form ...
						self.pot_params[str.strip(conf.key_words['pair_potential'][0].attributes[l].name)]=True
		except KeyError:
			self.potential = 'none'
			self.pot_params={}
		try:
			self.aligner=conf.key_words['pair_align'][0].name 
			self.J=float(conf.key_words['pair_align'][0].attributes[0].val)
		except KeyError: 
			try:
				self.aligner=conf.key_words['external_align'][0].name 
				self.J=float(conf.key_words['external_align'][0].attributes[0].val)
			except KeyError: 
				self.aligner='none'
				self.J=0.0
		print "Potential: " + self.potential
		print "Parameters: " 
		print self.pot_params
		print "Aligner: " +self.aligner
		print "J: " + str(self.J)
		
		# Something for our friends the cells
		# As written, we only ever read the first potential ... generalise
		#pair_potential vp { K = 1.0; gamma = 1.0; lambda = -6.283184000}
		#pair_potential line_tension { lambda = 0.0 }
		#pair_potential soft {k = 10.0; a=0.5}
		if self.potential=='vp':
			self.kappa=float(conf.key_words['pair_potential'][0].attributes[0].val)
			self.gamma=float(conf.key_words['pair_potential'][0].attributes[1].val)
			self.lambdaval=float(conf.key_words['pair_potential'][0].attributes[2].val)
			
		# NVE integrator
		# Everything is based on the assumption that there is only one of these, currently ...
		nNVE=0
		if conf.key_words['integrator'][0].name=='nve':
			print "NVE integrator " 
			self.nstepsNVE= int(conf.key_words['run'][0].name)
			nNVE+=1
		else:
			self.nstepsNVE=0
		print "NVE steps: " + str(self.nstepsNVE)
		print "NVE integrators: " + str(nNVE)
				
		# Type-wise pair potentials and aligners (careful: types and groups don't have to match!)
		# square lists of lists of dictionaries or names
		# Default: initialize with the defaults (including the 'none' if applicable)
		# Minor leap of faith: types are numbered 0 1 2 etc., and all of them are part of *some* group
		if self.ngroups>1:
			self.ntypes=max(self.grouptypes)
			print "Number of types: " + str(self.ntypes)
			if self.ntypes>1:
				# Potentials
				self.type_potential=[[self.potential for u in range(self.ntypes)] for u in range(self.ntypes)]
				self.type_pot_params=[[self.pot_params for u in range(self.ntypes)] for u in range(self.ntypes)]
				try:
					for l in range(len(conf.key_words['pair_param'])):
						type1=int(conf.key_words['pair_param'][l].attributes[0].val)-1
						type2=int(conf.key_words['pair_param'][l].attributes[1].val)-1
						potential=conf.key_words['pair_param'][l].name
						self.type_potential[type1][type2]=potential
						self.type_pot_params[type1][type2]={}
						for m in range(2,len(conf.key_words['pair_param'][l].attributes)):
							try:
								self.type_pot_params[type1][type2][str.strip(conf.key_words['pair_param'][l].attributes[m].name)]=float(conf.key_words['pair_param'][l].attributes[m].val)
							except:
								try:
									self.type_pot_params[type1][type2][str.strip(conf.key_words['pair_param'][l].attributes[m].name)]=float(conf.key_words['pair_param'][l].attributes[m].val)
								except: # use_particle_radii
									self.type_pot_params[type1][type2][str.strip(conf.key_words['pair_param'][l].attributes[m].name)]=True
				except KeyError:
					pass
				# Aligners
				self.type_aligner=[[self.aligner for u in range(self.ntypes)] for u in range(self.ntypes)]
				self.type_J=[[self.J for u in range(self.ntypes)] for u in range(self.ntypes)]
				try:
					for l in range(len(conf.key_words['align_param'])):
						type1=int(conf.key_words['align_param'][l].attributes[0].val)-1
						type2=int(conf.key_words['align_param'][l].attributes[1].val)-1
						aligner=conf.key_words['align_param'][l].name
						self.type_aligner[type1][type2]=aligner
						self.type_J[type1][type2]=float(conf.key_words['align_param'][l].attributes[2].val)
				except:
					pass
				print "Type potentials: " 
				print self.type_potential
				print "Type potential parameters" 
				print self.type_pot_params
				print "Type aligners: " 
				print self.type_aligner
				print "Type J" 
				print self.type_J
		else:
			self.ntypes=1
				
		# Main integrator(s)
		# Define straightforward paramters for the most common Brownian one
		# First: distinguish between groups and no groups
		self.one_integrator=False
		if self.ngroups==1:
			self.integrator=conf.key_words['integrator'][nNVE].name
			print "Main integrator: " + self.integrator
			self.int_params={}
			for l in range(len(conf.key_words['integrator'][nNVE].attributes)):
				try:
					self.int_params[str.strip(conf.key_words['integrator'][nNVE].attributes[l].name)]=str.strip(conf.key_words['integrator'][nNVE].attributes[l].val)
				except: # some odd thermal ones are effectively boolean
					self.int_params[str.strip(conf.key_words['integrator'][nNVE].attributes[l].name)]=True
			done = self.oneInt(conf)
		else:
			self.group_integrator=['none' for u in range(self.ngroups)] 
			self.group_int_params=[{} for u in range(self.ngroups)]
			nintegrator=len(conf.key_words['integrator'])
			print "Found " + str(nintegrator) + " intergrators!"
			for k in range(nNVE,nintegrator): # Excluding the NVE here
				int_params={}
				for l in range(len(conf.key_words['integrator'][k].attributes)):
					try:
						int_params[str.strip(conf.key_words['integrator'][k].attributes[l].name)]=str.strip(conf.key_words['integrator'][k].attributes[l].val)
					except:
						int_params[str.strip(conf.key_words['integrator'][k].attributes[l].name)]=True
				# now sort them into groups
				# first in case it's all of them
				# I don't care if some idiot has added more integrators on top of it. That's their problem.
				try:
					mygroup = int_params['group']
				except KeyError:
					mygroup='all'
				if mygroup =='all':
					self.integrator=conf.key_words['integrator'][k].name
					print "Main integrator: " + self.integrator
					self.int_params=int_params
					done = self.oneInt(conf)
				else:
					groupidx=self.groupnames.index(mygroup)
					self.group_integrator[groupidx]=conf.key_words['integrator'][k].name
					self.group_int_params[groupidx]=int_params
					if (nintegrator-nNVE)==1: #only one moving group, for example
						self.integrator=conf.key_words['integrator'][k].name
						print "Main integrator: " + self.integrator
						self.int_params=int_params
						done = self.oneInt(conf)
					else:
						if k==(nintegrator-1):
							print "Warning: multiple complex integrators "
							print self.group_integrator 
							print " for groups " 
							print self.groupnames
							print "Parameters are stored in the dictionary self.group_int_params:"
							print self.group_int_params
		
		# Population control
		# MISSING: The fade-in options 
		# Since this is very underdeveloped, just store the name(s), and the various options in a dictionary
		# Warning: the values of these parameters are strings, even the ones that should be int or double
		try:
			self.npopulation = len(conf.key_words['population'])
			if self.npopulation>0:
				print "Number of populations: " + str(self.npopulation)
				self.population=[]
				self.pop_params=[{} for k in range(self.npopulation)]
				for k in range(self.npopulation):
					self.population.append(conf.key_words['population'][k].name)
					for l in range(len(conf.key_words['population'][k].attributes)):
						try:
							self.pop_params[k][str.strip(conf.key_words['population'][k].attributes[l].name)]=float(conf.key_words['population'][k].attributes[l].val)
						except:
							self.pop_params[k][str.strip(conf.key_words['population'][k].attributes[l].name)]=str.strip(conf.key_words['population'][k].attributes[l].val)
				print "Populations: "
				print self.population
				print "Population parameters: "
				print self.pop_params
		except KeyError:
			self.npopulation=0
			pass
		
		self.nsteps = int(conf.key_words['run'][nNVE].name)
		print 'Simulation time steps'
		print self.nsteps
	
	
	def oneInt(self,conf):
		self.one_integrator=True
		if self.integrator=='brownian':
                        # In case it's one of the newer ones where the dt is on its own
                        # If the except doesn't work either, we are fucked in any case
                        try:
                            self.dt =float(self.int_params['dt'])
                        except:
                            #print conf.key_words['timestep']
                            #self.dt = float(conf.key_words['timestep'])
                            # Oh yes, that syntax is incompatible with the parser
                            # And leads to all kind of BS if I don't stay on top of it
                            self.dt=0.01
                        print "Time step: " + str(self.dt)
			self.seed = self.int_params['seed']
			print "Dynamics seed: " + self.seed
			self.mu = self.int_params['mu']
			print "Mobility: " + str(self.mu)
			# Again, the stupid v0 as external aligner type
			try:
                            self.v0 = float(conf.key_words['external'][0].attributes[0].val)
                        except:
                            self.v0 = self.int_params['v0']
			print "v0: " + str(self.v0)
			self.nu = self.int_params['nu']
			print "Noise strength: " + str(self.nu)
			self.nematic=False
			try:
				dmp=self.int_params['nematic']
				self.nematic=True
				self.tau_flip=self.int_params['tau']
				print "Nematic system with flip time " + str(self.tau_flip)
			except:
				pass
			if self.potential=='rod':
                                try:
                                    self.mur = self.int_params['mur']
                                    print "Rod rotational mobility " + str(self.mur)
                                except:
                                    self.mu = self.int_params['mu']
                                    print "Rod rotational mobility " + str(self.mu)
				
			self.thermal=False
			try:
				self.thermal_type=self.int_params['temperature_control']
				self.thermal=True
				self.kT=self.int_params['min_val'] # screw this, put in ramps only once we need them
				self.kT_steps=self.int_params['steps']
				print "Thermal brownian with " + self.thermal + " temperature " + self(kT) + " and steps " + str(self.kT_steps)
			except:
				pass
			self.movegroup='all'
			try:
				self.movegroup = self.int_params['group']
			except:
				pass
			print "Moving group: " + self.movegroup
		elif self.integrator=='vicsek':
			self.dt =float(self.int_params['dt']) 
			self.seed = self.int_params['seed']
			self.mu = self.int_params['mu']
			self.v0 = self.int_params['v0']
			self.nu = self.int_params['eta']
		else:
			self.dt =float(self.int_params['dt']) 
			print "Warning: unknown integrator type " + self.integrator + ". Parameters are stored in the dictionary self.int_params."
			return 1
		return 0
	