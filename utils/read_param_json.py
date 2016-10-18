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
# *   (c) 2015, 2016
# * 
# *   Author: Silke Henkes
# * 
# *   Department of Physics 
# *   Institute for Complex Systems and Mathematical Biology
# *   University of Aberdeen  
# * 
# *   (c) 2014, 2015, 2016
# *  
# *   This program cannot be used, copied, or modified without
# *   explicit written permission of the authors.
# * 
# * ***************************************************************

import os
import json

# Deciphering json syntax:
# data['integrator']
# data['integrator']['brownian']
# float(data['integrator']['brownian']['v0'])

# This is so far incomplete - written to deal with a bunch of single-type cells with a boundary
# Expand later
# A bit messy since it needs to be compatible with the other parameter class
class Param:
	def __init__(self,filename):
		file=open(filename)
		self.data=json.load(dummy)
		print "Loaded json configuration file ", filename
		self.unpack_integrator()
		self.unpack_geometry()
		self.unpack_potentials()
		
		
        def unpack_integrator(self):
                # Pray that the order never changes ...
                # This is very much the case when there is only one integrator
                self.integrator=data['integrator'].keys()[0]
                self.dt=float(data['integrator']['dt'])
                print "Time step: " + str(self.dt)
                self.movegroup=data['integrator']['group']
                print "Moving group: " + self.movegroup
                if self.integrator=='brownian':
                        self.nu=float(data['integrator']['nu'])
                        print "Noise strength: " + str(self.nu)
                        self.mu=float(data['integrator']['mu'])
                        print "Mobility: " + str(self.mu)
                        self.seed=float(data['integrator']['seed'])
                        print "Dynamics seed: " + self.seed
			self.dt =float(self.int_params['dt'])
			try:
                                activemotion=data['potential']['external'].keys()
                                if activemotion[0]=='self_propulsion':
                                    self.v0=float(data['potential']['external']['self_propulsion']['alpha'])
                        except:
                                self.v0=0.0
			self.nematic=False
			# Un-mothball these once I actually use them and see how the code has been restructured
			#try:
				#dmp=self.int_params['nematic']
				#self.nematic=True
				#self.tau_flip=self.int_params['tau']
				#print "Nematic system with flip time " + str(self.tau_flip)
			#except:
				#pass
			#if self.potential=='rod':
				#self.mur = self.int_params['mur']
				#print "Rod rotational mobility " + str(self.mur)
			#self.thermal=False
			#try:
				#self.thermal_type=self.int_params['temperature_control']
				#self.thermal=True
				#self.kT=self.int_params['min_val'] # screw this, put in ramps only once we need them
				#self.kT_steps=self.int_params['steps']
				#print "Thermal brownian with " + self.thermal + " temperature " + self(kT) + " and steps " + str(self.kT_steps)
			#except:
				#pass
	
	# Again, here assume that we are dealing with a cell-type potential, which contains vertex, soft repulsion, line tension and boundary bending
	# This is only partially backwards compatible with the other parameter class: Before, we did not allow for more than one potential at a time
	# Each list element behaves like the previous version of potential data
        def unpack_potentials(self):
            self.potentials=data['potential']['pair'].keys()
            self.pot_params=[]
            for pot in self.potentials:
                # Add the generic bit of dictionary
                self.pot_params.append(data['potential']['pair'][pot])
                # And now define some values that are useful in general
                if pot=='vertex_particle':
                    self.kappa=float(data['potential']['pair'][pot]['K'])
                    self.gamma=float(data['potential']['pair'][pot]['gamma'])
                    self.lambdaval=float(data['potential']['pair'][pot]['lambda'])
                    
        # get the box size, the constraints, and whether we are periodic or not
        def unpack_geometry(self):
            consvals=data['constraint'].keys()
            try:
                consvals.remove('max_iter')
            except:
                pass
            try:
                consvals.remove('tol')
            except:
                pass
            try:
                self.constraintgroup=consvals('group')
                consvals.remove('group')
            except:
                pass
            self.constraint=consvals[0]
            print self.constraint
            # We *need* a box. Define it 100x100x100 unless specified otherwise
            self.box=[100.0,100.0,100.0]
            # We also need to know about periodicity. Set default to false.
            self.periodic='False'
            if self.constraint=='plane':
                self.box[0]=float(data['constraint']['plane']['lx'])
                self.box[1]=float(data['constraint']['plane']['ly'])
                self.box[2]=10.0
        
        # Stuff like the name of the dump files, the message files etc
        # later ... these are still part of the AnalyzeCellConfiguration for now
        #def unpack_auxiliary()
                    
                
                