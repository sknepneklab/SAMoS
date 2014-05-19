# ################################################################
#  
#    Active Particles on Curved Spaces (APCS)
#    
#    Author: Silke Henkes
#
#    ICSMB, Department of Physics
#    University of Aberdeen 
#  
#    (c) 2014
#    
#    This program cannot be used, copied, or modified without
#    explicit permission of the authors.
#  
# ################################################################

import os
from read_conf import *

# Reads in parameter files, and places them in an intuitively readable class; using read_conf.py

class Param:
  def __init__(self,folder):
    for file in os.listdir(folder):
      if file.endswith(".conf"):
        print file
        self.filename=folder + file
        print self.filename
	
	# Reading in the different parameter values
	# See read_conf for an explanation
	# Yes, this is hamfisted ...
	# Configuration file key line examples:
	# box fixed { lx = 100;  ly = 100;  lz = 100 }
	# dump sphere_v0_0.02_j_0.1 { type=full; start=0; freq=10000; multi; coordinate; velocity; director; header  }
	# constraint sphere { r = 28.209479177387816 } (for a sphere)
	# constraint plane { lx = 100; ly = 100 } (for a plane)
	# pair_potential soft { k = 1 }
	# pair_align polar { j = 0.1 }
	# integrator nve { dt=0.001; limit=0.0001  }
	# integrator brownian { dt= 0.001; seed = 7;  nu = 0.002; mu = 1.0;  v0 = 0.02 }
	# run 10000 (for the relaxation step)
	# run 10000000 (note that's the second occurence; the one refering *NOT* to the relaxation step
	
	conf = ReadConf(self.filename)
	print conf.key_words['constraint'][0].name 
	self.simtype = conf.key_words['constraint'][0].name 
	if self.simtype == 'sphere':
	  self.r=float(conf.key_words['constraint'][0].attributes[0].val)
	  print 'Radius'
	  print self.r
	elif self.simtype == 'plane':
	  self.lx=float(conf.key_words['constraint'][0].attributes[0].val)
	  self.ly=float(conf.key_words['constraint'][0].attributes[1].val)
	  print 'Lx and Ly'
	  print self.lx
	  print self.ly
	else:
	  print 'Unknown simulation type (not sphere or plane)!'
	
	self.k=float(conf.key_words['pair_potential'][0].attributes[0].val)
	print 'Stiffness'
	print self.k
	self.j=float(conf.key_words['pair_align'][0].attributes[0].val)
	print 'Alignment strength'
	print self.j
	
	print "That's the first integrator:"
	print conf.key_words['integrator'][0].name
	print "That's the second integrator:"
	print conf.key_words['integrator'][1].name
	self.dt = float(conf.key_words['integrator'][1].attributes[0].val)
	self.seed = int(conf.key_words['integrator'][1].attributes[1].val)
	self.nu = float(conf.key_words['integrator'][1].attributes[2].val)
	self.mu = float(conf.key_words['integrator'][1].attributes[3].val)
	self.v0 = float(conf.key_words['integrator'][1].attributes[4].val)
	print 'Time step, seed, rotational noise, mobility, self-propulsion speed'
	print self.dt
	print self.seed
	print self.nu
	print self.mu
	print self.v0

	self.nsteps = int(conf.key_words['run'][1].name)
	print 'Simulation time steps'
	print self.nsteps
    


