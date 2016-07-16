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

# Utility code for generating initial configuration for simple particle simulation.
# This code places N particles ins a planar simulation box of size L

import sys
import argparse
import numpy as np
from random import uniform
from datetime import *
import math as m
from CellList import *



parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", type=str, default='particles.dat', help="output file name")
parser.add_argument("-d", "--density", type=float, default=0.5, help="number density")
parser.add_argument("-N", "--num", type=int, default=100, help="number of particles")
parser.add_argument("-m", "--min_dist", type=float, default=2.0, help="minium distance between particles")
parser.add_argument("-v", "--init_vel", type=float, default=1.0, help="magnitude of initial velocity")
args = parser.parse_args()

L = (args.num/args.density)**(1.0/3.0)

print
print "\tSoft Actve Matter on Surfaces (SAMoS)"
print "\tGenerates particles in a cubic box"
print 
print "\tRastko Sknepnek"
print "\tUniversity of Dundee"
print "\t(c) 2016"
print "\t----------------------------------------------"
print
print "\tOutput files : ", args.output
print "\tNumber density : ", args.density
print "\tNumber of particles : ", args.num
print "\tBox size : ", L
print "\tMinimum distance between paticles : ", args.min_dist
print "\tMagnitude of the initial velocity : ", args.init_vel
print



start = datetime.now()

cl = CellList([L,L,L],2*args.min_dist)

particles = []
i = 0
while i < args.num:
  x, y, z = uniform(-0.5*L,0.5*L), uniform(-0.5*L,0.5*L), uniform(-0.5*L,0.5*L)
  cid = cl.get_cell_idx((x,y,z))
  can_add = True
  for nb in cl.cell_list[cid].neighbors:
    for idx in cl.cell_list[nb].indices:
      xi, yi, zi = particles[idx]
      dx, dy, dz = x-xi, y-yi, z - zi
      if dx*dx + dy*dy + dz*dz < args.min_dist**2: 
        can_add = False
        break
    if not can_add: 
      break
  if can_add:
    print "Successfully added particle : ", i
    particles.append((x,y,z))
    cl.add_particle((x,y,z),i)
    i += 1

out = open(args.output,'w')
# Careful! New output format / input format for code with keys
# Note: everything not specified is set to default values
# in this case nz=0 and z =0 
# Need to give initial values for the normal to the surface, nvx, nvy and nvy because constraint comes in later
# needs to be there to do tesselation and compute proper surfaces
# Last is native area ("polydispersity" can be introduced this way)
out.write('keys:  id  x y z vx vy vz \n')
for i in range(len(particles)):
  x,y,z = particles[i]
  vx, vy, vz = uniform(-1.0,1.0), uniform(-1.0,1.0), uniform(-1.0,1.0)
  v = np.sqrt(vx*vx + vy*vy + vz*vz)
  vx *= args.init_vel/v
  vy *= args.init_vel/v
  vz *= args.init_vel/v
  out.write('%4d  %f  %f  %f  %f  %f  %f\n' % (i,x,y,z, vx, vy, vz))
out.close()


end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print
