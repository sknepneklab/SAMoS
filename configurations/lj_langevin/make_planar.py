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
from CellList2D import *



parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", type=str, default='particles.dat', help="output file name")
parser.add_argument("-L", "--length", type=float, default=50.0, help="simulation box size")
parser.add_argument("-N", "--num", type=int, default=100, help="number of particles")
parser.add_argument("-m", "--min_dist", type=float, default=2.0, help="minium distance between particles")
parser.add_argument("-v", "--init_vel", type=float, default=1.0, help="magnitude of initial velocity")
args = parser.parse_args()

print
print "\tSoft Actve Matter on Surfaces (SAMoS)"
print "\tGenerates particles in a square box"
print 
print "\tRastko Sknepnek"
print "\tUniversity of Dundee"
print "\t(c) 2016"
print "\t----------------------------------------------"
print
print "\tOutput files : ", args.output
print "\tBox size : ", args.length
print "\tNumber of particles : ", args.num
print "\tMinimum distance between paticles : ", args.min_dist
print "\tMagnitude of the initial velocity : ", args.init_vel
print "\tPacking fraction : ", (0.25*args.num*np.pi)/(args.length**2)  # we assume unit particle diamter 
print

start = datetime.now()

L = args.length
cl = CellList2D([L,L],2*args.min_dist)

particles = []
i = 0
while i < args.num:
  x, y = uniform(-0.5*L,0.5*L), uniform(-0.5*L,0.5*L)
  cid = cl.get_cell_idx((x,y))
  can_add = True
  for nb in cl.cell_list[cid].neighbors:
    for idx in cl.cell_list[nb].indices:
      xi, yi = particles[idx]
      dx, dy = x-xi, y-yi
      if dx*dx + dy*dy < args.min_dist**2: 
        can_add = False
        break
    if not can_add: 
      break
  if can_add:
    print "Successfully added particle : ", i
    particles.append((x,y))
    cl.add_particle((x,y),i)
    i += 1

out = open(args.output,'w')
# Careful! New output format / input format for code with keys
# Note: everything not specified is set to default values
# in this case nz=0 and z =0 
# Need to give initial values for the normal to the surface, nvx, nvy and nvy because constraint comes in later
# needs to be there to do tesselation and compute proper surfaces
# Last is native area ("polydispersity" can be introduced this way)
out.write('keys:  id  x y vx vy \n')
for i in range(len(particles)):
  x,y = particles[i]
  phi = uniform(0,2*m.pi)
  out.write('%4d  %f  %f  %f  %f\n' % (i,x,y, args.init_vel*m.cos(phi),args.init_vel*m.sin(phi)))
out.close()


end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print
