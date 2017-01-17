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

# Utility code for generating intial configuration for cell simulations.
# This code places N cells in a patch of radius R keeing in mind that the 
# minimum distance between two cells shold be greater than a certain value.

import sys
import argparse
import numpy as np
from random import uniform
from datetime import *
import math as m
from CellList2D import *



parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", type=str, default='patch.dat', help="output file name")
parser.add_argument("-R", "--radius", type=float, default=20.0, help="patch radius")
parser.add_argument("-N", "--num", type=int, default=100, help="number of particles")
parser.add_argument("-m", "--min_dist", type=float, default=1.5, help="minium distance between particles")
parser.add_argument("-A", "--A0", type=float, default=m.pi, help="native cell area")
args = parser.parse_args()

print
print "\tSoft Actve Matter on Surfaces (SAMoS)"
print "\tGenerates a circial cell patch"
print 
print "\tRastko Sknepnek"
print "\tUniversity of Dundee"
print "\t(c) 2015"
print "\t----------------------------------------------"
print
print "\tOutput files : ", args.output
print "\tPatch radius : ", args.radius
print "\tNumber of cells : ", args.num
print "\tMinimum distance between cells : ", args.min_dist
print

start = datetime.now()

R = args.radius
cl = CellList2D([2.2*R,2.2*R],2*args.min_dist)

particles = []
i = 0
while i < args.num:
  x, y = uniform(-R,R), uniform(-R,R)
  if (x**2 + y**2 < R**2):
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
out.write('keys:  id  x y nx ny nvx nvy nvz area\n')
for i in range(len(particles)):
  x,y = particles[i]
  phi = uniform(0,2*m.pi)
  out.write('%4d  %f  %f  %f  %f  %f  %f  %f  %f\n' % (i,x,y, m.cos(phi),m.sin(phi), 0, 0, 1.0, args.A0))
out.close()


end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print