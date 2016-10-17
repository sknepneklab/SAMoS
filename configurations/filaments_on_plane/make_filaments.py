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

# Utility code for building an input file for multiple filaments in plane

from datetime import *
from random import uniform, seed
from math import *
import argparse
import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument("-N", "--N", type=int, default=20, help="numper of particles in a polymer")
parser.add_argument("-M", "--M", type=int, default=50, help="numper of polymers")
parser.add_argument("-o", "--output", type=str, default='filaments.input', help="output file")
parser.add_argument("-i", "--input", type=str, default='filaments.xyz', help="output file")
args = parser.parse_args()

X = []
Y = []
with open(args.input,'r') as inp:
  lines = inp.readlines()
  for i in range(2,len(lines)):
    line = lines[i].strip()
    x, y = map(float,line.split()[1:3])
    X.append(x)
    Y.append(y)

bonds = []
for i in range(len(X)):
  if i == 0 or (i+1) % args.N != 0:
    bonds.append((i,i+1))

angles = []
for j in range(args.M):
  for i in range(1,args.N-1):
    angles.append((j*args.N+i-1,j*args.N+i,j*args.N+i+1))

out = open(args.output,'w')
out.write('keys: id  molecule type  x  y  nx  ny\n')
mol = 0
for i in range(len(X)):
  phi = uniform(0,2*pi)
  nx, ny = cos(phi), sin(phi)
  if i > 0 and i % args.N == 0: mol += 1
  out.write('%d  %d  %d %f  %f  %f  %f\n' % (i,mol,1,X[i],Y[i],nx,ny))
out.close()

out = open(args.output.split('.')[0]+'.bonds','w')
idx = 0
for (i,j) in bonds:
  out.write('%d  %d  %d  %d\n' % (idx,1,i,j))
  idx += 1
out.close()

out = open(args.output.split('.')[0]+'.angles','w')
idx = 0
for (i,j,k) in angles:
  out.write('%d  %d  %d  %d  %d\n' % (idx,1,i,j,k))
  idx += 1
out.close()


