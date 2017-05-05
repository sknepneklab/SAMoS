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
# *   (c) 2014, 2015, 2016
# *  
# *   This program cannot be used, copied, or modified without
# *   explicit written permission of the authors.
# * 
# * ***************************************************************

# Reads in data files together with radius information and sticks them back together into a single file (for the equilibration)

import gzip
import sys
import argparse

from read_data import * 

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="input file ")
parser.add_argument("-r", "--radius", type=str, help="radius file")
parser.add_argument("-d", "--directory", type=str, help="input directory")
args = parser.parse_args()

# Read the two files into data structures
# Names
fileinput=args.directory + args.input
fileradius=args.directory + args.radius
# User ReadData to get them
data_input = ReadData(fileinput)
data_radius = ReadData(fileradius)

# get the radius out of the second, and update the first with the correct keyword
print data_radius.keys
radpos = data_radius.keys['radius']
print radpos
		  
if not data_input.keys.has_key('radius'): 
  print "will do the merger surgery now"
  print data_input.keys
  nkeys=len(data_input.data)
  print nkeys
  drad = {'radius':nkeys}
  print drad
  # I suppose it's easiest to just glue it to the end?
  data_input.keys.update(drad)
  print data_input.keys
  data_input.data.append(data_radius.data[radpos])
  #print data_input.data[nkeys]
  print data_input.keys.keys()
  nkeys+=1
else:
  print "Error: This file actually already has radii! Nothing to be done!"
  

keyval=[]
keylabel=[]
for u in range(nkeys):
	# this is admittedly backwards, but it has to be done: find the position of the key of that particular line number
	keyval.append(data_input.keys.values().index(u))
	keylabel.append(data_input.keys.keys()[keyval[u]])
print keyval
print keylabel

# First write the header labels
# Now we just need to write it to file again
f = open(args.directory + 'equilini.dat', 'w')
f.write('# ')
for u in range(nkeys):
	lpos=keyval.index(u)
	f.write(keylabel[u]+ '      ')
f.write('\n')
# and then the rest of the data ...
N=len(data_input.data[0])
print N
for i in range(N):
	for u in range(nkeys):
		f.write(str(data_input.data[u][i]) + '  ')
	f.write('\n')
f.close()

