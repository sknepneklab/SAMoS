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

# Reads in data files. Base handler for data analysis

import gzip
from collections import OrderedDict

class ReadData:
  
  def __init__(self, filename):
	if filename.split('.')[-1] == 'gz':
	  self.lines = gzip.open(filename,'rb').read()
	else:
	  self.lines = open(filename,'r').read()
	self.__read_data()
	self.N = len(self.data[0])
      
  def __read_data(self):
	lines = self.lines.split('\n')
	lines = map(lambda x: x.strip(), lines)
	if lines[0][0] == '#':
	  self.has_header = True
	elif lines[0][0:5] == 'keys:':
	  self.has_header = True
	else:
	  self.has_header = False
	if self.has_header:
	  # Determine number of header lines and read the correct keys
	  nheader=0
	  read_keys=False
          self.keys = OrderedDict()
	  while read_keys==False:
		if lines[nheader][0] == '#':
		  header = lines[nheader].strip()[1:]
		  keys = header.split()
		  for k in keys: self.keys[k] = keys.index(k)
		  #print keys
		  nheader+=1
		  if (self.keys.has_key('id')) or (self.keys.has_key('type')) or (self.keys.has_key('radius')) or(self.keys.has_key('x')):
			read_keys=True
		elif lines[0][0:5] == 'keys:':
		  header = lines[nheader].strip()[6:]
		  keys = header.split()
		  for k in keys: self.keys[k] = keys.index(k)
		  nheader+=1
		  if (self.keys.has_key('id')) or (self.keys.has_key('type')) or (self.keys.has_key('radius')) or(self.keys.has_key('x')):
			read_keys=True
		else:
		  print "Error: There are no keys for these file, just comments! Cannot determine type of input!"
		  #print "Found " + str(nheader) + " lines of comments and keys "
		lines = lines[nheader:]
		self.data = [[] for i in range(len(keys))]
	else:
	  self.data = [[] for i in range(len(lines[0].split()))]
	for line in lines:
	  data_line = map(float, line.split())
	  for i in range(len(data_line)):
		self.data[i].append(data_line[i])

