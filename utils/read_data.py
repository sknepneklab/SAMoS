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

# Reads in data files. Base handler for data analysis

from collections import OrderedDict
import gzip

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
	  while read_keys==False:
		if lines[nheader][0] == '#':
		  header = lines[nheader].strip()[1:]
		  keys = header.split()
		  self.keys = OrderedDict()
		  for k in keys: self.keys[k] = keys.index(k)
		  #print keys
		  nheader+=1
		  if (self.keys.has_key('id')) or (self.keys.has_key('type')) or (self.keys.has_key('radius')) or(self.keys.has_key('x')):
			read_keys=True
		elif lines[0][0:5] == 'keys:':
		  header = lines[nheader].strip()[6:]
		  keys = header.split()
		  self.keys = OrderedDict()
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

