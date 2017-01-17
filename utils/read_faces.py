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

import gzip

class ReadFaces:
  
  def __init__(self, filename):
	if filename.split('.')[-1] == 'gz':
	  self.lines = gzip.open(filename,'rb').read()
	else:
	  self.lines = open(filename,'r').read()
	self.__read_data()
	
      
  def __read_data(self):
	lines = self.lines.split('\n')
	lines = map(lambda x: x.strip(), lines)
	self.Nfaces = len(lines)
	# No keys here, just connections
	self.Faces = []
	i=0
	for line in lines:
	  data_line = map(int, line.split())
	  #print data_line[1:]
	  self.Faces.append(data_line[1:])
	

