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
	

