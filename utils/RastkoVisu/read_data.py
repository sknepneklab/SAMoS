# ################################################################
#  
#    Active Particles on Curved Spaces (APCS)
#    
#    Author: Rastko Sknepnek
#   
#    Division of Physics
#    School of Engineering, Physics and Mathematics
#    University of Dundee
#    
#    (c) 2013
#    
#    This program cannot be used, copied, or modified without
#    explicit permission of the author.
#  
# ################################################################

# Reads in data files. Base handler for data analysis

class ReadData:
  
  def __init__(self, filename):
    with open(filename,'r') as self.inp:
      self.__read_data()
      self.inp.close()
    self.N = len(self.data[0])
      
  def __read_data(self):
    lines = self.inp.readlines()
    lines = map(lambda x: x.strip(), lines)
    if lines[0][0] == '#':
      self.has_header = True
    else:
      self.has_header = False
    if self.has_header:
      header = lines[0].strip()[1:]
      keys = header.split()
      self.keys = {}
      for k in keys: self.keys[k] = keys.index(k)
      lines = lines[1:]
      self.data = [[] for i in range(len(keys))]
    else:
      self.data = [[] for i in range(len(lines[0].split()))]
    for line in lines:
      data_line = map(float, line.split())
      for i in range(len(data_line)):
        self.data[i].append(data_line[i])

