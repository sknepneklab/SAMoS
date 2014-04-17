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
#    (c) 2013, 2014
#    
#    This program cannot be used, copied, or modified without
#    explicit permission of the author.
#  
# ################################################################

# Reads in configuration files

import re

class Attribute:
  
  def __init__(self,attrib):
    ln = attrib.split('=')
    self.name = ln[0]
    if len(ln) > 1:
      self.val = ln[1]
    else:
      self.val = None

class Parameter:
  
  def __init__(self,line):
    ln = line.split()
    if ln[0] != '{':
      self.name = ln[0]
    else:
      self.name = None
    self.attributes = []
    if len(line.split()) > 1:
      res = re.search('{(.*)}',line).group(1)
      attribs = res.split(';')
      for attrib in attribs:
        self.attributes.append(Attribute(attrib))

class ReadConf:
  
  def __init__(self, file_name):
    with open(file_name) as inp:
      conf = map(lambda x: x.strip(), inp.readlines())
      self.commands = filter(lambda x: len(x) > 0 and x[0] != '#', conf)
      self.key_words = {}
      for comm_line in self.commands:
        comm = comm_line.split()[0]
        rest = ' '.join(comm_line.split()[1:])
        if self.key_words.has_key(comm):
          self.key_words[comm].append(Parameter(rest))
        else:
          self.key_words[comm] = [Parameter(rest)]