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
#  
# * ***************************************************************
#
#    Commented by Silke Henkes
# 
#    This is an example of a parsing tree (of level 2).                                                                                           	
#    The ReadConf class has a dictionary called key_words. These keywords refer to the 
#    first words in non-commented lines, of the type "pair_align polar { j = 0.1 }"                                                                                          
#    so here this would be the keyword 'pair_align'. Certain keywords, like for example
#    'dump', can have multiple lines associated to them. Ater creating an instance of the
#    class, e.g. conf=ReadConf(conffile), they can be found under conf.key_words,                                  
#    and accessed like conf.key_words['pair_align'] Each keyword has a list of parameters,
#    which refer to the rest of the lines, e.g. here 'polar', 'j' and 0.1. These are stored
#    in a name list, where 'polar goes, and an attribute list, where both the name of the 
#    variable 'j' and its value '0.1' go. For example, conf.key_words['pair_align'][0].name 
#    returns 'polar'. The attributes are at conf.key_words['pair_align'][0].attributes. 
#    Again, these might have multiple variables associated to them, if there are multiple
#    variables defined on each lines. In our example, there is only one, and we have                                                        
#    conf.key_words['pair_align'][0].attributes[0].name = ' j '                                                                                      
#    conf.key_words['pair_align'][0].attributes[0].val = '1'                                                                                         
#    Note that to get numerical values, we still need to do a conversion, e.g. float(conf.key_words['pair_align'][0].attributes[0].val = '1')        
#    (This appears to return double precision values on my desktop)                                                                                  
#
#    Reads in configuration files

import re

class Attribute:
  
  def __init__(self,attrib):
    ln = attrib.split('=')
    self.name = ln[0].strip()
    if len(ln) > 1:
      # strips the leading and trailing spaces and leaves me with the value (as a string still)
      self.val = ln[1].strip() 
    else:
      self.val = None

# 	This bit parses the 'rest' part of each line. If there are no {}, the parameter name is just the first element (word)
# 	Otherwise, we'll have at least one attribute. Look for anything within double braces, and split them along the ; signs. Each ' j = 1 ' bit then 
# 	Goes to make an attribute, with a name made from the first part and a value from the second
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
      # strip strips junk like tabs, empty spaces etc. map applies (lambda) function to every element of the list
      conf = map(lambda x: x.strip(), inp.readlines()) 
      # only lines longer than zero and don't start with '#'
      self.commands = filter(lambda x: len(x) > 0 and x[0] != '#', conf) 
      self.key_words = {}    #   keyword dictionary
      for comm_line in self.commands: 
        comm = comm_line.split()[0]            # split string at spaces into words; first part goes into comm
        rest = ' '.join(comm_line.split()[1:]) # all the rest goes here, joined up again
        if self.key_words.has_key(comm):
          self.key_words[comm].append(Parameter(rest))
        else:
          self.key_words[comm] = [Parameter(rest)]
