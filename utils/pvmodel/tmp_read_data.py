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

# Ok we still have no immediate plans to merge my python code with Silkes even
# though I haven't changed Silke code so merging it is trivial.
# I could copy this class and subclass it and use that but screw it. I'll just copy read_data.py and edit it
# to my liking.

# Most but not all of the changes to ReadData where to deal with environment particles 
#  at the same time as cells and boundary particles.
# If there no environment particles then the original ReadData object should just about work

import gzip
from collections import OrderedDict
import copy
import sys

class ReadData:

    def __init__(self, filename, isenv=False):
        if filename.split('.')[-1] == 'gz':
            self.lines = gzip.open(filename,'rb').read()
        else:
            self.lines = open(filename,'r').read()
        self.isenv = isenv  # flag for reading the environment
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

        # new object for the environment
        # Just looks the same as the data object though
        self.env = copy.deepcopy(self.data)
        # assume that the tissue precedes the environment

        # how to tell if the particle is in the tissue (?)
        if 'in_tissue' in self.keys:
            def in_tissue(dl): return dl[self.keys['in_tissue']] == 1 
        else:
            def in_tissue(dl): return True

        for line in lines:
            if not line: continue
            #should really cast data to the appropriate type, i.e. particle ids should be 'int'
            data_line = map(float, line.split())
            if in_tissue(data_line):
            #if True:
                # reading a tissue particle
                for i in range(len(data_line)):
                    self.data[i].append(data_line[i])
            elif not in_tissue(data_line) and self.isenv:
                # reading an environment particle
                for i in range(len(data_line)):
                    self.env[i].append(data_line[i])
