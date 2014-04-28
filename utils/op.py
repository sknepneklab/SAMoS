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

# Utility code for computing order parameter in planar or spherical geometry
# for a single frame

import sys
from math import *
import numpy as np


class OP:
  
  def __init__(self, data, inp_type):
    self.data = data
    if inp_type == "sphere":
      self.sphere = True
    else:
      self.sphere = False
    
  def compute(self):
    if self.data.has_header:
      vx, vy, vz = np.array(self.data.data[self.data.keys['vx']]), np.array(self.data.data[self.data.keys['vy']]), np.array(self.data.data[self.data.keys['vz']])
      v = np.vstack((vx,vy,vz)).transpose()
      if self.sphere:
        x, y, z = np.array(self.data.data[self.data.keys['x']]), np.array(self.data.data[self.data.keys['y']]), np.array(self.data.data[self.data.keys['z']])
        r = np.vstack((x,y,z)).transpose()
        len_r = np.apply_along_axis(np.linalg.norm,1,r).reshape((len(r),1))
        len_v = np.apply_along_axis(np.linalg.norm,1,v).reshape((len(v),1))
        r = r/len_r
        v = v/len_v
        loc_op = np.cross(r,v)
      else:
        loc_op = np.apply_along_axis(sum,0,v)/self.data.N
    return loc_op