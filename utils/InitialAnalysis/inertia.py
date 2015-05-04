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

# Utility code for computing moment of inertia of the system
# Note: Current implementation only works for non-periodic boundary conditions

import sys
from math import *
import numpy as np


class Inertia:
  
  def __init__(self, data):
    self.data = data
        
  def compute(self):
    if self.data.has_header:
      x, y, z = np.array(self.data.data[self.data.keys['x']]), np.array(self.data.data[self.data.keys['y']]), np.array(self.data.data[self.data.keys['z']])
      I = np.array([[sum(x**2),sum(x*y),sum(x*z)],[sum(y*x),sum(y**2),sum(y*z)],[sum(z*x),sum(z*y),sum(z**2)]])
    return np.linalg.eig(I)