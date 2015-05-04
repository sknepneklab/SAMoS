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

# Base class for Python particle

class Particle:
  
  def __init__(self, idx, tp=1, R=1.0, l=1.0):
    self.idx = idx
    self.tp = tp
    self.R = R
    self.r = [0.0,0.0,0.0]
    self.v = [0.0,0.0,0.0]
    self.f = [0.0,0.0,0.0]
    self.n = [0.0,0.0,0.0]
    self.l = l
    self.omega = 0.0
    
  