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

# Utility code for computing potential energy on each particle and
# colouring it accordingly. Output will be in a POVRay file

from read_data import *

from datetime import *
from random import uniform 
from math import *
import numpy as np
from vtktools import *
import argparse

class StressEnergy:
  
  def __init__(self,data,box=None):
    if not (data.keys.has_key('x') and data.keys.has_key('y') and data.keys.has_key('z')):
      raise "Particle coordinate not specified in the input data."
    self.data = data
    self.box = box
    self.X, self.Y, self.Z = self.data.data[self.data.keys['x']], self.data.data[self.data.keys['y']], self.data.data[self.data.keys['z']]
    if not self.data.keys.has_key('r'): self.A = [1.0 for i in range(self.data.N)]
    else: self.A = self.data.data[self.data.keys['r']]
    if self.box != None:
      self.lx, self.ly, self.lz = box
    
  def flip_boundary(self, dx, dy, dz):
    if self.box != None:
      if (dx > 0.5*self.lx): dx -= self.lx
      elif (dx < -0.5*self.lx): dx += self.lx
      if (dy > 0.5*self.ly): dy -= self.ly
      elif (dy < -0.5*self.ly): dy += self.ly
      if (dz > 0.5*self.lz): dz -= self.lz
      elif (dz < -0.5*self.lz): dz += self.lz
      
  def compute_harmonic(self, k = 1.0):
    eng = [0.0 for i in range(self.data.N)]
    idx = 0
    for (x0,y0,z0,a0) in zip(self.X,self.Y,self.Z,self.A):
      for (x1,y1,z1,a1) in zip(self.X,self.Y,self.Z,self.A):
        dx, dy, dz = x1-x0, y1-y0, z1-z0
        self.flip_boundary(dx,dy,dz)
        dr = sqrt(dx*dx+dy*dy+dz*dz)
        if (dr > 0.0):
          if dr < a0+a1:
            diff = a0+a1-dr
            eng[idx] += 0.5*k*diff*diff
      idx += 1
    return eng
  
  def compute_pressure(self, k = 1.0):
    press = [0.0 for i in range(self.data.N)]
    idx = 0
    for (x0,y0,z0,a0) in zip(self.X,self.Y,self.Z,self.A):
      for (x1,y1,z1,a1) in zip(self.X,self.Y,self.Z,self.A):
		dx, dy, dz = x1-x0, y1-y0, z1-z0
		self.flip_boundary(dx,dy,dz)
		dr = sqrt(dx*dx+dy*dy+dz*dz)
		if (dr > 0.0):
		  if dr < a0+a1:
			diff = a0+a1-dr
			press[idx] += 0.5*k*diff*dr
      idx += 1
    return press

  def compute_sigmaij(self, k = 1.0):
    sigmaij = [[[0.0 for i in range(3)] for i in range(3)] for i in range(self.data.N)]
    idx = 0
    for (x0,y0,z0,a0) in zip(self.X,self.Y,self.Z,self.A):
      for (x1,y1,z1,a1) in zip(self.X,self.Y,self.Z,self.A):
		dx, dy, dz = x1-x0, y1-y0, z1-z0
		self.flip_boundary(dx,dy,dz)
		dr = sqrt(dx*dx+dy*dy+dz*dz)
		if (dr > 0.0):
		  if dr < a0+a1:
			diff = a0+a1-dr
			rvec=[dx, dy, dz]
			for i in range(3):
			  for j in range(3):
				sigmaij[i][j][idx] += 0.5*k*diff*rvec[i]*rvec[j]/dr
      idx += 1
    return sigmaij

class POVPrint:
 
  low_colour = [0,0,1.0]
  hi_colour = [1,0,0]
  sphere_radius = 10.0
 
  def __init__(self,outfilename, data, eng):
    self.outfilename = outfilename
    self.data = data
    self.eng = eng
  
  def write(self):
    self.out = open(self.outfilename,'w')
    self.out.write('global_settings{ assumed_gamma 1.0 }\n')
    self.out.write('#default{ finish{ ambient 0.1 diffuse 0.9 }}\n')
    self.out.write('#include "colors.inc"\n')
    self.out.write('#include "textures.inc"\n')
    self.out.write('#include "glass.inc"\n')
    self.out.write('#include "metals.inc"\n')
    self.out.write('#include "golds.inc"\n')
    self.out.write('#include "stones.inc"\n')
    self.out.write('#include "woods.inc"\n')
    self.out.write('#declare Camera_0 = camera { perspective\n')
    self.out.write('                             angle 11\n')
    self.out.write('                             right     x*image_width/image_height\n')
    self.out.write('                             location  < 0.00, 0.00,%f>\n' % (-16.66666*self.sphere_radius))
    self.out.write('                             look_at   < 0.00, 0.00, 0.00>\n')
    self.out.write('                            }\n')
    self.out.write('camera{Camera_0}\n')
    self.out.write('light_source{<1500, 500,-2500> color White}\n')
    self.out.write('sky_sphere{ pigment{ gradient <0,1,0>\n')
    self.out.write('                 color_map{ [0   color rgb<1,1,1>         ]//White\n')
    self.out.write('                            [0.4 color rgb<0.14,0.14,0.56>]//~Navy\n')
    self.out.write('                            [0.6 color rgb<0.14,0.14,0.56>]//~Navy\n')
    self.out.write('                            [1.0 color rgb<1,1,1>         ]//White\n')
    self.out.write('                          }\n')
    self.out.write('                 scale 2 }\n')
    self.out.write('       } // end of sky_sphere\n')
    self.__write_particles()
    self.out.close()

  def __write_particles(self):
    min_eng, max_eng = min(self.eng), max(self.eng)
    eng_range = max_eng - min_eng
    X, Y, Z = self.data.data[self.data.keys['x']], self.data.data[self.data.keys['y']], self.data.data[self.data.keys['z']]
    if not self.data.keys.has_key('r'): R = [1.0 for i in range(self.data.N)]
    else: R = self.data.data[self.data.keys['r']]
    for (x,y,z,r,e) in zip(X,Y,Z,R,self.eng):
      tt = (e-min_eng)/eng_range
      c_r = (1.0-tt)*self.low_colour[0] + tt*self.hi_colour[0]
      c_g = (1.0-tt)*self.low_colour[1] + tt*self.hi_colour[1]
      c_b = (1.0-tt)*self.low_colour[2] + tt*self.hi_colour[2]
      self.out.write('sphere { <%f,%f,%f>,%f\n' % (x,y,z,r)) 
      self.out.write('          texture\n')
      self.out.write('          {\n')
      self.out.write('            pigment\n')
      self.out.write('            {\n')         
      self.out.write('              color rgb<%f,%f,%f>\n' % (c_r,c_g,c_b))
      self.out.write('            }\n')
      self.out.write('            finish\n')
      self.out.write('            { \n')
      self.out.write('              reflection 0.05\n')
      self.out.write('              specular 0.05 \n')
      self.out.write('              ambient 0.5 \n')
      self.out.write('              diffuse 0.5 \n')
      self.out.write('            }\n')                                          
      self.out.write('          }\n')
      self.out.write('       }\n')

class XYZCPrint:
 
  def __init__(self,outfilename, data, eng):
    self.outfilename = outfilename
    self.data = data
    self.eng = eng
  
  def write(self):
    self.out = open(self.outfilename,'w')
    self.out.write('%d\n' % self.data.N)
    self.out.write('#Generate by pot_energy_plot.py\n')
    self.__write_particles()
    self.out.close()

  def __write_particles(self):
    min_eng, max_eng = min(self.eng), max(self.eng)
    eng_range = max_eng - min_eng
    X, Y, Z = self.data.data[self.data.keys['x']], self.data.data[self.data.keys['y']], self.data.data[self.data.keys['z']]
    for (x,y,z,e) in zip(X,Y,Z,self.eng):
      self.out.write('A %f  %f  %f  %f\n' % (x,y,z,e))

# Scripting version: Only execute if this is called as a script. Otherwise, it attempts to go through here when loading as a module 
# and throws errors because some arguments aren't defined
if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", type=str, help="Input file with particle velocity field")
	parser.add_argument("-o", "--output", type=str, default="StressEnergy", help="Output file (POV-Ray scene script)")
	parser.add_argument("-k", "--k", type=float, default=1.0, help="soft potential strength")
	parser.add_argument("-R", "--sphere_r", type=float, default=10.0, help="radius of sphere for spherical system")
	parser.add_argument("-L", "--low_colour", type=float, nargs=3, default=[0.0,0.0,1.0], help="lowest energy colour")
	parser.add_argument("-H", "--hi_colour", type=float, nargs=3, default=[1.0,0.0,0.0], help="highest energy colour")
	parser.add_argument("-c", "--connectivity", type=str, default=None, help="Connectivity file (xyzl format produced by stripack)")
	args = parser.parse_args()

	print
	print "\tActive Particles on Curved Spaces (APCS)"
	print "\tPotential energy distribution"
	print 
	print "\tRastko Sknepnek"
	print "\tUniversity of Dundee"
	print "\t(c) 2013"
	print "\t----------------------------------------------"
	print 
	print "\tInput : ", args.input
	print "\tOutput : ", args.output
	print "\tLow colour : ", args.low_colour
	print "\tHi colour : ", args.hi_colour
	if args.connectivity != None:
	  print "\tConnectivity file : ", args.connectivity
	print 

	start = datetime.now()

	print "Reading data..."
	data = ReadData(args.input)
	str_eng = StressEnergy(data)

	print "Computing potential energy and stresses ..."
	engs = str_eng.compute_harmonic(args.k)
	press = str_eng.compute_pressure(args.k)
	sigmaij = str_eng.compute_pressure(args.k)

	# later ...
	#print "Generating POV Ray output..."
	#p = POVPrint(args.output+'.pov',data,eng)
	#p.sphere_radius = args.sphere_r
	#p.hi_colour = args.hi_colour
	#p.low_colour = args.low_colour
	#p.write()

	#print "Generating XYZC output..."
	#xyzc = XYZCPrint(args.output+'.xyzc',data,engs)
	#xyzc.write()

	print "Generating VTU output..."
	x, y, z = np.array(data.data[data.keys['x']]), np.array(data.data[data.keys['y']]), np.array(data.data[data.keys['z']])
	vtk_writer = VTK_XML_Serial_Unstructured()

	nneigh = []
	if args.connectivity != None:
	  nneigh = [0 for i in range(data.N)]
	  print "Generating connectivities..."
	  con = open(args.connectivity,'r')
	  edges = [f for f in map(lambda x: int(x.strip()), con.readlines()) if f != -1]
	  edge_pairs = zip(edges[::2],edges[1::2])
	  for (i,j) in edge_pairs:
		nneigh[i-1] += 0.5
		nneigh[j-1] += 0.5
	  
	if len(nneigh) == 0:
	  vtk_writer.snapshot(args.output+'.vtu',x,y,z,energies=engs,pressures=press)
	else:
	  vtk_writer.snapshot(args.output+'.vtu',x,y,z,energies=engs,pressures=press,nneigh=nneigh)

	print "Generating XYZ file for SRTIPACK triangulation..."
	out = open(args.output+'.xyz','w')
	for (xx,yy,zz) in zip(x,y,z):
	  out.write('%f  %f  %f\n' % (xx,yy,zz))
	out.close()

	print "Lowest energy : ", min(engs)
	print "Highest energy : ", max(engs)
	print "Lowest pressure : ", min(press)
	print "Highest pressure : ", max(press)



	end = datetime.now()

	total = end - start

	print 
	print "  *** Completed in ", total.total_seconds(), " seconds *** "
	print
 

    