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

# Utility code for converting velocity file into POV-Ray script

from datetime import *
from random import uniform 
from math import *
import argparse
from read_data import *


class POVPrint:
  
  base_radius = 0.05
  cone_base_radius = 2.0  # multiple of the base radius
  base_hight = 0.7  # fraction of the total vector length that is taken by the cylindrical part 
  vec_scale = 1.0  #vector scale
  colour = [1.0,0.0,0.0]
  sphere_colour = [1.0,1.0,1.0]
  sphere_radius = 10.0
  particle_radius = 0.5
  
  
  def __init__(self,outfilename, data, to_plot='velocity'):
    self.outfilename = outfilename
    if not (data.keys.has_key('x') and data.keys.has_key('y') and data.keys.has_key('z')):
      raise Exception('Missing coordinates in the input data file.')
    if to_plot.lower() == 'velocity':
      if not (data.keys.has_key('vx') and data.keys.has_key('vy') and data.keys.has_key('vz')):
        raise Exception('Missing velocities in the input data file.')
      else:
        self.to_plot = ['vx','vy','vz']
    elif to_plot.lower() == 'director':
      if not (data.keys.has_key('nx') and data.keys.has_key('ny') and data.keys.has_key('nz')):
        raise Exception('Missing directors in the input data file.')
      else:
        self.to_plot = ['nx','ny','nz']
    else:
      raise Exception('Unknown plot type.')
    self.data = data
  
  def write(self,glass=True):
    self.out = open(self.outfilename,'w')
    self.out.write('global_settings{ assumed_gamma 1.0 }\n')
    #self.out.write('global_settings{ max_trace_level 30.0 }\n')
    self.out.write('#default{ finish{ ambient 0.1 diffuse 0.9 }}\n')
    self.out.write('#include "colors.inc"\n')
    self.out.write('#include "textures.inc"\n')
    self.out.write('#include "glass.inc"\n')
    self.out.write('#include "metals.inc"\n')
    self.out.write('#include "golds.inc"\n')
    self.out.write('#include "stones.inc"\n')
    self.out.write('#include "woods.inc"\n')
    self.out.write('#declare ParticleSolid = texture\n') 
    self.out.write('                     {\n')
    self.out.write('                       pigment\n') 
    self.out.write('                       {\n')
    self.out.write('                         color rgb<%f,%f,%f,0.15>\n' % (self.sphere_colour[0],self.sphere_colour[1],self.sphere_colour[2]))
    self.out.write('                      }\n')
    self.out.write('                       finish\n')
    self.out.write('                       { \n')
    self.out.write('                         reflection 0.05\n')
    self.out.write('                         specular 0.05 \n')
    self.out.write('                         ambient 0.5 \n')
    self.out.write('                         diffuse 0.5 \n')
    self.out.write('                         //roughness 1.0\n')
    self.out.write('                       }\n')
    self.out.write('                     }\n\n')
    self.out.write('#declare ArrowSolid = texture\n') 
    self.out.write('                     {\n')
    self.out.write('                       pigment\n') 
    self.out.write('                       {\n')
    self.out.write('                         color rgb<%f,%f,%f>\n' % (self.colour[0],self.colour[1],self.colour[2]))
    self.out.write('                       }\n')
    self.out.write('                       finish\n')
    self.out.write('                       { \n')
    self.out.write('                         reflection 0.05\n')
    self.out.write('                         specular 0.05 \n')
    self.out.write('                         ambient 0.5 \n')
    self.out.write('                         diffuse 0.5 \n')
    self.out.write('                         //roughness 0.02\n')
    self.out.write('                       }\n')
    self.out.write('                     }\n\n')
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
    if glass:
      R = self.sphere_radius
      self.out.write('sphere { <0,0,0>,%f\n' % R)
      self.out.write('         pigment { rgbf <0.8,0.9,1,0.9> } // A blue-tinted glass\n')
      self.out.write('         finish { phong 0.9 phong_size 40  // A highlight\n')
      self.out.write('                  reflection 0.2  // Glass reflects a bit\n')
      self.out.write('                }\n')
      self.out.write('         interior { ior 1.0  }// Glass refraction\n')
      self.out.write('        }\n')
    self.__write_arrows()
    self.out.close()

  def __write_arrows(self):
    X, Y, Z = self.data.data[self.data.keys['x']], self.data.data[self.data.keys['y']], self.data.data[self.data.keys['z']]
    VX, VY, VZ = self.data.data[self.data.keys[self.to_plot[0]]], self.data.data[self.data.keys[self.to_plot[1]]], self.data.data[self.data.keys[self.to_plot[2]]]
    for (x,y,z,vx,vy,vz) in zip(X,Y,Z,VX,VY,VZ):
      vvx, vvy, vvz = self.vec_scale*vx, self.vec_scale*vy, self.vec_scale*vz
      x1, y1, z1 = x - 0.5*vvx, y - 0.5*vvy, z - 0.5*vvz 
      x2, y2, z2 = x + 0.5*vvx, y + 0.5*vvy, z + 0.5*vvz 
      cx = (1.0 - self.base_hight)*x1 + self.base_hight*x2
      cy = (1.0 - self.base_hight)*y1 + self.base_hight*y2
      cz = (1.0 - self.base_hight)*z1 + self.base_hight*z2
      hight = sqrt((x2-x1)**2+(y2-y1)**2+(y2-y1)**2)
      if self.particle_radius > 0:
        self.out.write('sphere { <%f,%f,%f>,%f material {texture {ParticleSolid} interior { ior 1.0 } } }\n' % (x,y,z,self.particle_radius))      
      if (hight != 0.0):
        self.out.write('object {\n')
        self.out.write('       union {\n')
        self.out.write('              cylinder { <%f,%f,%f>,<%f,%f,%f>,%f texture { ArrowSolid } }\n' % (x1,y1,z1,cx,cy,cz,self.base_radius))
        self.out.write('              cone { <%f,%f,%f>,%f,<%f,%f,%f>,%f texture { ArrowSolid } }\n' % (cx,cy,cz,self.base_radius*self.cone_base_radius,x2,y2,z2,0.0))
        self.out.write('             }\n')
        self.out.write('        }\n')

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="Input file with particle data")
parser.add_argument("-o", "--output", type=str, default="out.pov", help="Output file (POV-Ray scene script)")
parser.add_argument("-s", "--scale", type=float, default=1.0, help="velocity scale")
parser.add_argument("-r", "--base_radius", type=float, default=0.05, help="arrow base radius")
parser.add_argument("-p", "--particle_radius", type=float, default=0.5, help="radius of each particle")
parser.add_argument("-c", "--cone", type=float, default=2.0, help="cone base size (multiple of arrow base radius)")
parser.add_argument("-H", "--hight", type=float, default=0.7, help="fraction of the arrow in the base (arrow tip will be (1-hight) long) ")
parser.add_argument("-C", "--cone_colour", type=float, nargs=3, default=[1.0,0.0,0.0], help="arrow colour ")
parser.add_argument("-P", "--particle_colour", type=float, nargs=3, default=[1.0,1.0,1.0], help="particle colour ")
parser.add_argument("-S", "--sphere", action='store_true', default=False, help="include glass sphere")
parser.add_argument("-R", "--sphere_r", type=float, default=3.0, help="radius of sphere for spherical system")
parser.add_argument("-T", "--to_plot", type=str, default='director', help="type of output quantity (velocity or director)")
args = parser.parse_args()

print
print "\tActive Particles on Curved Spaces (APCS)"
print "\tConvert velocity field into POV-Ray scene"
print 
print "\tRastko Sknepnek"
print "\tUniversity of Dundee"
print "\t(c) 2013"
print "\t----------------------------------------------"
print 
print "\tInput : ", args.input
print "\tOutput : ", args.output
print "\tScale velocity : ", args.scale
print "\tBase arrow radius : ", args.base_radius
print "\tCone base radius : ", args.base_radius*args.cone
print "\tArrow base fraction : ", args.hight
print "\tArrow colour : ", args.cone_colour
print "\tParticle colour : ", args.particle_colour
print "\tParticle radius : ", args.particle_radius
print "\tOutput particle : ", args.to_plot
print 

start = datetime.now()

data = ReadData(args.input)
p = POVPrint(args.output,data,args.to_plot)
p.vec_scale = args.scale
p.base_radius = args.base_radius
p.cone_base_radius = args.cone
p.base_hight = args.hight
p.colour = args.cone_colour
p.sphere_radius = args.sphere_r
p.particle_radius = args.particle_radius
p.sphere_colour = args.particle_colour
p.write(args.sphere)


end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print
