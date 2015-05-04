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

# Utility code for computing time correlation function for
# velocity (or director) an a sphere for a series of velocity/director files.

import sys
import argparse
from read_data import *
import math
from glob import glob
import numpy as np
import numpy.linalg as la
from datetime import *


def rotate_matrix_vectorial(axis,theta):
    axlen=np.sqrt(axis[:,0]**2+axis[:,1]**2+axis[:,2]**2)
    #print axlen
    axis[:,0]=axis[:,0]/axlen
    axis[:,1]=axis[:,1]/axlen
    axis[:,2]=axis[:,2]/axlen
    a=np.cos(theta/2)
    b=-axis[:,0]*np.sin(theta/2)
    c=-axis[:,1]*np.sin(theta/2)
    d=-axis[:,2]*np.sin(theta/2)
    rotmat=np.empty((len(axis[:,0]),3,3))
    rotmat[:,0,0]=a*a+b*b-c*c-d*d
    rotmat[:,0,1]=2*(b*c-a*d)
    rotmat[:,0,2]=2*(b*d+a*c)
    rotmat[:,1,0]=2*(b*c+a*d)
    rotmat[:,1,1]=a*a+c*c-b*b-d*d
    rotmat[:,1,2]=2*(c*d-a*b)
    rotmat[:,2,0]=2*(b*d-a*c)
    rotmat[:,2,1]=2*(c*d+a*b)
    rotmat[:,2,2]=a*a+d*d-b*b-c*c
    return rotmat
  
def rotate_vectorial(v,n,phi):
    vrot=np.empty(np.shape(v))
    np.shape(vrot)
    rotmat=rotate_matrix_vectorial(n,phi)
    np.shape(rotmat)
    vrot[:,0]=rotmat[:,0,0]*v[:,0]+rotmat[:,0,1]*v[:,1]+rotmat[:,0,2]*v[:,2]
    vrot[:,1]=rotmat[:,1,0]*v[:,0]+rotmat[:,1,1]*v[:,1]+rotmat[:,1,2]*v[:,2]
    vrot[:,2]=rotmat[:,2,0]*v[:,0]+rotmat[:,2,1]*v[:,1]+rotmat[:,2,2]*v[:,2]
    return vrot
    

def rotation_matrix(axis,theta):
    axis = axis/math.sqrt(np.dot(axis,axis))
    a = math.cos(theta/2)
    b,c,d = -axis*math.sin(theta/2)
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])

def rotate(v,n,phi):
    return np.dot(rotation_matrix(n,phi),v)
  
def rotate2(v,n,phi):
    x, y, z = v
    c = math.cos(phi)
    s = math.sin(phi)
    len_n = math.sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2])
    U, V, W = n[0]/len_n, n[1]/len_n, n[2]/len_n
    vx = U*(U*x+V*y+W*z)*(1-c) + x*c + (-W*y+V*z)*s
    vy = V*(U*x+V*y+W*z)*(1-c) + y*c + ( W*x-U*z)*s
    vz = W*(U*x+V*y+W*z)*(1-c) + z*c + (-V*x+U*y)*s
    return [vx, vy, vz]


parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", type=str, default='time_correl.dat', help="contains time correlation data")
parser.add_argument("-i", "--input", type=str, help="base name of the input files")
parser.add_argument("-s", "--skip", type=int, default=0, help="skip this many samples")
parser.add_argument("-d", "--director", action="store_true", help="compute time correlations of the director (velocity is default)")
args = parser.parse_args()

print
print "\tActive Particles on Curved Spaces (APCS)"
print "\tComputes time correlation function of velocity"
print "\tor director field on a sphere"
print 
print "\tRastko Sknepnek"
print "\tUniversity of Dundee"
print "\t(c) 2014"
print "\t----------------------------------------------"
print
print "\tInput files : ", args.input
print "\tOutput files : ", args.output
print "\tSkip frames : ", args.skip
if args.director:
  print "\tCompute time correlation of the director field."
print

start = datetime.now()

files = sorted(glob(args.input+'*.dat'))[args.skip:]

if args.director:
  qx, qy, qz = 'nx', 'ny', 'nz'
else:
  qx, qy, qz = 'vx', 'vy', 'vz'

vals = []
idx = 0
for f in files:
  print "Processing file : ", f
  data = ReadData(f)
  frame = []
  if data.has_header:
    for (x, y, z, ax, ay, az) in zip(data.data[data.keys['x']],data.data[data.keys['y']],data.data[data.keys['z']],data.data[data.keys[qx]],data.data[data.keys[qy]],data.data[data.keys[qz]]):
      frame.append([x,y,z,ax,ay,az])
    vals.append([idx,frame])  
  idx += 1

time_file = open(args.output,'w')

print "Computing averages..."


for tau in range(1,len(vals)-1):
  print "Processing tau : ", tau
  avg_tau = 0.0
  for t in range(0,len(vals)-tau):
    frame1 = np.array(vals[t][1])
    frame2 = np.array(vals[t+tau][1])
    r1=frame1[:,0:3]
    r2=frame2[:,0:3]
    a1=frame1[:,3:6]
    a2=frame2[:,3:6]
    r2_x_r1=np.cross(r2,r1)
    len_r2_x_r1=np.sqrt(r2_x_r1[:,0]**2+r2_x_r1[:,1]**2+r2_x_r1[:,2]**2)
    lenr1=np.sqrt(r1[:,0]**2+r1[:,1]**2+r1[:,2]**2)
    lenr2=np.sqrt(r2[:,0]**2+r2[:,1]**2+r2[:,2]**2)
    dot_r1r2=r1[:,0]*r2[:,0]+r1[:,1]*r2[:,1]+r1[:,2]*r2[:,2]
    n=np.empty(np.shape(r1))
    n[:,0] = r2_x_r1[:,0]/len_r2_x_r1
    n[:,1] = r2_x_r1[:,1]/len_r2_x_r1
    n[:,2] = r2_x_r1[:,2]/len_r2_x_r1
    phi = np.arccos(dot_r1r2/(lenr1*lenr2))
    a2trans=rotate_vectorial(a2,n,-phi)
    avg=np.sum(a1[:,0]*a2trans[:,0]+a1[:,1]*a2trans[:,1]+a1[:,2]*a2trans[:,2])
    avg/=len(r1[:,0])
    #print avg
    avg_tau+=avg
  time_file.write('%d   %f\n' % (tau,avg_tau/(len(vals)-tau)))
  time_file.flush()
    
time_file.close()
  


end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print
