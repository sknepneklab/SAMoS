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

# Utility code for computing potential energy profile averaged in the 
# azimuthal direction 


from read_data import *
from op import *
from inertia import *
from glob import glob
from datetime import *
from random import uniform 
from math import *
import numpy as np
import argparse
import scipy.spatial.distance as sd

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


def appendSpherical_np(xyz):
    R = np.sqrt(xyz[:,0]**2 + xyz[:,1]**2 + xyz[:,2]**2)
    angle = np.arccos(xyz[:,2]/R)
    return np.degrees(angle)
    #ptsnew = np.hstack((xyz, np.zeros(xyz.shape)))
    #xy = xyz[:,0]**2 + xyz[:,1]**2
    #ptsnew[:,3] = np.sqrt(xy + xyz[:,2]**2)
    #ptsnew[:,4] = np.arctan2(np.sqrt(xy), xyz[:,2]) # for elevation angle defined from Z-axis down
    ##ptsnew[:,4] = np.arctan2(xyz[:,2], np.sqrt(xy)) # for elevation angle defined from XY-plane up
    #ptsnew[:,5] = np.arctan2(xyz[:,1], xyz[:,0])
    #return ptsnew

def compute_energy_and_pressure(r,k):
  eng = np.zeros(len(r))
  press = np.zeros(len(r))
  overlap = np.zeros(len(r))
  a = np.ones(len(r))
  dist = sd.cdist(r,r)
  for i in range(len(r)):
    for j in range(i+1,len(r)):
      dr = dist[i,j]
      if dr < 2.0:
        diff = 2.0-dr
        fact = 0.5*k*diff
        eng_val = fact*diff
        press_val = fact*dr
        eng[i] += eng_val
        eng[j] += eng_val
        press[i] += press_val
        press[j] += press_val
        overlap[i] += 1
        overlap[j] += 1
  return [eng, press, overlap]



parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="Input file with particle velocity field")
parser.add_argument("-o", "--output", type=str, default="out", help="Output file (POV-Ray scene script)")
parser.add_argument("-k", "--k", type=float, default=1.0, help="soft potential strength")
parser.add_argument("-R", "--sphere_r", type=float, default=10.0, help="radius of sphere for spherical system")
parser.add_argument("-n", "--bin", type=int, default=25, help="number of bins for angle average")
parser.add_argument("-s", "--skip", type=int, default=0, help="skip this many samples")
args = parser.parse_args()

print
print "\tActive Particles on Curved Spaces (APCS)"
print "\tPotential energy profile"
print 
print "\tRastko Sknepnek"
print "\tUniversity of Dundee"
print "\t(c) 2013"
print "\t----------------------------------------------"
print 
print "\tInput : ", args.input
print "\tOutput : ", args.output
print "\tSpring constant : ", args.k
print "\tRadius of the sphere : ", args.sphere_r
print "\tNumber of angle average bins : ", args.bin
print "\tSkip frames : ", args.skip
print 


start = datetime.now()

files = sorted(glob(args.input+'*.dat'))[args.skip:]
print files
#print "Reading data..."
#data = ReadData(args.input)

ez = np.array([0,0,1])  # lab frame z-axis

tot_avg_theta = [0 for i in range(args.bin+1)]
tot_avg_v = [0 for i in range(args.bin+1)]
tot_avg_press = [0 for i in range(args.bin+1)]
tot_avg_overlap = [0 for i in range(args.bin+1)]
tot_nval = [0 for i in range(args.bin+1)]
tot = 0
for f in files:
  print "Processing file : ", f
  data = ReadData(f)
  inertia = Inertia(data)
  I = inertia.compute()   # Compute moment of inertia
  direction = I[1][:,0]   # Presumably largest component is along z-axis

  # rotate the system such that the principal direction of the moment of inertia
  # corresponding to the largest eiqenvalue align with the lab z axis
  axis = np.cross(direction,ez)
  axis = axis/np.linalg.norm(axis)
  rot_angle = acos(np.dot(direction,ez))
  x, y, z = np.array(data.data[data.keys['x']]), np.array(data.data[data.keys['y']]), np.array(data.data[data.keys['z']])
  vx, vy, vz = np.array(data.data[data.keys['vx']]), np.array(data.data[data.keys['vy']]), np.array(data.data[data.keys['vz']])
  r = np.column_stack((x,y,z))
  n = np.empty(np.shape(r))
  n[:,0] = axis[0]
  n[:,1] = axis[1]
  n[:,2] = axis[2]
  vrot = rotate_vectorial(r,n,-rot_angle)
  #ptsnew = appendSpherical_np(vrot)
  theta = appendSpherical_np(vrot)
  print "Computing potential energy..."
  eng, press, overlap = compute_energy_and_pressure(vrot,args.k)
  print "Computing velocity magnitude..."
  vel = np.sqrt(vx**2 + vy**2 + vz**2)
  print "Computing angular average..."
  #theta = ptsnew[:,4]
  t_max, t_min = 180.0, 0.0#max(theta), min(theta)
  dtheta = (t_max-t_min)/args.bin
  print 't_min = ', t_min
  print 't_max = ', t_max
  print 'dt = ', dtheta
 
  avg_theta = [0 for i in range(args.bin+1)]
  avg_v = [0 for i in range(args.bin+1)]
  avg_press = [0 for i in range(args.bin+1)]
  avg_overlap = [0 for i in range(args.bin+1)]
  nval = [0 for i in range(args.bin+1)]
  for i in range(len(eng)):
    #idx = int(round(((theta[i]-min(theta))/dtheta)))
    idx = int(round((theta[i]-t_min)/dtheta))
    avg_theta[idx] += eng[i]
    avg_v[idx] += vel[i]
    avg_press[idx] += press[i]
    avg_overlap[idx] += overlap[i]
    nval[idx] += 1
    
  
  for idx in xrange(len(avg_theta)):
    if nval[idx] != 0: 
      tot_avg_theta[idx] += avg_theta[idx]/float(nval[idx])
      tot_avg_v[idx] += avg_v[idx]/float(nval[idx])
      tot_avg_press[idx] += avg_press[idx]/float(nval[idx])
      tot_avg_overlap[idx] += avg_overlap[idx]/float(nval[idx])
    tot_nval[idx] += nval[idx]  
    
  tot += 1

out = open(args.output,'w')
out.write('# Total potential energy as a function of polar angle (measured from equator)\n')
out.write('# Generated on: %s\n' % str(datetime.now()))
out.write('# angle(deg)   <E>  <v>  <P>  <density>  <A>  <overlap>\n')

for i in range(len(nval)):
  E = tot_avg_theta[i]/tot
  V = tot_avg_v[i]/tot
  P = tot_avg_press[i]/tot
  O = tot_avg_overlap[i]/tot
  theta = t_min+i*dtheta
  A = 2.0*pi*args.sphere_r**2 * sin(radians(theta))*radians(dtheta)
  if A != 0.0:
    rho = (tot_nval[i]/float(tot))/A
  else:
    rho = 0.0
  out.write('%f %f  %f  %f  %f  %f  %f\n' % (theta-90, E, V, P, rho, A, O))
  #out.write('%f %f  %f\n' % (degrees(t_min+i*dtheta-0.5*pi), tot_avg_theta[i]/tot, tot_avg_v[i]/tot))

out.close()

end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print

