# ################################################################
#  
#    Active Particles on Curved Spaces (APCS)
#    
#	 Author: Silke Henkes
#   
#    ICSMB, Department of Physics
#    University of Aberdeen
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
#from inertia import *
from glob import glob
from datetime import *
from random import uniform 
from math import *
import numpy as np
import argparse
import scipy.spatial.distance as sd

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
#from matplotlib import rc
import matplotlib
from mpl_toolkits.mplot3d import Axes3D

# setting global parameters
#matplotlib.rcParams['text.usetex'] = 'true'
matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['axes.linewidth'] = 2
matplotlib.rcParams['xtick.major.size'] = 8
matplotlib.rcParams['ytick.major.size'] = 8
#matplotlib.rcParams['font.size']=40.0
#matplotlib.rcParams['legend.fontsize']=22.0
matplotlib.rcParams['font.size']=28
matplotlib.rcParams['legend.fontsize']=14


cdict = {'red':   [(0.0,  0.25, 0.25),
				   (0.3,  1.0, 1.0),
                   (0.5,  0.4, 0.0),
                   (1.0,  0.0, 0.0)],

         'green': [(0.0,  0.0, 0.0),
				   (0.25,  0.0, 0.5),
                   (0.5, 1.0, 1.0),
                   (0.75,  0.5, 0.0),
                   (1.0,  0.0, 0.0)],

         'blue':  [(0.0,  0.0, 0.0),
                   (0.5,  0.0, 0.0),
                   (0.7, 1.0, 1.0),
                   (1.0,  0.25, 0.25)]}


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
  
def compute_energy_and_pressure_rastko(r,k):
  eng = np.zeros(len(r))
  press = np.zeros(len(r))
  a = np.ones(len(r))
  dist = sd.cdist(r,r)
  for i in range(len(r)):
    for j in range(i+1,len(r)):
      dr = dist[i,j]
      if dr < 2:
        diff = 2.0-dr
        fact = 0.5*k*diff
        eng_val = fact*diff
        press_val = fact*dr
        eng[i] += eng_val
        eng[j] += eng_val
        press[i] += press_val
        press[j] += press_val
  return [eng, press]

def compute_energy_and_pressure(r,k,sigma,L):
	eng = np.zeros(len(r))
	press = np.zeros(len(r))
	stress = np.zeros((len(r),3,3))
	#dist = sd.cdist(r,r)
	dmax=4*sigma**2
	for i in range(len(r)):
	#for i in range(10):
		drval=r-r[i,:]
		# Periodic boundary conditions
		#x1=x1-self.L*np.round((x1-x0)/self.L)
		drval-=L*np.round(drval/L)
		dist=np.sum((drval)**2,axis=1)
		neighbours=[index for index,value in enumerate(dist) if value <dmax]
		neighbours.remove(i)
		dr=np.sqrt(dist[neighbours])
		diff=2.0-dr
		fact = 0.5*k*diff
		eng_val = fact*diff
		press_val = fact*dr
		# Stress (force moment) has to be element by element) r_a F_b = -k r_a dist_b 
		drvec=drval[neighbours]
		Fvec=k*((diff/dr).transpose()*(drvec).transpose()).transpose()
		for u in range(3):
			for v in range(3):
				stress[neighbours,u,v]+=0.5*drvec[:,u]*Fvec[:,v]
		eng[neighbours]+=eng_val
		press[neighbours]+=press_val
	return [eng, press, stress]

def getStats(f,stiffness,sigma,L,debug=False):
	
	print "Processing file : ", f
	data = ReadData(f)

	# rotate the system such that the principal direction of the moment of inertia
	# corresponding to the largest eiqenvalue align with the lab z axis
	x, y, z = np.array(data.data[data.keys['x']]), np.array(data.data[data.keys['y']]), np.array(data.data[data.keys['z']])
	vx, vy, vz = np.array(data.data[data.keys['vx']]), np.array(data.data[data.keys['vy']]), np.array(data.data[data.keys['vz']])
	nx, ny, nz = np.array(data.data[data.keys['nx']]), np.array(data.data[data.keys['ny']]), np.array(data.data[data.keys['nz']])
	rval = np.column_stack((x,y,z))
	vval = np.column_stack((vx,vy,vz))
	nval = np.column_stack((nx,ny,nz))
	ez = np.array([0,0,1])  # lab frame z-axis
	# The order parameter with v_0 still in it. Normalize in final polish
	orderparV=np.sum(vval,axis=0)/len(vval)
	orderpar=np.sum(nval,axis=0)/len(nval)
	print orderpar
	print orderparV
	direction = orderpar/np.linalg.norm(orderpar)
	directionV = orderparV/np.linalg.norm(orderparV)
	axisorth= np.cross(direction,directionV)
	axisval=np.linalg.norm(axisorth)
	alpha=np.arcsin(axisval)
	axisorth=axisorth/axisval
	axisnorm=np.cross(ez,directionV)
	axisnorm/=np.linalg.norm(axisnorm)
	
	print directionV
	print axisorth
	
	vel = np.sqrt(vval[:,0]**2 + vval[:,1]**2 + vval[:,2]**2)
	velnorm=((vval).transpose()/(vel).transpose()).transpose()
	
	eng, press,stress = compute_energy_and_pressure(rval,stiffness,sigma,L)
	print np.shape(stress)
	# Project the stresses into the e,theta,phi components. The rr component hast to be 0, and the r cross components
	# belong to the projection. So they are not all that interesting. 
	# We want the theta theta, theta phi, phi theta ant phi phi components (implicitly testing symmetries ...)
	# I give up on the notation. Stress is (N,3,3), the axes are (N,3). We want e_i sigma_ij e_j
	s_tt=np.sum(axisnorm*np.einsum('kij,j->ki',stress,axisnorm),axis=1)
	s_tp=np.sum(axisnorm*np.einsum('kij,j->ki',stress,directionV),axis=1)
	s_pt=np.sum(directionV*np.einsum('kij,j->ki',stress,axisnorm),axis=1)
	s_pp=np.sum(directionV*np.einsum('kij,j->ki',stress,directionV),axis=1)
	print np.shape(s_tt)
	# Mean density really makes no sense? Determined by the initial conditions in periodic boundary conditions.
	# I do not wish to set up artificial bins in a translationally invariant system
	vel_av=np.mean(vel)
	eng_av=np.mean(eng)
	press_av=np.mean(press)
	s_tt_av=np.mean(s_tt)
	s_tp_av=np.mean(s_tp)
	s_pt_av=np.mean(s_pt)
	s_pp_av=np.mean(s_pp)
	
	# Debugging output
	if debug==True:
		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')
		ax.scatter(rval[:,0], rval[:,1], rval[:,2], zdir='z', c='b')
		
	return [vel_av,eng_av,press_av,s_tt_av,s_tp_av,s_pt_av,s_pp_av,alpha,direction,directionV,orderpar,orderparV]


# Scripting version: Only execute if this is called as a script. Otherwise, it attempts to go through here when loading as a module 
# and throws errors because some arguments aren't defined
if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", type=str, help="Input file with particle velocity field")
	parser.add_argument("-o", "--output", type=str, default="profiles.dat", help="Output file (text file)")
	parser.add_argument("-k", "--k", type=float, default=1.0, help="soft potential strength")
	parser.add_argument("-L", "--system_l", type=float, default=28.2094791, help="system size for the flat system")
	parser.add_argument("-r", "--particle_r", type=float, default=1.0, help="radius of particle ")
	parser.add_argument("-s", "--skip", type=int, default=0, help="skip this many samples")
	args = parser.parse_args()

	print
	print "\tActive Particles on Curved Spaces (APCS)"
	print "\tPotential energy profile"
	print 
	print "\tSilke Henkes"
	print "\tUniversity of Aberdeen"
	print "\t(c) 2014"
	print "\t----------------------------------------------"
	print 
	print "\tInput : ", args.input
	print "\tOutput : ", args.output
	print "\tSpring constant : ", args.k
	print "\tSystem size : ", args.system_l
	print "\tSkip frames : ", args.skip
	print 

	start = datetime.now()
	
	files = sorted(glob(args.input+'*.dat'))[args.skip:]
	nini=800
	nfin=840
	usesnap=np.arange(nini,nfin)
	
	vel = np.zeros((len(files),)) 
	eng = np.zeros((len(files),)) 
	press = np.zeros((len(files),)) 
	s_tt = np.zeros((len(files),)) 
	s_tp = np.zeros((len(files),)) 
	s_pt = np.zeros((len(files),)) 
	s_pp = np.zeros((len(files),)) 
	alpha = np.zeros((len(files),)) 
	axis = np.zeros((len(files),3)) 
	orderpar = np.zeros((len(files),3)) 
	axisV = np.zeros((len(files),3)) 
	orderparV = np.zeros((len(files),3)) 
	tot = 0
	for f in files[nini:nfin]:
		[vel[tot],eng[tot],press[tot],s_tt[tot],s_tp[tot],s_pt[tot],s_pp[tot],alpha[tot],axis[tot,:],axisV[tot,:],orderpar[tot,:],orderparV[tot,:]]=getStats(f,args.k,args.particle_r,args.system_l)
		tot +=1
	
	#np.savetxt(args.output,  np.transpose(np.array([theta_bin,rho_profile,vel_profile,eng_profile,press_profile,alpha_profile,alpha_v_profile])),fmt='%12.6g', header='theta rho vel energy pressure alpha alpha_v')   # x,y,z equal sized 1D arrays
	
	end = datetime.now()

	total = end - start
	print "  *** Completed in ", total.total_seconds(), " seconds *** "
	
	plt.figure()
	plt.plot(usesnap,vel[:(nfin-nini)],'.-')
	plt.title('velocity')
	
	plt.figure()
	plt.plot(usesnap,eng[:(nfin-nini)],'.-')
	plt.title('energy')
	
	plt.figure()
	plt.plot(usesnap,press[:(nfin-nini)],'.-k')
	plt.plot(usesnap,s_tt[:(nfin-nini)],'.-r')
	plt.plot(usesnap,s_tp[:(nfin-nini)],'.-g')
	plt.plot(usesnap,s_pp[:(nfin-nini)],'.-b')
	plt.title('pressure')
	
	plt.figure()
	plt.plot(usesnap,alpha[:(nfin-nini)],'.-')
	plt.title('alpha')

	plt.show()