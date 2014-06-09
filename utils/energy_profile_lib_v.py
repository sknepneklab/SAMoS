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

def compute_energy_and_pressure(r,k,sigma):
	eng = np.zeros(len(r))
	press = np.zeros(len(r))
	stress = np.zeros((len(r),3,3))
	#dist = sd.cdist(r,r)
	dmax=4*sigma**2
	for i in range(len(r)):
	#for i in range(10):
		dist=np.sum((r-r[i,:])**2,axis=1)
		neighbours=[index for index,value in enumerate(dist) if value <dmax]
		neighbours.remove(i)
		dr=np.sqrt(dist[neighbours])
		diff=2.0-dr
		fact = 0.5*k*diff
		eng_val = fact*diff
		press_val = fact*dr
		# Stress (force moment) has to be element by element) r_a F_b = -k r_a dist_b 
		drvec=r[neighbours,:]-r[i,:]
		Fvec=k*((diff/dr).transpose()*(drvec).transpose()).transpose()
		for u in range(3):
			for v in range(3):
				stress[neighbours,u,v]+=0.5*drvec[:,u]*Fvec[:,v]
		eng[neighbours]+=eng_val
		press[neighbours]+=press_val
	return [eng, press, stress]

def getProfiles(f,nbin,radius,stiffness,sigma,debug=False):
	
	print "Processing file : ", f
	data = ReadData(f)
	#inertia = Inertia(data)
	#I = inertia.compute()   # Compute moment of inertia
	#direction = I[1][:,0]   # Presumably smallest component is along z-axis. 

	# rotate the system such that the principal direction of the moment of inertia
	# corresponding to the largest eiqenvalue align with the lab z axis
	x, y, z = np.array(data.data[data.keys['x']]), np.array(data.data[data.keys['y']]), np.array(data.data[data.keys['z']])
	vx, vy, vz = np.array(data.data[data.keys['vx']]), np.array(data.data[data.keys['vy']]), np.array(data.data[data.keys['vz']])
	nx, ny, nz = np.array(data.data[data.keys['nx']]), np.array(data.data[data.keys['ny']]), np.array(data.data[data.keys['nz']])
	rval = np.column_stack((x,y,z))
	vval = np.column_stack((vx,vy,vz))
	nval = np.column_stack((nx,ny,nz))
	ez = np.array([0,0,1])  # lab frame z-axis
	# Simply get the axis as the mean crossproduct or r and v; assuming alignment. This should also not flip.
	direction=np.sum(np.cross(rval,vval),axis=0)
	orderpar=direction/len(x)
	print orderpar
	direction = direction/np.linalg.norm(direction)
	axis = np.cross(direction,ez)
	axis = axis/np.linalg.norm(axis)
	rot_angle = acos(np.dot(direction,ez))
	axis0 = np.empty(np.shape(rval))
	axis0[:,0] = axis[0]
	axis0[:,1] = axis[1]
	axis0[:,2] = axis[2]
	rval = rotate_vectorial(rval,axis0,-rot_angle)
	vval = rotate_vectorial(vval,axis0,-rot_angle)
	nval = rotate_vectorial(nval,axis0,-rot_angle)
	nval=((nval).transpose()/(np.sqrt(np.sum(nval**2,axis=1))).transpose()).transpose()
	rhat=((rval).transpose()/(np.sqrt(np.sum(rval**2,axis=1))).transpose()).transpose()
	vel = np.sqrt(vval[:,0]**2 + vval[:,1]**2 + vval[:,2]**2)
	velnorm=((vval).transpose()/(vel).transpose()).transpose()
	
	# Determine the Euler angles, essentially. Find theta and phi for each particle, and use it to compute alpha and stress components
	# Angle theta with the z axis. arccos is between 0 and pi, so that's ok already
	theta=np.arccos(rhat[:,2])
	# From the euler angles: rx = sin theta cos phi
	# Choosing correct quadrant through the sign of ry=sin theta sin phi
	phi=np.sign(rhat[:,1]/(np.sin(theta)))*np.arccos(rhat[:,0]/(np.sin(theta)))
	# The other two of our trio of local coordinate vectors
	etheta = np.empty(np.shape(rval))
	etheta[:,0]=np.cos(theta)*np.cos(phi)
	etheta[:,1]=np.cos(theta)*np.sin(phi)
	etheta[:,2]=-np.sin(theta)
	ephi=np.empty(np.shape(rval))
	ephi[:,0]=-np.sin(phi)
	ephi[:,1]=np.cos(phi)
	ephi[:,2]=0
	# Alpha, the angle between the local polarity and the equator; here represented by ephi
	alpha=-np.arcsin(np.sum(nval*etheta, axis=1))
	# Same thing for the velocity
	# No - add pi/2 to get something that does not add up to zero 
	alpha_v=np.arccos(np.sum(velnorm*etheta, axis=1))
	
	eng, press,stress = compute_energy_and_pressure(rval,stiffness,sigma)
	# Project the stresses into the e,theta,phi components. The rr component hast to be 0, and the r cross components
	# belong to the projection. So they are not all that interesting. 
	# We want the theta theta, theta phi, phi theta ant phi phi components (implicitly testing symmetries ...)
	# I give up on the notation. Stress is (N,3,3), the axes are (N,3). We want e_i sigma_ij e_j
	s_tt=np.sum(etheta*np.einsum('kij,kj->ki',stress,etheta),axis=1)
	s_tp=np.sum(etheta*np.einsum('...ij,...j->...i',stress,ephi),axis=1)
	s_pt=np.sum(ephi*np.einsum('...ij,...j->...i',stress,etheta),axis=1)
	s_pp=np.sum(ephi*np.einsum('...ij,...j->...i',stress,ephi),axis=1)
	
	# Setting up the binning. I changed this to go from -pi/2 to pi/2 consistently. This maybe makes less pretty pictures,
	# but the edges are going to be a lot cleaner. Also only one bin to handle accross multiple v0/J.
	# Can always rebin to less resolution if necessary
	# Position angle with the z axis
	theta_bin=np.linspace(0,np.pi,nbin+1)
	dtheta=theta_bin[1]-theta_bin[0]
	theta_out=theta_bin[:nbin]+dtheta/2-np.pi/2
	
	rho_profile, bin_edges = np.histogram(theta, bins=theta_bin,density=True)
	isdata=[index for index,value in enumerate(rho_profile) if (value >0)]
	normz=2*np.pi*radius*abs(np.cos(theta_out))
	rho_profile[isdata]=rho_profile[isdata]/normz[isdata]
	rho_profile/=np.mean(rho_profile)
	vel_profile=np.zeros(np.shape(rho_profile))
	eng_profile=np.zeros(np.shape(rho_profile))
	press_profile=np.zeros(np.shape(rho_profile))
	s_tt_profile=np.zeros(np.shape(rho_profile))
	s_tp_profile=np.zeros(np.shape(rho_profile))
	s_pt_profile=np.zeros(np.shape(rho_profile))
	s_pp_profile=np.zeros(np.shape(rho_profile))
	alpha_profile=np.zeros(np.shape(rho_profile))
	alpha_v_profile=np.zeros(np.shape(rho_profile))
	for idx in range(nbin):
		inbin=[index for index,value in enumerate(theta) if (value >= theta_bin[idx]  and value<=theta_bin[idx+1])]
		#print len(inbin)
		if len(inbin)>0:
			vel_profile[idx]=np.mean(vel[inbin])
			eng_profile[idx]=np.mean(eng[inbin])
			press_profile[idx]=np.mean(press[inbin])
			s_tt_profile[idx]=np.mean(s_tt[inbin])
			s_tp_profile[idx]=np.mean(s_tp[inbin])
			s_pt_profile[idx]=np.mean(s_pt[inbin])
			s_pp_profile[idx]=np.mean(s_pp[inbin])
			alpha_profile[idx]=np.mean(alpha[inbin])
			alpha_v_profile[idx]=np.mean(alpha_v[inbin])
	
	
	
	# Debugging output
	if debug==True:
		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')
		ax.scatter(rval[:,0], rval[:,1], rval[:,2], zdir='z', c='b')
		
	return [theta_out,rho_profile,vel_profile,eng_profile,press_profile,s_tt_profile,s_tp_profile,s_pt_profile,s_pp_profile,alpha_profile,alpha_v_profile,direction,orderpar]


# Scripting version: Only execute if this is called as a script. Otherwise, it attempts to go through here when loading as a module 
# and throws errors because some arguments aren't defined
if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", type=str, help="Input file with particle velocity field")
	parser.add_argument("-o", "--output", type=str, default="profiles.dat", help="Output file (text file)")
	parser.add_argument("-k", "--k", type=float, default=1.0, help="soft potential strength")
	parser.add_argument("-R", "--sphere_r", type=float, default=28.2094791, help="radius of sphere for spherical system")
	parser.add_argument("-r", "--particle_r", type=float, default=1.0, help="radius of particle ")
	parser.add_argument("-n", "--bin", type=int, default=100, help="number of bins for angle average")
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
	
	nbin=args.bin
	files = sorted(glob(args.input+'*.dat'))[args.skip:]
	
	rho_profile =np.zeros((nbin,)) 
	vel_profile =np.zeros((nbin,)) 
	eng_profile = np.zeros((nbin,)) 
	press_profile = np.zeros((nbin,)) 
	s_tt_profile = np.zeros((nbin,)) 
	s_tp_profile = np.zeros((nbin,)) 
	s_pt_profile = np.zeros((nbin,)) 
	s_pp_profile = np.zeros((nbin,)) 
	alpha_profile = np.zeros((nbin,)) 
	alpha_v_profile = np.zeros((nbin,)) 
	axis = np.zeros((len(files),3)) 
	orderpar = np.zeros((len(files),3)) 
	iscount = np.zeros((nbin,)) 
	tot = 0
	for f in files[800:801]:
		[theta_bin,rho_profile0,vel_profile0,eng_profile0,press_profile0,s_tt_profile0,s_tp_profile0,s_pt_profile0,s_pp_profile0,alpha_profile0,alpha_v_profile0,axis[tot,:],orderpar[tot,:]]=getProfiles(f,args.bin,args.sphere_r,args.k,args.particle_r)	
		isparticles=np.array([index for index,value in enumerate(rho_profile0) if (value >0)])
		iscount[isparticles]+=1
		rho_profile[isparticles]+=rho_profile0[isparticles]
		vel_profile[isparticles]+=vel_profile0[isparticles]
		eng_profile[isparticles]+=eng_profile0[isparticles]
		press_profile[isparticles]+=press_profile0[isparticles]
		s_tt_profile[isparticles]+=s_tt_profile0[isparticles]
		s_tp_profile[isparticles]+=s_tp_profile0[isparticles]
		s_pt_profile[isparticles]+=s_pt_profile0[isparticles]
		s_pp_profile[isparticles]+=s_pp_profile0[isparticles]
		alpha_profile[isparticles]+=alpha_profile0[isparticles]
		alpha_v_profile[isparticles]+=alpha_v_profile0[isparticles]
		tot +=1
	issomething=[index for index,value in enumerate(iscount) if (value >0)]
	rho_profile[issomething]/=iscount[issomething]
	vel_profile[issomething]/=iscount[issomething]
	eng_profile[issomething]/=iscount[issomething]
	press_profile[issomething]/=iscount[issomething]
	press_profile[issomething]/=iscount[issomething]
	s_tt_profile[issomething]/=iscount[issomething]
	s_tp_profile[issomething]/=iscount[issomething]
	s_pt_profile[issomething]/=iscount[issomething]
	s_pp_profile[issomething]/=iscount[issomething]
	alpha_profile[issomething]/=iscount[issomething]
	alpha_v_profile[issomething]/=iscount[issomething]
	
	#np.savetxt(args.output,  np.transpose(np.array([theta_bin,rho_profile,vel_profile,eng_profile,press_profile,alpha_profile,alpha_v_profile])),fmt='%12.6g', header='theta rho vel energy pressure alpha alpha_v')   # x,y,z equal sized 1D arrays
	
	end = datetime.now()

	total = end - start
	print "  *** Completed in ", total.total_seconds(), " seconds *** "
	
	plt.figure()
	plt.plot(theta_bin[issomething],rho_profile[issomething],'.-')
	plt.title('density')
	
	plt.figure()
	plt.plot(theta_bin[issomething],vel_profile[issomething],'.-')
	plt.title('velocity')
	
	plt.figure()
	plt.plot(theta_bin[issomething],eng_profile[issomething],'.-')
	plt.title('energy')
	
	plt.figure()
	plt.plot(theta_bin[issomething],press_profile[issomething],'.-k')
	plt.plot(theta_bin[issomething],s_tt_profile[issomething],'.-r')
	plt.plot(theta_bin[issomething],s_tp_profile[issomething],'.-g')
	plt.plot(theta_bin[issomething],s_pp_profile[issomething],'.-b')
	plt.title('pressure')
	
	plt.figure()
	plt.plot(theta_bin[issomething],alpha_profile[issomething],'.-')
	plt.title('alpha')
	
	plt.figure()
	plt.plot(theta_bin[issomething],alpha_v_profile[issomething],'.-')
	plt.title('alpha_v')

	plt.show()