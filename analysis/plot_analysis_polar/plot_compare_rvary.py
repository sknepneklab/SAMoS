# ################################################################
#  
#    Active Particles on Curved Spaces (APCS)
#    
#    Author: Silke Henkes
#   
#    ICSMB, Department of Physics
#    University of Aberdeen
#    
#    (c) 2013, 2014
#    
#    This program cannot be used, copied, or modified without
#    explicit permission of the author.
#  
# ################################################################

#! /usr/bin/python

import sys, os, glob
import cPickle as pickle
import numpy as np
import scipy as sp
from scipy.io import savemat
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
#from matplotlib import rc
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import minimize
from scipy.optimize import fmin_cg

# setting global parameters
#matplotlib.rcParams['text.usetex'] = 'true'
matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['axes.linewidth'] = 2
matplotlib.rcParams['xtick.major.size'] = 8
matplotlib.rcParams['ytick.major.size'] = 8
matplotlib.rcParams['font.size']=20
matplotlib.rcParams['legend.fontsize']=14


cdict = {'red':   [(0.0,  0.75, 0.75),
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
                   

 # Computes spring part of the energy and the external energy, i.e. the integrated version of the active forceing (normal part)
# Spring energy:
def eng_extr(theta,k,v0):
    # The spring part, including *all* interactions
    EngSpring=0
    for l in range(len(theta)):
		theta_diff=2.0-R*np.abs(theta-theta[l])
		neighbours=[index for index,value in enumerate(theta_diff) if value >0]
		#print neighbours
		neighbours.remove(l)
		EngSpring += 0.5*k/R*sum(theta_diff[neighbours]**2)
    # The active forcing part
    EngSpring/=2.0
    EngExt = -slope*v0*sum(np.cos(theta))
    #print("Spring energy " + str(EngSpring) + " Active energy " + str(EngExt))
    return (EngSpring+EngExt)

# Gradient part of the energies. That is equivalent to the equations of motion in the paper
def eng_grad(theta,k,v0):
    #N = len(theta)
    der = np.zeros_like(theta)
    for l in range(len(theta)):
		dist=R*np.abs(theta-theta[l])
		neighbours=[index for index,value in enumerate(dist) if value <2.0]
		neighbours.remove(l)
		dr=dist[neighbours]
		diff=2.0-dr
		drvec=R*(theta[neighbours]-theta[l])
		#Fvec=k*((diff/dr).transpose()*(drvec).transpose()).transpose()
		#Fvec=k*(diff/dr)*drvec
		der[l]=k*np.sum(diff*drvec/dr)
    der += v0*np.sin(slope*theta)
    #print("Gradient vector" + str(der))
    return der
  
# Produce results as a function of v_0. I am using simulation parameters here

RList=[8,12,16,20,40,60]
testmap=LinearSegmentedColormap('test',cdict,N=len(RList))
v0=0.5
k = 1
sigma=1
# alpha-theta slope
slope=0.55
#slope=1.25
# Target initial overlap
delta=0.5
# Leading to spacing:
dphi=2*sigma-delta
# hence number of particles
w=0
plt.figure(figsize=(10,7),linewidth=2.0)
ax=plt.gca()
for R in RList:
	N=int(round(R*np.pi/dphi))
	print N
	
	# Initial passive pressure purely due to overlap (check dimensions)
	pini=2*sigma*k*delta
	# Arrange particles on a regular grid
	theta0=np.linspace(-np.pi/2,np.pi/2,N)
	#print theta0
	# Restrict positions to this half-circle
	bounds = [(-np.pi/2, np.pi/2)] * N
	#plt.figure(figsize=(10,7),linewidth=2.0)
	#ax=plt.gca() 
	# http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html
	# Minimize the energy associated to these theta positions, using BFGS. 
	# jac=eng_grad associates our gradient function to the jacobian, instead of a numerical estimate
	# The output is:
	# res.x: The solution array, that is the new positions of the particles, I assume?
	# For the rest: see here
	# http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.OptimizeResult.html#scipy.optimize.OptimizeResult
	# Plain minimization with gradient: only works for simple 1st neighbor energy
	#res = minimize(eng_extr, theta0, args=(k,v0), method='BFGS', jac=eng_grad)
	# Adding the constraints. Again, gradients and full energy do not work
	#res = minimize(eng_extr, theta0, args=(k,v0), method='L-BFGS-B', jac=eng_grad,bounds=bounds)
	# We also tried:
	#scipy.optimize.fmin_cg(f, x0, fprime=None, args=(), gtol=1e-05, norm=inf, epsilon=1.4901161193847656e-08, maxiter=None, full_output=0, disp=1, retall=0, callback=None)
	#res = fmin_cg(eng_extr,theta0,fprime=eng_grad,args=(k,v0))
	# Which works better (if slower) without the gradient; but doesn't allow for constraints
	# Using this here, very simple
	res = minimize(eng_extr, theta0, args=(k,v0), method='L-BFGS-B',bounds=bounds)
	theta=res.x 
	# Displacement field
	u=res.x-theta0
	# We are following the particles, not the absolute coordinate. This is a Lagrangian strain tensor
	# u_s = 2 du/dR + (du/dR)^2
	#dudR=(res.x[1:]-res.x[:-1])/dphi
	# or equivalently
	dudR=(u[1:]-u[:-1])/dphi
	strain0=2*dudR+dudR*dudR
	
	press=np.zeros((N,))

	for l in range(len(theta)):
		dist=R*np.abs(theta-theta[l])
		neighbours=[index for index,value in enumerate(dist) if value <2.0]
		neighbours.remove(l)
		dr=dist[neighbours]
		diff=2.0-dr
		drvec=R*(theta[neighbours]-theta[l])
		#Fvec=k*((diff/dr).transpose()*(drvec).transpose()).transpose()
		#Fvec=k*(diff/dr)*drvec
		press[l]=sigma*np.sum(np.abs(k*diff*drvec/dr))
	press/=2.0
	#plt.plot(theta[:,w],press,'--',color=testmap(w), linestyle='solid',label=vList[w])
	plt.plot(theta,press,'--',color=testmap(w), linestyle='--')
	w+=1
ax.set_xlabel(r'$\theta$')
#ax.set_ylabel(r'$p(\vartheta)$')
ax.set_ylabel('pressure')
#plt.legend(loc=1,ncol=2)
#plt.ylim(0,2)
#plt.xlim(-np.pi/2,np.pi/2)
plt.xlim(-1.5,1.5)
#plt.title('Pressure - chain calculation')
# This is the structured data file hierarchy. Replace as appropriate (do not go the Yaouen way and fully automatize ...)

basedir='/media/drogon/home/silke/Documents/Curved/Runs_rebuttal/'


JList=['1']
vList=['0.5']
#JList=['10']
#vList=['1']
nbin=180
rval=28.2094791
#RList=['5','8','12','16','20','40','60']
RList=['8','12','16','20','40','60']
nuList=['0.0', '0.01','0.05', '0.1','0.2', '0.3','0.5', '0.7','1', '2']

testmap=LinearSegmentedColormap('test',cdict,N=len(RList))
testmap2=LinearSegmentedColormap('test',cdict,N=len(nuList))
nstep=20000000
nsave=10000
nsnap=int(nstep/nsave)
skip=int(nsnap/3)
dt=0.001


# Profiles
# Set column to plot
usecolumn=4

profList=[r'$\theta$',r'$\rho$',r'$\sqrt{\langle v^2 \rangle}/v_0$','energy','pressure',r'$\Sigma_{\theta \theta}$',r'$\Sigma_{\theta \phi}$',r'$\Sigma_{\phi \theta}$',r'$\Sigma_{\phi \phi}$',r'$\alpha$',r'$\alpha_v$']
profName=['theta','rho','vrms','energy','pressure','stt','stp','spt','spp','alpha','alpha_v']

for i in range(len(vList)):
	#plt.figure(figsize=(10,7),linewidth=2.0)
	for j in range(len(RList)):
		print vList[i],RList[j]
		ax=plt.gca()
		outfile=basedir+'/profiles_v0' + vList[i] + '_R' + RList[j] + '.dat'
		outfile2=basedir + '/axis_v0'  + vList[i] + '_R' + RList[j] + '.dat'
		# header='theta rho vel energy pressure alpha alpha_v'
		profiles=sp.loadtxt(outfile, unpack=True)[:,:] 
		isdata=[index for index,value in enumerate(profiles[1,:]) if (value >0)]
		# Corrected
		## Forgot the normalization of rho by the angle band width
		#if usecolumn==1:
			#normz=2*np.pi*rval*abs(np.cos(profiles[0,:]))
			#profiles[1,isdata]=profiles[1,isdata]/normz[isdata]
			#profiles[1,:]/=np.mean(profiles[1,:])
		if usecolumn==2:
			plt.plot(profiles[0,isdata],profiles[usecolumn,isdata]/float(vList[i]),color=testmap(j), linestyle='solid',label=RList[j])
		else:
			plt.plot(profiles[0,isdata],profiles[usecolumn,isdata],color=testmap(j), linestyle='solid',label=RList[j])
	#if usecolumn<=8:
		#plt.ylim(0,1.25*profiles[usecolumn,nbin/2])
	if usecolumn==9:
		plt.plot(2*profiles[0,isdata],0.45*2*profiles[0,isdata],'k--')
		plt.text(0.5,0.1,'slope 0.45')
		
	#plt.xlim(-np.pi/2,np.pi/2)
	#plt.xlim(-1,1)
	plt.xlabel(profList[0]) 
	plt.ylabel(profList[usecolumn]) 
	plt.legend(loc=2,ncol=2)
	#plt.title('Interaction strength ' + r'$R=$' + RList[j])
	
	#filename=picsfolder + '/profile_' + profName[usecolumn] + '_J' + JList[j] +'.pdf'
	#plt.savefig(filename)
	
#for i in range(len(vList)):
    #plt.figure(figsize=(10,7),linewidth=2.0)
    #for j in range(len(nuList)):
        #print vList[i],nuList[j]
        #ax=plt.gca()
        #outfile=basedir+'/profiles_v0' + vList[i] + '_nu' + nuList[j] + '.dat'
        #outfile2=basedir + '/axis_v0'  + vList[i] + '_nu' + nuList[j] + '.dat'
        ## header='theta rho vel energy pressure alpha alpha_v'
        #profiles=sp.loadtxt(outfile, unpack=True)[:,:] 
        #isdata=[index for index,value in enumerate(profiles[1,:]) if (value >0)]
        ## Corrected
        ### Forgot the normalization of rho by the angle band width
        ##if usecolumn==1:
            ##normz=2*np.pi*rval*abs(np.cos(profiles[0,:]))
            ##profiles[1,isdata]=profiles[1,isdata]/normz[isdata]
            ##profiles[1,:]/=np.mean(profiles[1,:])
        #if usecolumn==2:
            #plt.plot(profiles[0,isdata],profiles[usecolumn,isdata]/float(vList[i]),color=testmap2(j), linestyle='solid',label=nuList[j])
        #else:
            #plt.plot(profiles[0,isdata],profiles[usecolumn,isdata],color=testmap2(j), linestyle='solid',label=nuList[j])
    ##if usecolumn<=8:
        ##plt.ylim(0,1.25*profiles[usecolumn,nbin/2])
    ##if usecolumn==9:
        ##if j==1:
            ##plt.plot(2*profiles[0,isdata],1.25*2*profiles[0,isdata],'k--')
            ##plt.text(0.5,0.75,'slope 1.25')
        ##if j==0:
            ##plt.plot(2*profiles[0,isdata],0.1*2*profiles[0,isdata],'k--')
            ##plt.text(0.5,0.3,'slope 0.1')
            ##plt.ylim(-0.4,0.4)
    ##plt.xlim(-np.pi/2,np.pi/2)
    #plt.xlim(-1.5,1.5)
    #plt.xlabel(profList[0]) 
    #plt.ylabel(profList[usecolumn]) 
    #plt.legend(loc=2,ncol=2)
    ##plt.title('Interaction strength ' + r'$R=$' + RList[j])
    
    ##filename=picsfolder + '/profile_' + profName[usecolumn] + '_J' + JList[j] +'.pdf'
    ##plt.savefig(filename)
	
plt.show()	
		
		
		
		
		
		