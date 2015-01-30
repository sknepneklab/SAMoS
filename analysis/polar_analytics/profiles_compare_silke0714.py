# -*- coding: utf-8 -*-
"""
Created on Wed Jun  4 13:40:04 2014

@author: sknepnek
"""
# Extended and used for comparison to numerics starting Sat 28 Jun 2014
# Silke Henkes

#! /usr/bin/python
import numpy as np
from scipy.optimize import minimize
from scipy.optimize import fmin_cg
import matplotlib.pyplot as plt
import math
from matplotlib.colors import LinearSegmentedColormap
#from matplotlib import rc
import matplotlib

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
R = 28.2094791774
k = 1
sigma=1
# alpha-theta slope
slopeList=[0.25,0.6,1.25]
JList=[10,1,0.1]
# Target initial overlap
delta=0.5
# Leading to spacing:
dphi=2*sigma-delta
# hence number of particles
N=int(round(R*np.pi/dphi))
print N
v0List=[0.1,0.5]
testmap=LinearSegmentedColormap('test',cdict,N=len(JList))
# Initial passive pressure purely due to overlap (check dimensions)
pini=2*sigma*k*delta
# Arrange particles on a regular grid
theta0=np.linspace(-np.pi/2,np.pi/2,N)
#print theta0
# Restrict positions to this half-circle
bounds = [(-np.pi/2, np.pi/2)] * N
#theta=np.zeros((N,len(vList)))
#strain=np.zeros((N-1,len(vList)))

plt.figure(figsize=(10,7),linewidth=2.0)
ax=plt.gca() 
for v0 in v0List:
	w=0
	for slope in slopeList:
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
plt.xlim(-np.pi/2,np.pi/2)
#plt.title('Pressure - chain calculation')



#plt.figure(figsize=(10,7),linewidth=2.0)
#ax=plt.gca()
#w=0
#for v0 in vList:
	#u=theta[:,w]-theta0
	## We are following the particles, not the absolute coordinate. This is a Lagrangian strain tensor
	## u_s = 2 du/dR + (du/dR)^2
	#dudR=(u[1:]-u[:-1])/dphi
	#xval=theta[1:,w]-0.5*dphi*dudR
	#plt.plot(xval,strain[:,w],'.-',color=testmap(w), linestyle='solid',label=vList[w])
	#w+=1
#ax.set_xlabel(r'$\vartheta$')
#ax.set_ylabel(r'$u_s(\vartheta)$')
#plt.legend(loc=1,ncol=2)
#plt.title('Strain')

## Might as well go for the pressure
## p = -k*u_s, simply
## Well, except if there is a pre-stress:
## p = p0 -k u_s, p_0 is the initial stress on the particles
#print pini
#plt.figure(figsize=(10,7),linewidth=2.0)
#ax=plt.gca()
#w=0
#for v0 in vList:
	#u=theta[:,w]-theta0
	## We are following the particles, not the absolute coordinate. This is a Lagrangian strain tensor
	## u_s = 2 du/dR + (du/dR)^2
	#dudR=(u[1:]-u[:-1])/dphi
	#xval=theta[1:,w]-0.5*dphi*dudR
	#plt.plot(xval,pini-5*k*strain[:,w],'.-',color=testmap(w), linestyle='solid',label=vList[w])
	#w+=1
#ax.set_xlabel(r'$\vartheta$')
#ax.set_ylabel(r'$p(\vartheta)$')
#plt.legend(loc=1,ncol=2)
##plt.ylim(0,2)
##plt.xlim(-np.pi/2,np.pi/2)
#plt.title('Pressure - continuum and pre-stress')





# Going for broke: The simulation results, right here.
outfolder= '/home/silke/Documents/CurrentProjects/Rastko/analysis/data/SilkeJvar/'
outfolder_var= '/home/silke/Documents/CurrentProjects/Rastko/analysis/data/SilkeJ1/'

JList=['10','1','0.1']
#vList=['0.005','0.01','0.02','0.05','0.1','0.2','0.5','1']
vList=['0.1','0.5']
#JList=['10']
#vList=['1']
nbin=180
rval=28.2094791

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

for j in range(len(JList)):
	#plt.figure(figsize=(10,7),linewidth=2.0)
	for i in range(len(vList)):
		print vList[i],JList[j]
		#ax=plt.gca()
		if j==1:
			outfile=outfolder_var+'profilesV_v0' + vList[i] + '_j' + JList[j] + '.dat'
			outfile2=outfolder_var + 'axisV_v0' + vList[i] + '_j' + JList[j] + '.dat'
		else:
			outfile=outfolder+'profilesV_v0' + vList[i] + '_j' + JList[j] + '.dat'
			outfile2=outfolder + 'axisV_v0' + vList[i] + '_j' + JList[j] + '.dat'
		# header='theta rho vel energy pressure alpha alpha_v'
		profiles=np.loadtxt(outfile, unpack=True)[:,:] 
		isdata=[index for index,value in enumerate(profiles[1,:]) if (value >0)]
		# Corrected
		## Forgot the normalization of rho by the angle band width
		#if usecolumn==1:
			#normz=2*np.pi*rval*abs(np.cos(profiles[0,:]))
			#profiles[1,isdata]=profiles[1,isdata]/normz[isdata]
			#profiles[1,:]/=np.mean(profiles[1,:])
		if usecolumn==2:
			plt.plot(profiles[0,isdata],profiles[usecolumn,isdata]/float(vList[j]),color=testmap(i), linestyle='solid',label=JList[j])
		else:
			plt.plot(profiles[0,isdata],profiles[usecolumn,isdata],color=testmap(j), linestyle='solid',label=JList[j])
	#if usecolumn<=8:
		#plt.ylim(0,1.25*profiles[usecolumn,nbin/2])
	if usecolumn==9:
		plt.plot(profiles[0,isdata],profiles[0,isdata],'k-')
		plt.ylim(-0.5,0.5)
	plt.xlim(-np.pi/2,np.pi/2)
	#plt.xlabel(profList[0]) 
	#plt.ylabel(profList[usecolumn]) 
	plt.legend(loc=2,ncol=2)
	#plt.title('Pressure - simulation')
	

	
	
plt.show()