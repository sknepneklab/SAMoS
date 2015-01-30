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
def eng_extr(theta,k,v0):
    N = len(theta)
    two_sig = 2.0*np.ones(N-1)
    theta_diff = two_sig - R*(theta[1:]-theta[:-1])
    EngSpring = 0.5*k/R*sum(theta_diff**2)
    EngExt = -v0*sum(np.cos(theta))
    return (EngSpring+EngExt)

# Gradient part of the energies. That is equivalent to the equations of motion in the paper
def eng_grad(theta,k,v0):
    N = len(theta)
    der = np.zeros_like(theta)
    two_sig = 2.0*np.ones(N-1)
    theta_p1 = theta[2:]
    theta_m1 = theta[:-2]
    der[0] = k*(2.0 - R*(theta[1]-theta[0]))
    der[1:-1] = -k*R*(theta_p1 - 2*theta[1:-1] + theta_m1)
    der[-1] = -k*(2.0 - R*(theta[-1]-theta[-2]))
    der += v0*np.sin(theta)
    return der
   
# Produce results as a function of v_0. I am using simulation parameters here
R = 28.2094791774
k = 1
sigma=1
# Initial state: just touching parameters ...
N=44
# real N for dphi= 2 \sigma/R (just touching initial condition)
#N=int(round(R*np.pi/(2*sigma)))
print N
#vList=[0.005,0.01,0.02,0.05,0.1,0.2,0.5,1.0]
vList=[0.005,0.01,0.02,0.05,0.1,0.2,0.5]
testmap=LinearSegmentedColormap('test',cdict,N=len(vList))
# Initial condition: String of pearls from -pi/2 to pi/2. Note that dphi is *not* 2 \sigma /R.
dphi = math.pi/float(N)
#theta0 = np.array([-0.5*math.pi+(i+0.5)*dphi for i in xrange(N)])
theta0=np.linspace(-np.pi/2,np.pi/2,N)
#print theta0
bounds = [(-np.pi/2, np.pi/2)] * N
theta=np.zeros((N,len(vList)))
strain=np.zeros((N-1,len(vList)))
w=0
plt.figure(figsize=(10,7),linewidth=2.0)
ax=plt.gca()
for v0 in vList:
	# http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html
	# Minimize the energy associated to these theta positions, using BFGS. 
	# jac=eng_grad associates our gradient function to the jacobian, instead of a numerical estimate
	# The output is:
	# res.x: The solution array, that is the new positions of the particles, I assume?
	# For the rest: see here
	# http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.OptimizeResult.html#scipy.optimize.OptimizeResult
	#res = minimize(eng_extr, theta0, args=(k,v0), method='BFGS', jac=eng_grad)
	res = minimize(eng_extr, theta0, args=(k,v0), method='L-BFGS-B', jac=eng_grad,bounds=bounds)
	theta[:,w]=res.x 
	# Displacement field
	u=res.x-theta0
	# We are following the particles, not the absolute coordinate. This is a Lagrangian strain tensor
	# u_s = 2 du/dR + (du/dR)^2
	#dudR=(res.x[1:]-res.x[:-1])/dphi
	# or equivalently
	dudR=(u[1:]-u[:-1])/dphi
	strain0=2*dudR+dudR*dudR
	strain[:,w]=strain0
	#plt.plot(theta[:,w],strain[:,w],'.-',color=testmap(w), linestyle='solid',label=vList[w])
	plt.plot(theta[:,w],u,'.-',color=testmap(w), linestyle='solid',label=vList[w])
	w+=1
ax.set_xlabel(r'$\vartheta$')
ax.set_ylabel(r'$u(\vartheta)$')
plt.legend(loc=1,ncol=2)
plt.title('Displacements')


plt.figure(figsize=(10,7),linewidth=2.0)
ax=plt.gca()
w=0
for v0 in vList:
	u=theta[:,w]-theta0
	# We are following the particles, not the absolute coordinate. This is a Lagrangian strain tensor
	# u_s = 2 du/dR + (du/dR)^2
	#dudR=(res.x[1:]-res.x[:-1])/dphi
	# or equivalently
	dudR=(u[1:]-u[:-1])/dphi
	xval=theta[1:,w]-0.5*dphi*dudR
	plt.plot(xval,strain[:,w],'.-',color=testmap(w), linestyle='solid',label=vList[w])
	w+=1
ax.set_xlabel(r'$\vartheta$')
ax.set_ylabel(r'$u_s(\vartheta)$')
plt.legend(loc=1,ncol=2)
plt.title('Strain')


# Might as well go for the pressure
# p = -k*u_s, simply
# Well, except if there is a pre-stress:
# p = p0 -k u_s, p_0 is the initial stress on the particles
pini=2*sigma*np.sum(2.0*sigma - R*(theta0[1:]-theta0[:-1]))/(np.pi*R)
#pini=0
print pini
plt.figure(figsize=(10,7),linewidth=2.0)
ax=plt.gca()
w=0
for v0 in vList:
	u=theta[:,w]-theta0
	# We are following the particles, not the absolute coordinate. This is a Lagrangian strain tensor
	# u_s = 2 du/dR + (du/dR)^2
	#dudR=(res.x[1:]-res.x[:-1])/dphi
	# or equivalently
	dudR=(u[1:]-u[:-1])/dphi
	xval=theta[1:,w]-0.5*dphi*dudR
	plt.plot(xval,pini-k*strain[:,w],'.-',color=testmap(w), linestyle='solid',label=vList[w])
	w+=1
ax.set_xlabel(r'$\vartheta$')
ax.set_ylabel(r'$p(\vartheta)$')
plt.legend(loc=1,ncol=2)
plt.ylim(0,1.5)
#plt.xlim(-np.pi/2,np.pi/2)
plt.title('Pressure')

## DOH.
## Might as well go for the pressure
## p = -k*u_s, simply
## Well, except if there is a pre-stress:
## p = p0 -k u_s, p_0 is the initial stress on the particles
#pini=2*sigma*np.sum(2.0*sigma - R*(theta0[1:]-theta0[:-1]))/(np.pi*R)
#print pini
#plt.figure(figsize=(10,7),linewidth=2.0)
#ax=plt.gca()
#w=0
#for v0 in vList:
	#u=theta[:,w]-theta0
	## We are following the particles, not the absolute coordinate. This is a Lagrangian strain tensor
	## u_s = 2 du/dR + (du/dR)^2
	##dudR=(res.x[1:]-res.x[:-1])/dphi
	## or equivalently
	#dudR=(u[1:]-u[:-1])/dphi
	#xval=theta[1:,w]-0.5*dphi*dudR
	#press=2*sigma*(2.0*sigma - R*(theta[1:,w]-theta[:-1,w]))/(np.pi*R)
	#plt.plot(xval,press,'.-',color=testmap(w), linestyle='solid',label=vList[w])
	#w+=1
#ax.set_xlabel(r'$\vartheta$')
#ax.set_ylabel(r'$p(\vartheta)$')
#plt.legend(loc=1,ncol=2)
##plt.ylim(0,2)
##plt.xlim(-np.pi/2,np.pi/2)
#plt.title('Pressure')

# Analytical results, as far as can be done by hand
# This is the Lagrangian result, to 0th order (the only one doable by hand)
# It has the same generic form, with a solution for theta, and then a strain
plt.figure(figsize=(10,7),linewidth=2.0)
ax=plt.gca()
w=0
for v0 in vList:
	# alpha = R/(2*sigma)^2 v_0/(mu k)
	# beta = v_0/(2 \sigma mu k)
	alpha=R/(2*sigma)**2*v0/k
	print alpha
	beta=v0/(2*sigma*k)
	print beta
	u=-alpha*np.sin(theta0)-beta*theta0
	us=-2*alpha*np.cos(theta0)-2*beta+alpha**2*np.cos(theta0)**2+2*alpha*beta*np.cos(theta0)+beta**2
	plt.plot(theta0+u,-k*us,'-',color=testmap(w), linestyle='solid',label=vList[w])
	#plt.plot(theta0+u,u,'-',color=testmap(w), linestyle='solid',label=vList[w])
	w+=1
ax.set_xlabel(r'$\vartheta$')
ax.set_ylabel(r'$p(\vartheta)$')
plt.ylim(0,1.5)
plt.xlim(-np.pi/2,np.pi/2)
plt.legend(loc=1,ncol=2)
plt.title('Pressure, Lagrange 0th order')

# Eulerian result, 0th order
# From Mathematica: the numerical results for thetaM, for our exact parameters, as a function of v_0
thetaMList=[1.53373, 1.50041, 1.44261, 1.31511, 1.18064, 1.02631, 0.818654,0.673751]
plt.figure(figsize=(10,7),linewidth=2.0)
ax=plt.gca()
w=0
for v0 in vList:
	# alpha = R/(2*sigma)^2 v_0/(mu k)
	# beta = v_0/(2 \sigma mu k)
	alpha=R/(2*sigma)**2*v0/k
	print alpha
	beta=v0/(2*sigma*k)
	print beta
	tm=thetaMList[w]
	theta=np.linspace(-tm,tm,N)
	c=alpha*np.cos(tm)-beta*np.sin(tm)
	u=-alpha*np.sin(theta)+c*theta
	us=-2*alpha*np.cos(theta)+2*c-alpha**2*np.cos(theta)**2+2*alpha*c*np.cos(theta)-c**2
	plt.plot(theta,-k*us,'-',color=testmap(w), linestyle='solid',label=vList[w])
	#plt.plot(theta,u,'-',color=testmap(w), linestyle='solid',label=vList[w])
	w+=1
ax.set_xlabel(r'$\vartheta$')
ax.set_ylabel(r'$p(\vartheta)$')
plt.ylim(0,1.5)
plt.xlim(-np.pi/2,np.pi/2)
plt.legend(loc=1,ncol=2)
plt.title('Pressure, Euler 0th order')


plt.show()
