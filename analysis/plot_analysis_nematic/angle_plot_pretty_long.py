# * *************************************************************
# *  
# *   Soft Active Mater on Surfaces (SAMoS)
# *   
# *   Author: Rastko Sknepnek
# *  
# *   Division of Physics
# *   School of Engineering, Physics and Mathematics
# *   University of Dundee
# *   
# *   (c) 2013, 2014
# * 
# *   School of Science and Engineering
# *   School of Life Sciences 
# *   University of Dundee
# * 
# *   (c) 2015
# * 
# *   Author: Silke Henkes
# * 
# *   Department of Physics 
# *   Institute for Complex Systems and Mathematical Biology
# *   University of Aberdeen  
# * 
# *   (c) 2014, 2015
# *  
# *   This program cannot be used, copied, or modified without
# *   explicit written permission of the authors.
# * 
# * ***************************************************************

import sys, os, glob
import cPickle as pickle
import numpy as np
import scipy as sp
#from scipy.io import savemat
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
#from matplotlib import rc
import matplotlib
#from mpl_toolkits.mplot3d import Axes3D
import sys
#import argparse
from read_data import *
#import numpy as np
import numpy.linalg as lin
#import matplotlib.pyplot as plt
import math as m
from datetime import *


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
                   

# This is the structured data file hierarchy. Replace as appropriate (do not go the Yaouen way and fully automatize ...)
basefolder= '/home/silke/Documents/CurrentProjects/Rastko/nematic/data/'
basefolder='/home/silke/Documents/CurrentProjects/Rastko/nematic/data/R_30.0_J0.1/'
#vList=['0.2','0.3','0.5','0.7','1.0','1.5','2.0','3.0','5.0','7.0','10.0']
#vList=['0.1','1.0','1.25','1.5','1.75','2.0','3.0']
#vList=['1.0','1.25','1.5','1.75','2.0','3.0']
vList=['0.2','0.3','0.4', '0.5','0.75']
#JList=['0.01','0.05']
#vList=['0.2','1.0','5.0']
#RList=['5.0','6.0','7.0','8.0','9.0','10.0','12.0','14.0','16.0','18.0','20.0','25.0','30.0','40.0']
phiList=['1.0']
R='30'
J='0.1'

vmap=LinearSegmentedColormap('test',cdict,N=len(vList))
Rmap=LinearSegmentedColormap('test',cdict,N=len(phiList))

#vList=['0.2','0.3','0.4']

#nsnap=20020
nskip=1000
nskip=500
step=5000
dt=0.001

r=0
#Jval=np.zeros((len(JList),))
#avtot=np.zeros((len(phiList),len(JList)))
#davtot=np.zeros((len(phiList),len(JList)))

bins=np.linspace(0,180,45)
dbin=bins[1]-bins[0]

#for J in JList:	
#plt.figure(figsize=(10,7),linewidth=2.0)
w=0
for v in vList:	
	#plt.figure(figsize=(10,7),linewidth=2.0)
	#print J
	print v
	#/home/silke/Documents/CurrentProjects/Rastko/nematic/data/R_16.0_long/defects_J_0.01_R_16_long.dat
	#infile= basefolder +'/R_' + R+'.0_long/defects_J_' + J + '_R_' + R +'_long.dat'
	infile= basefolder +'defects_J_' + J + 'v0_' + v +'_long.dat'
	print infile
	# header='theta rho vel energy pressure alpha alpha_v'
	datamat=(sp.loadtxt(infile, unpack=True)[:,(nskip+1):]).T 
	
	nsnap=len(datamat[:,0])+nskip
	defectmat=np.empty((nsnap-nskip,4,3))
	angles=np.zeros((nsnap-nskip,6))
	avang=np.zeros((nsnap-nskip))
	ang_correl=np.zeros((nsnap-nskip))
	stdang=np.zeros((nsnap-nskip))

	ndefect=datamat[:,0]
	
	# This includes potential zeros. Be very careful in the subsequent analysis
	# Currently divides by the absolute value of the first element. Which should be fine, but a bit imprecise ...
	defectmat[:,0,:] = datamat[:,1:4]/lin.norm(datamat[0,1:4])
	defectmat[:,1,:] = datamat[:,4:7]/lin.norm(datamat[0,4:7])
	defectmat[:,2,:] = datamat[:,7:10]/lin.norm(datamat[0,7:10])
	defectmat[:,3,:] = datamat[:,10:14]/lin.norm(datamat[0,10:14])
		
	# Defined for two defects
	angles[:,0] = np.degrees(np.arccos(np.sum(defectmat[:,0,:]*defectmat[:,1,:],axis=1)))
	# Defined for three defects
	angles[:,1] = np.degrees(np.arccos(np.sum(defectmat[:,0,:]*defectmat[:,2,:],axis=1)))
	angles[:,2] = np.degrees(np.arccos(np.sum(defectmat[:,1,:]*defectmat[:,2,:],axis=1)))
	# Defined for four defects
	angles[:,3] = np.degrees(np.arccos(np.sum(defectmat[:,0,:]*defectmat[:,3,:],axis=1)))
	angles[:,4] = np.degrees(np.arccos(np.sum(defectmat[:,1,:]*defectmat[:,3,:],axis=1)))
	angles[:,5] = np.degrees(np.arccos(np.sum(defectmat[:,2,:]*defectmat[:,3,:],axis=1)))
	
	huh=np.isnan(angles)
	angles[huh]=109
	# First stop: Correlations on the mean angle
	# The whole mean angle story only makes proper sense for four of them
	# Start with that, at least
	isfour = [index for index,value in enumerate(ndefect) if value >=4]
	isthree = [index for index,value in enumerate(ndefect.astype(int)) if value ==3]
	istwo = [index for index,value in enumerate(ndefect.astype(int)) if value ==2]
	# One defect is useless, physically impossible and I can't calculate angles anyway
	isnil = [index for index,value in enumerate(ndefect) if value <=1]
	print "System spends a fraction " + str(len(isfour)/(1.0*(nsnap-nskip))) +" with four defects, " + str(len(isthree)/(1.0*(nsnap-nskip))) + " with three defects, "
	print str(len(istwo)/(1.0*(nsnap-nskip))) + " with two defects and " + str(len(isnil)/(1.0*(nsnap-nskip))) + " with one defect or less."
	
	# only fill where defined
	avang[isfour]=np.mean(angles[isfour,:],axis=1)
	stdang[isfour]=np.std(angles[isfour,:],axis=1)
	avang[isthree]=np.mean(angles[isthree,0:2],axis=1)
	stdang[isthree]=np.std(angles[isthree,0:2],axis=1)
	avang[istwo]=angles[istwo,0]
	
	#if len(isfour)>0:
		#avtot[r,v]=np.mean(avang[isfour])
		## Mean of the variance: how far out ar the other angles?
		##davtot[r,v]=np.mean(stdang[isfour])
		## Variance of the mean: how much fluctuations are there?
		#davtot[r,v]=np.std(avang[isfour])
	
	# Plotting as desired
	time=np.linspace(nskip*dt*step,nsnap*dt*step,nsnap-nskip)
	#plt.plot(time,avang,'-',color=vmap(v),label=JList[v])
	#plt.plot(time[isfour],avang[isfour],'-',color=vmap(w),label=vList[w])
	
	# Calculate correlation function - but only if were are on four angles throughout. 
	# Deal with the rest later ... complex situation
	# Same thing with the histogram - the whole flat vs. tetrahedron discussion is for those anyway
	#if len(isfour)>(nsnap-nskip)/2.0:
	# Angular correlation
	# Nope, this misses the normalization ...
	avshift=avang-np.mean(avang)
	for u in range(len(avshift)):
		if np.isnan(avshift[u]):
			avshift[u]=0.0
	#ang_correl = np.correlate(avang-np.mean(avang),avang-np.mean(avang),mode='full')
	#ang_correl = ang_correl[:(nsnap-nskip)]
	#ang_correl /= ang_correl[0]
	for u in range(nsnap-nskip):
		ang_correl[u]=np.mean(avshift[0:(nsnap-nskip-u)]*avshift[u:])
	ang_correl /= ang_correl[0]
	time=np.linspace(0,(nsnap-nskip)*dt*step,nsnap-nskip)
	# Plotting for both (comment out as desired)
	plt.plot(time,ang_correl,'-',color=vmap(w),label=vList[w])
	# Histogram
	ang_hist,bin_edges =np.histogram(np.ravel(angles[isfour,:]),bins=bins,density=True)
	#plt.plot(bins[:44]+dbin/2,ang_hist,'.-',color='r')
	#plt.plot(bins[:44]+dbin/2,ang_hist,'.-',color=vmap(w),label=vList[w])

	
	## angles
	#plt.xlabel('time') 
	#plt.ylabel('angle') 
	#plt.legend(loc=3,ncol=2)
	## correlations
	plt.xlabel('time') 
	plt.ylabel('C(t)') 
	#plt.legend(loc=3,ncol=2)
	plt.xlim(0,50000)
	plt.ylim(-0.5,1.05)
	## Histograms
	plt.xlabel('angle') 
	plt.ylabel('P(angle)') 
	
	
	plt.title('Alignment ' + J)
	#plt.title('Velocity ' + str(v))
	w+=1

#plt.figure(figsize=(10,7),linewidth=2.0)
#for r in range(len(RList)):
	#plt.errorbar(vval,avtot[r,:],yerr=davtot[r,:],color=Rmap(r),marker='o',label='R=' + RList[r])
#plt.plot(vval,109.47*vval/vval,'k--')
#plt.xlim(0,2)
#plt.ylim(0,140)
#plt.xlabel('v_0') 
#plt.ylabel('angle') 
#plt.legend(loc=3,ncol=1)
plt.legend(loc=1,ncol=2)
	
plt.show()


