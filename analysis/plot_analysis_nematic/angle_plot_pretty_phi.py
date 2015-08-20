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
#vList=['0.2','0.3','0.5','0.7','1.0','1.5','2.0','3.0','5.0','7.0','10.0']
JList=['0.01','0.05','0.1','0.5','5.0','10.0']
#vList=['0.2','1.0','5.0']
#RList=['5.0','6.0','7.0','8.0','9.0','10.0','12.0','14.0','16.0','18.0','20.0','25.0','30.0','40.0']
phiList=['0.75','1.0']
R='16'

vmap=LinearSegmentedColormap('test',cdict,N=len(JList))
Rmap=LinearSegmentedColormap('test',cdict,N=len(phiList))

nsnap=2020
nskip=1000
nskip=500
step=5000
dt=0.001

r=0
Jval=np.zeros((len(JList),))
avtot=np.zeros((len(phiList),len(JList)))
davtot=np.zeros((len(phiList),len(JList)))
defectmat=np.empty((nsnap-nskip,4,3))
angles=np.zeros((nsnap-nskip,6))
avang=np.zeros((nsnap-nskip))
ang_correl=np.zeros((nsnap-nskip))
stdang=np.zeros((nsnap-nskip))
bins=np.linspace(0,180,45)
dbin=bins[1]-bins[0]
for phi in phiList:
	v=0
	plt.figure(figsize=(10,7),linewidth=2.0)
	for J in JList:	
		print phi,J
		#/home/silke/Documents/CurrentProjects/Rastko/nematic/data/phi_0.75/defects_phi0.75/defects_J_0.01_R_16.dat
		infile= basefolder +'/phi_' + phi+'/defects_phi' + phi + '/defects_J_' + J + '_R_' + R +'.dat'
		# header='theta rho vel energy pressure alpha alpha_v'
		datamat=(sp.loadtxt(infile, unpack=True)[:,(nskip+1):]).T 
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
		
		if len(isfour)>0:
			avtot[r,v]=np.mean(avang[isfour])
			# Mean of the variance: how far out ar the other angles?
			#davtot[r,v]=np.mean(stdang[isfour])
			# Variance of the mean: how much fluctuations are there?
			davtot[r,v]=np.std(avang[isfour])
		
		# Plotting as desired
		time=np.linspace(nskip*dt*step,nsnap*dt*step,nsnap-nskip)
		#plt.plot(time,avang,'-',color=vmap(v),label=JList[v])
		#plt.plot(time[isfour],avang[isfour],'-',color=vmap(v),label=JList[v])
		
		# Calculate correlation function - but only if were are on four angles throughout. 
		# Deal with the rest later ... complex situation
		# Same thing with the histogram - the whole flat vs. tetrahedron discussion is for those anyway
		if len(isfour)>(nsnap-nskip)/2.0:
			# Angular correlation
			# Nope, this misses the normalization ...
			avshift=avang-np.mean(avang)
			#ang_correl = np.correlate(avang-np.mean(avang),avang-np.mean(avang),mode='full')
			#ang_correl = ang_correl[:(nsnap-nskip)]
			#ang_correl /= ang_correl[0]
			for u in range(nsnap-nskip):
				ang_correl[u]=np.mean(avshift[0:(nsnap-nskip-u)]*avshift[u:])
			ang_correl /= ang_correl[0]
			time=np.linspace(0,(nsnap-nskip)*dt*step,nsnap-nskip)
			# Plotting for both (comment out as desired)
			plt.plot(time,ang_correl,'-',color=vmap(v),label=JList[v])
		# Histogram
		ang_hist,bin_edges =np.histogram(np.ravel(angles[isfour,:]),bins=bins,density=True)
		#plt.plot(bins[:44]+dbin/2,ang_hist,'.-',color=vmap(v),label=JList[v])
	
		Jval[v]=float(JList[v])
		v+=1
		## angles
		plt.xlabel('time') 
		plt.ylabel('angle') 
		## correlations
		#plt.xlabel('time') 
		#plt.ylabel('C(t)') 
		#plt.legend(loc=3,ncol=2)
		#plt.xlim(0,500)
		#plt.ylim(-0.5,1.05)
		# Histograms
		#plt.xlabel('angle') 
		#plt.ylabel('P(angle)') 
		plt.legend(loc=1,ncol=2)
		#plt.xlim(0,500)
		#plt.ylim(-0.5,1.05)
		
		plt.title('Density ' + str(phi))
	
	r+=1

#plt.figure(figsize=(10,7),linewidth=2.0)
#for r in range(len(RList)):
	#plt.errorbar(vval,avtot[r,:],yerr=davtot[r,:],color=Rmap(r),marker='o',label='R=' + RList[r])
#plt.plot(vval,109.47*vval/vval,'k--')
#plt.xlim(0,2)
#plt.ylim(0,140)
#plt.xlabel('v_0') 
#plt.ylabel('angle') 
#plt.legend(loc=3,ncol=1)
		
	
plt.show()


