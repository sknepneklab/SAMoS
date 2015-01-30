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
basefolder = '/home/silke/Documents/CurrentProjects/Rastko/Runs/'
outfolder= '/home/silke/Documents/CurrentProjects/Rastko/analysis/data/Rastko/'
outpics='/home/silke/Documents/CurrentProjects/Rastko/analysis/pics/Rastko/'
# Slice 1
#JList=['0.1','1.0']
#vList=['0.05','0.5','5.0']
# Slice 2
JList=['0.01','0.10','1.00','10.00']
vList=['0.01','0.03','0.1','0.3','1.0']
nbin=180
rval=28.2094791
testmap=LinearSegmentedColormap('test',cdict,N=len(vList))
testmap2=LinearSegmentedColormap('test',cdict,N=len(JList))
nstep=20200000
nsave=10000
nsnap=int(nstep/nsave)
skip=500
dt=0.001


# Profiles
# Set column to plot
usecolumn=4

profList=[r'$\theta$',r'$\rho$',r'$\sqrt{\langle v^2 \rangle}/v_0$','energy','pressure',r'$\Sigma_{\theta \theta}$',r'$\Sigma_{\theta \phi}$',r'$\Sigma_{\phi \theta}$',r'$\Sigma_{\phi \phi}$',r'$\alpha$',r'$\alpha_v$']
profName=['theta','rho','vrms','energy','pressure','stt','stp','spt','spp','alpha','alpha_v']

for j in range(len(JList)):
	plt.figure(figsize=(10,7),linewidth=2.0)
	for i in range(len(vList)):
		print vList[i],JList[j]
		ax=plt.gca()
		outfile=outfolder+'profilesV_v0' + vList[i] + '_j' + JList[j] + '.dat'
		outfile2=outfolder + 'axisV_v0' + vList[i] + '_j' + JList[j] + '.dat'
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
			plt.plot(profiles[0,isdata],profiles[usecolumn,isdata]/float(vList[i]),color=testmap(i), linestyle='solid',label=vList[i])
		else:
			plt.plot(profiles[0,isdata],profiles[usecolumn,isdata],color=testmap(i), linestyle='solid',label=vList[i])
	#if usecolumn<=8:
		#plt.ylim(0,1.25*profiles[usecolumn,nbin/2])
	if usecolumn==9:
		plt.plot(profiles[0,isdata],profiles[0,isdata],'k--')
		#plt.ylim(-0.5,0.5)
	plt.xlim(-np.pi/2,np.pi/2)
	plt.xlabel(profList[0]) 
	plt.ylabel(profList[usecolumn]) 
	plt.legend(loc=2,ncol=2)
	plt.title('Interaction strength ' + r'$J=$' + JList[j])
	
	filename=outpics + '/profile_' + profName[usecolumn] + '_J' + JList[j] +'_slice2.pdf'
	plt.savefig(filename)
	
	
#for i in range(len(vList)):	
	#plt.figure(figsize=(10,7),linewidth=2.0)
	#for j in range(len(JList)):
		#print vList[i],JList[j]
		#ax=plt.gca()
		#outfile=outfolder+'/profilesV_v0' + vList[i] + '_j' + JList[j] + '.dat'
		#outfile2=outfolder + '/axisV_v0' + vList[i] + '_j' + JList[j] + '.dat'
		## header='theta rho vel energy pressure alpha alpha_v'
		#profiles=sp.loadtxt(outfile, unpack=True)[:,:] 
		#isdata=[index for index,value in enumerate(profiles[1,:]) if (value >0)]
		#if j==0:
			#plt.plot(profiles[0,isdata],profiles[4,isdata],color='k', linestyle='solid',label=r'$Tr(\Sigma)$')
			#plt.plot(profiles[0,isdata],profiles[5,isdata],color='r', linestyle='solid',label=r'$\Sigma_{\theta \theta}$')
			#plt.plot(profiles[0,isdata],profiles[6,isdata],color='g', linestyle='solid',label=r'$\Sigma_{\theta \phi}$')
			#plt.plot(profiles[0,isdata],profiles[8,isdata],color='b', linestyle='solid',label=r'$\Sigma_{\phi \phi}$')
		#else:
			#plt.plot(profiles[0,isdata],profiles[4,isdata],color='k', linestyle='solid')
			#plt.plot(profiles[0,isdata],profiles[5,isdata],color='r', linestyle='solid')
			#plt.plot(profiles[0,isdata],profiles[6,isdata],color='g', linestyle='solid')
			#plt.plot(profiles[0,isdata],profiles[8,isdata],color='b', linestyle='solid')
	##plt.xlim(-np.pi/2,np.pi/2)
	#plt.xlabel(profList[0]) 
	#plt.ylabel('stresses') 
	#plt.legend(loc=2,ncol=2)
	#plt.title('Velocity ' + r'$v_0=$' + vList[i])
	
	##filename=outfolder + 'stress_profiles_v' + vList[i] +'.pdf'
	##plt.savefig(filename)
	
## Axis in 3d
#fig=plt.figure(figsize=(10,7),linewidth=2.0)
#ax = fig.add_subplot(111, projection='3d')
#for i in range(len(vList)):	
	#for j in range(len(JList)):
		#print vList[i],JList[j]
		#ax=plt.gca()
		#outfile=outfolder+'data_Rastko/profilesV_v0' + vList[i] + '_j' + JList[j] + '.dat'
		#outfile2=outfolder + 'data_Rastko/axisV_v0' + vList[i] + '_j' + JList[j] + '.dat'
		#axis=sp.loadtxt(outfile2, unpack=True)[:,:] 
		#ax.scatter(axis[0,:], axis[1,:], axis[2,:], zdir='z',color=testmap(i), linestyle='solid',label=vList[i])
#plt.xlabel('x') 
#plt.ylabel('y') 
##plt.xlim(-1,1)
##plt.ylim(-1,1)
##plt.zlim(-1,1)
##ax.zlabel('z') 
#ax.legend()
#plt.legend()
#plt.title('Velocity ' + r'$v_0=$' + vList[i])
	
	#filename=outfolder + 'pics/axisV_' + profName[usecolumn] + '_v' + vList[i] +'.pdf'
	#plt.savefig(filename)
	
## Order parameter n = |\frac{1}{N R v_0} \sum r \times v|
#fig=plt.figure(figsize=(10,7),linewidth=2.0)
#ax=plt.gca()
#for i in range(len(vList)):	
	#for j in range(len(JList)):
		#print vList[i],JList[j]	
		#outfile=outfolder+'data_Rastko/profilesV_v0' + vList[i] + '_j' + JList[j] + '.dat'
		#outfile2=outfolder + 'data_Rastko/axisV_v0' + vList[i] + '_j' + JList[j] + '.dat'
		#axis=sp.loadtxt(outfile2, unpack=True)[:,:] 
		#orderpar=np.sqrt(axis[3,:]**2+axis[4,:]**2+axis[5,:]**2)
		#xval=np.linspace(skip*nsave*dt,nstep*dt,(nsnap-skip+2))
		#plt.plot(xval,orderpar,color=testmap(i), linestyle='solid',label=vList[i])
#plt.xlabel('x') 
#plt.ylabel('y') 
#plt.ylim(0,1)
#plt.legend(loc=3,ncol=2)
	
	#filename=outfolder + 'pics/orderpar_v' + vList[i] +'.pdf'
	#plt.savefig(filename)
	
## Order parameter n = |\frac{1}{N R v_0} \sum r \times v|
#fig=plt.figure(figsize=(10,7),linewidth=2.0)
#ax=plt.gca()
#orderpar=np.zeros((len(vList),))
#dorder=np.zeros((len(vList),))
#vval=np.zeros((len(vList),))
#for j in range(len(JList)):
	#for i in range(len(vList)):	
		#print vList[i],JList[j]
		#outfile=outfolder+'profilesV_v0' + vList[i] + '_j' + JList[j] + '.dat'
		#outfile2=outfolder + 'axisV_v0' + vList[i] + '_j' + JList[j] + '.dat'
		#axis=sp.loadtxt(outfile2, unpack=True)[:,:] 
		#orderpar0=np.sqrt(axis[3,:]**2+axis[4,:]**2+axis[5,:]**2)
		#orderpar[i]=np.mean(orderpar0)
		#dorder[i]=np.std(orderpar0) 
		#vval[i]=np.log10(float(vList[i]))
	#plt.errorbar(vval,orderpar,yerr=dorder,color=testmap2(j), linestyle='solid',label='J='+JList[j])
#plt.xlabel(r'$\log_{10}v_0$') 
#plt.ylabel('n') 
#plt.ylim(0,1)
#plt.legend(loc=2,ncol=1)
#plt.title('Order parameter')
	
	#filename=outfolder + 'pics/orderpar_v' + vList[i] +'.pdf'
	#plt.savefig(filename)
	


	
	
plt.show()	
		
		
		
		
		
		