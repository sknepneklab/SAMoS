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
                   

	
baseJvarR='/home/silke/Documents/CurrentProjects/Rastko/analysis/data/Rastko/'
baseJvarS='/home/silke/Documents/CurrentProjects/Rastko/analysis/data/SilkeJvar/'
baseJ1='/home/silke/Documents/CurrentProjects/Rastko/analysis/data/SilkeJ1/'
baseFlat='/home/silke/Documents/CurrentProjects/Rastko/analysis/data/SilkefFlat/'



## Order parameter n = |\frac{1}{N R v_0} \sum r \times v|
#for i in range(len(vList)):	
	#fig=plt.figure(figsize=(10,7),linewidth=2.0)
	#for j in range(len(JList)):
		#print vList[i],JList[j]
		#ax=plt.gca()
		#outfile2=baseJ1 + '/axisV_v0' + vList[i] + '_j' + JList[j] + '.dat'
		#axis=sp.loadtxt(outfile2, unpack=True)[:,:] 
		#orderpar=np.sqrt(axis[3,:]**2+axis[4,:]**2+axis[5,:]**2)
		#xval=np.linspace(skip*nsave*dt,nstep*dt,(nsnap-skip+2))
		#plt.plot(xval,orderpar,color=testmap2(j), linestyle='solid',label=JList[j])
	#plt.xlabel('x') 
	#plt.ylabel('y') 
	#plt.ylim(0,1)
	#plt.legend(loc=3,ncol=2)
	#plt.title('Velocity ' + r'$v_0=$' + vList[i])
	
# Order parameter n = |\frac{1}{N R v_0} \sum r \times v|
fig=plt.figure(figsize=(10,7),linewidth=2.0)
ax=plt.gca()
JList=['10', '1', '0.1', '0.01']
testmap2=LinearSegmentedColormap('test',cdict,N=len(JList))

JList=['0.1']
#vList=['0.05','0.5','5.0']
#vList=['0.01','0.05','0.1','0.2','0.5','1.0']
#vList=['0.05','0.5','5.0']
vList=['0.01','0.05','0.1','0.2','0.5','1.0','5.0']


orderpar=np.zeros((len(vList),))
dorder=np.zeros((len(vList),))
vval=np.zeros((len(vList),))
for j in range(len(JList)):
	for i in range(len(vList)):	
		print vList[i],JList[j]
		if vList[i]=='5.0':
			outfile2=baseJvarR + '/axisV_v0' + vList[i] + '_j' + JList[j] + '.dat'
		else:
			outfile2=baseJvarS + '/axisV_v0' + vList[i] + '_j' + JList[j] + '.dat'
		axis=sp.loadtxt(outfile2, unpack=True)[:,:] 
		orderpar0=np.sqrt(axis[3,:]**2+axis[4,:]**2+axis[5,:]**2)
		orderpar[i]=np.mean(orderpar0)
		dorder[i]=np.std(orderpar0) 
		vval[i]=np.log10(float(vList[i]))
	plt.errorbar(vval,orderpar,yerr=dorder,color=testmap2(1), linestyle='solid',marker='s',markersize=10,label='J='+JList[j])
	
JList=['10', '1', '0.1', '0.01']
vList=['0.005','0.01','0.02','0.05','0.1','0.2','0.5','1']
#vList=['0.005','0.05','0.1','0.2','0.5','1']
#JList=['10']
#vList=['1']

testmap2=LinearSegmentedColormap('test',cdict,N=len(JList))
nstep=10000000
nsave=10000
nsnap=int(nstep/nsave)
skip=int(nsnap/2)
#skip=0
dt=0.001


orderpar=np.zeros((len(vList),len(JList)))
order=np.zeros((len(vList),))
dorder=np.zeros((len(vList),))
vval=np.zeros((len(vList),))

for i in range(len(vList)):	
	for j in range(len(JList)):
		print vList[i],JList[j]
		outfile2=baseJ1 + '/axisV_v0' + vList[i] + '_j' + JList[j] + '.dat'
		axis=sp.loadtxt(outfile2, unpack=True)[:,:] 
		orderpar0=np.sqrt(axis[3,:]**2+axis[4,:]**2+axis[5,:]**2)
		orderpar[i,j]=np.mean(orderpar0)
	order[i]=np.mean(orderpar[i,:])
	dorder[i]=np.std(orderpar[i,:])/np.sqrt(len(JList))
	vval[i]=np.log10(float(vList[i]))
plt.errorbar(vval,order,yerr=dorder,color=testmap2(2), linestyle='solid',marker='s',markersize=10,label='J=1')

#orderpar=np.zeros((len(vList),))
#dorder=np.zeros((len(vList),))
#vval=np.zeros((len(vList),))
#for j in range(len(JList)):
	#for i in range(len(vList)):	
		#print vList[i],JList[j]
		#outfile2=baseJ1 + '/axisV_v0' + vList[i] + '_j' + JList[j] + '.dat'
		#axis=sp.loadtxt(outfile2, unpack=True)[:,:] 
		#orderpar0=np.sqrt(axis[3,:]**2+axis[4,:]**2+axis[5,:]**2)
		#orderpar[i]=np.mean(orderpar0)
		#dorder[i]=np.std(orderpar0) 
		#vval[i]=np.log10(float(vList[i]))
	#plt.errorbar(vval,orderpar,yerr=dorder,color=testmap2(j), linestyle='solid',marker='s',markersize=10,label='J='+JList[j])





nstep=20000000
nsave=10000
nsnap=int(nstep/nsave)
skip=int(nsnap/3)
vList=['0.01','0.05','0.1','0.2','0.5','1.0']
JList=['0.1','10']

orderpar=np.zeros((len(vList),))
dorder=np.zeros((len(vList),))
vval=np.zeros((len(vList),))
for j in range(1,len(JList)):
	for i in range(len(vList)):	
		print vList[i],JList[j]
		outfile2=baseJvarS + '/axisV_v0' + vList[i] + '_j' + JList[j] + '.dat'
		axis=sp.loadtxt(outfile2, unpack=True)[:,:] 
		orderpar0=np.sqrt(axis[3,:]**2+axis[4,:]**2+axis[5,:]**2)
		orderpar[i]=np.mean(orderpar0)
		dorder[i]=np.std(orderpar0) 
		vval[i]=np.log10(float(vList[i]))
	plt.errorbar(vval,orderpar,yerr=dorder,color=testmap2(3), linestyle='solid',marker='s',markersize=10,label='J='+JList[j])

#JList=['0.01','0.1','1', '10']
JList=['1']
vList=['0.005','0.01','0.02','0.05','0.1','0.2','0.5','1']
testmap=LinearSegmentedColormap('test',cdict,N=len(vList))
#testmap2=LinearSegmentedColormap('test',cdict,N=len(JList))
nstep=10000000
nsave=10000
nsnap=int(nstep/nsave)
skip=int(nsnap/2)
dt=0.001

orderpar=np.zeros((len(vList),))
dorder=np.zeros((len(vList),))
vval=np.zeros((len(vList),))
for j in range(0,len(JList)):
	for i in range(len(vList)):	
		print vList[i],JList[j]
		outfile2=baseFlat + '/flatstat_v0' + vList[i] + '_j' + JList[j] + '.dat'
		axis=sp.loadtxt(outfile2, unpack=True)[:,:] 
		orderpar0=np.sqrt(axis[17,skip:]**2+axis[18,skip:]**2+axis[19,skip:]**2)/float(vList[i])	
		orderpar[i]=np.mean(orderpar0)
		dorder[i]=np.std(orderpar0) 
		vval[i]=np.log10(float(vList[i]))
	plt.errorbar(vval,orderpar,yerr=dorder,color=testmap2(2), linestyle='solid',marker='o',markersize=10,label='flat, J='+JList[j])
hmm=np.linspace(-3,1,10)
plt.plot(hmm,hmm/hmm,linestyle='--',color='k')	
	
plt.xlabel(r'$\log v_0$') 
plt.ylabel('p,n') 
plt.ylim(0,1.1)
plt.xlim(-2.5,1)
plt.legend(loc=4,ncol=1)
#plt.title('Order parameter')	
	
plt.show()	
		
		
		
		
		
		