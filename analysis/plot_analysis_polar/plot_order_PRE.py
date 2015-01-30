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
                   

	
basedir='/media/drogon/home/silke/Documents/Curved/Runs_rebuttal/'

	
# Order parameter n = |\frac{1}{N R v_0} \sum r \times v|
fig=plt.figure(figsize=(10,7),linewidth=2.0)
ax=plt.gca()
RList=['5','8','12','16','20','40','60']
#testmap2=LinearSegmentedColormap('test',cdict,N=len(RList))

JList=['1']
vList=['0.5']

orderpar=np.zeros((len(RList),))
dorder=np.zeros((len(RList),))
rval=np.zeros((len(RList),))
for i in range(len(vList)):
	for j in range(len(RList)):	
		print vList[i],RList[j]
		outfile2=basedir + '/axis_v0'  + vList[i] + '_R' + RList[j] + '.dat'
		axis=sp.loadtxt(outfile2, unpack=True)[:,:] 
		orderpar0=np.sqrt(axis[3,:]**2+axis[4,:]**2+axis[5,:]**2)
		orderpar[j]=np.mean(orderpar0)
		dorder[j]=np.std(orderpar0) 
		rval[j]=float(RList[j])
	plt.errorbar(rval,orderpar,yerr=dorder,color=(0.8,0,0), linestyle='solid',marker='s',markersize=10,)
	
plt.xlabel(r'$R$') 
plt.ylabel('p') 
#plt.ylim(0,1.1)
plt.xlim(0,65)
#plt.legend(loc=4,ncol=1)
plt.title('Order parameter')    

	
basedir='/media/drogon/home/silke/Documents/Curved/Runs_nuslices/'
# Order parameter n = |\frac{1}{N R v_0} \sum r \times v|
fig=plt.figure(figsize=(10,7),linewidth=2.0)
ax=plt.gca()
nuList=['0.001', '0.1','0.2', '0.3','0.5', '0.7','1','1.5', '2', '2.5', '3', '4', '5']


JList=['1']
vList=['0.01','0.1','0.5']
testmap2=LinearSegmentedColormap('test',cdict,N=len(vList))

orderpar=np.zeros((len(nuList),))
dorder=np.zeros((len(nuList),))
nuval=np.zeros((len(nuList),))
for i in range(len(vList)):
    for j in range(len(nuList)): 
        print vList[i],nuList[j]
        outfile2=basedir + '/axis_v0'  + vList[i] + '_nu' + nuList[j] + '.dat'
        axis=sp.loadtxt(outfile2, unpack=True)[:,:] 
        orderpar0=np.sqrt(axis[3,:]**2+axis[4,:]**2+axis[5,:]**2)
        orderpar[j]=np.mean(orderpar0)
        dorder[j]=np.std(orderpar0) 
        nuval[j]=float(nuList[j])
    plt.errorbar(nuval,orderpar,yerr=dorder,color=testmap2(i), linestyle='solid',marker='s',markersize=10,)
	


#hmm=np.linspace(-3,1,10)
#plt.plot(hmm,hmm/hmm,linestyle='--',color='k')	
	
plt.xlabel(r'$\nu$') 
plt.ylabel('p') 
#plt.ylim(0,1.1)
#plt.xlim(0,2.2)
plt.legend(loc=4,ncol=1)
plt.title('Order parameter')	
	
plt.show()	
		
		
		
		
		
		