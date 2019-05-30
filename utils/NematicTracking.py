# ***************************************************************************
# *
# *  Copyright (C) 2013-2016 University of Dundee
# *  All rights reserved. 
# *
# *  This file is part of SAMoS (Soft Active Matter on Surfaces) program.
# *
# *  SAMoS is free software; you can redistribute it and/or modify
# *  it under the terms of the GNU General Public License as published by
# *  the Free Software Foundation; either version 2 of the License, or
# *  (at your option) any later version.
# *
# *  SAMoS is distributed in the hope that it will be useful,
# *  but WITHOUT ANY WARRANTY; without even the implied warranty of
# *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# *  GNU General Public License for more details.
# *
# *  You should have received a copy of the GNU General Public License
# *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
# *
# *****************************************************************************


import sys
sys.path.insert(0,  '/home/cpsmth/s01sh3/Documents/SAMoS/samos/utils/')

import argparse
import pickle

import numpy as np
from Geometry import *
#from read_param import *

import matplotlib
import matplotlib.pyplot as plt
#import matplotlib.patches as ptch
#import matplotlib.lines as lne
from matplotlib.colors import LinearSegmentedColormap

#import matplotlib
#matplotlib.rcParams['text.usetex'] = 'true'
matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['axes.linewidth'] = 2
matplotlib.rcParams['xtick.major.size'] = 8
matplotlib.rcParams['ytick.major.size'] = 8
matplotlib.rcParams['font.size']=16.0
matplotlib.rcParams['legend.fontsize']=14.0

matplotlib.rcParams['text.usetex'] = 'false'
matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['axes.linewidth'] = 2
matplotlib.rcParams['xtick.major.size'] = 8
matplotlib.rcParams['ytick.major.size'] = 8
matplotlib.rcParams['font.size']=16.0
matplotlib.rcParams['legend.fontsize']=14.0

cdict = {'red':   [(0.0,  0.5, 1.0),
		   		   (0.3,  1.0, 0.5),
                   (0.5,  0.5, 0.0),
                   (1.0,  0.0, 0.0)],

         'green': [(0.0,  0.0, 0.2),
				   (0.35,  0.0, 0.5),
                   (0.5, 1.0, 1.0),
                   (0.8,  0.5, 0.0),
                   (1.0,  0.0, 0.0)],

         'blue':  [(0.0,  0.0, 0.0),
                   (0.5,  0.0, 0.0),
                   (0.7, 0.5, 0.75),
                   (1.0,  1.0, 0.5)]}
         
# Hack to avoid the whole parameter read-in saga
# and use the geometry class anyways
class param:
    def __init__(self,R):
        self.r=R
        self.box=[100, 100, 100]
        self.periodic=False

# Set colors for the 4 defect tracks in the most common case
defectsmap=LinearSegmentedColormap('test',cdict,N=4) 
		
# This is very little data. No problem having a little bit of redundancy in analysis here ...
# defect tracking

class TrackDefects:
    
    def __init__(self,filename,R):
	data=pickle.load(open(filename))
	self.ndefects=data["numdefects_n"]
	self.defects=data["defects_n"]
	par=param(R)
        # Create geometry (coordinates) for this sphere
        self.geom=GeometrySphere(par)
	
    def sortDefects(self,verbose=False):	
	# Get the Euler angles for all the defects, and then plot, according to charge
	# These are not tracked as such, hence the 0 label
	rvalhalf0=[]
	rvalone0=[]
	# Temporal loop over defect list
	nsnaps=1.0*len(self.defects1)
	self.track2=[]
	self.track4=[]
	self.isBand=False
	self.isFour=False
	# purely for output reasons ...
	self.rvalhalf=0
	self.rvalone=0
	for u in range(0,len(self.defects1)):
	    # Dig out defects, associate position, charge and angles
	    ndef=len(self.defects1[u])
	    charge=[]
	    angles=np.zeros((ndef,2))
	    rval=np.zeros((ndef,3))
	    for n in range(ndef):
		charge.append(self.defects1[u][n][0])
		rval[n,:]=self.defects1[u][n][1:4]
	    # Sort them into defects of charge +1/2, -1/2 and 1.0
	    ishalf=[index for index,value in enumerate(charge) if value==0.5]
	    ismhalf=[index for index,value in enumerate(charge) if value==-0.5]
	    isone=[index for index,value in enumerate(charge) if value==1.0]
	    if verbose:
		print len(ishalf)
		print len(ismhalf)
		print len(isone)
	    # See if anything is left over:
	    ntypical=len(ishalf)+len(ismhalf)+len(isone)
	    if (ndef-ntypical)>0:
		print "There are " + str(ndef-ntypical) + " unusual defects in snapshot " + str(u)
		print charge
	    # Check if we are in one of two common scenarious: a) 4 +1/2 defects; or b) 2 +1 defects
	    if len(ishalf)==4:
		rvalhalf0.append(rval[ishalf,:])
		self.track4.append(u)
	    if len(isone)==2:
		rvalone0.append(rval[isone,:])
		self.track2.append(u)
	    if verbose:
		if (len(ishalf)!=4) and (len(isone)!=2):
		    print u
		    print "# of 1/2 defects: " + str(len(ishalf))
		    print "# of -1/2 defects: " + str(len(ismhalf))
		    print "# of 1 defects: " + str(len(isone))
	
	# Sort out the tracks if possible
	# If we have a substantial set of snapshots where there were 4 +1/2 defects
	if len(self.track4)>0.25*len(self.defects1):
	    #trackFourDefects(rvalhalf0,track4,ax,polar)
	    self.isFour=True
	    print "Detected a predominantly four defect structure! Tracking!"
	    self.trackFourDefects(rvalhalf0)
	# If we have a substantial set of snapshots where there are 2 +1 defects
	if len(self.track2)>0.25*len(self.defects1):
	    #trackTwoDefects(rvalone0,track2,ax,polar)
	    self.isBand=True
	    print "Detected a predominantly band structure! Tracking!"
	    self.trackTwoDefects(rvalone0)
	# In case this is not true, but we want to still attempt to really, really track stuff
	# This will do the whole hog and go for all identifiable defect traces
	#if ForceTrack:
	
	
    def plotDefects(self,title,merged,polar,verbose=False):
	# Initialize a polar or angular plot
	onefig=plt.figure()
	if polar:
	    ax = plt.subplot(111, polar=True)
	else:
	    ax1 = plt.subplot(211)
	    ax2 = plt.subplot(212)
	    
	if merged:
	    nsnaps=len(self.defects1)
	else:
	    nsnaps=len(self.defects)
	for u in range(0,nsnaps):
	    # Dig out defects, associate position, charge and angles
	    if merged:
		ndef=len(self.defects1[u])
	    else:
		ndef=len(self.defects[u])
	    charge=[]
	    angles=np.zeros((ndef,2))
	    rval=np.zeros((ndef,3))
	    if merged:
		for n in range(ndef):
		    charge.append(self.defects1[u][n][0])
		    rval[n,:]=self.defects1[u][n][1:4]
	    else:
		for n in range(ndef):
		    charge.append(self.defects[u][n][0])
		    rval[n,:]=self.defects[u][n][1:4]
	    # Sort them into defects of charge +1/2, -1/2 and 1.0
	    ishalf=[index for index,value in enumerate(charge) if value==0.5]
	    ismhalf=[index for index,value in enumerate(charge) if value==-0.5]
	    isone=[index for index,value in enumerate(charge) if value==1.0]    
	    # Get the angles associated to the defects
	    angles[:,0],angles[:,1],etheta,ephi=self.geom.TangentBundle(rval)  
	    if polar:
		ax.plot(angles[ishalf,1],angles[ishalf,0], color='k',linestyle='none',marker='.')
		ax.plot(angles[ismhalf,1],angles[ismhalf,0], color='r',linestyle='none',marker='.')
		ax.plot(angles[isone,1],angles[isone,0], color='g',linestyle='none',marker='.')
	    else:
		# first angle, then second angle
		pu=[u for k in ishalf]
		ax1.plot(pu,angles[ishalf,0], color='k',linestyle='none',marker='.')
		ax2.plot(pu,angles[ishalf,1], color='k',linestyle='none',marker='.')
		pu=[u for k in ismhalf]
		ax1.plot(pu,angles[ismhalf,0], color='r',linestyle='none',marker='.')
		ax2.plot(pu,angles[ismhalf,1], color='r',linestyle='none',marker='.')
		pu=[u for k in isone]
		ax1.plot(pu,angles[isone,0], color='g',linestyle='none',marker='.')
		ax2.plot(pu,angles[isone,1], color='g',linestyle='none',marker='.')
	# Add tracks in case they exist 
	if self.isFour:	
	    angles=np.zeros((len(self.track4),2))
	    for a in range(4):
		angles[:,0],angles[:,1],etheta,ephi=self.geom.TangentBundle(self.rvalhalf[a,:,:])
		if polar:
		    ax.plot(angles[:,1],angles[:,0], color=defectsmap(a),linestyle='-',marker='')
		else:
		    ax1.plot(self.track4,angles[:,0], color=defectsmap(a),linestyle='-',marker='')
		    ax2.plot(self.track4,angles[:,1], color=defectsmap(a),linestyle='-',marker='')
	if self.isBand:
	    angles=np.zeros((len(self.track2),2))
	    for a in range(2):
		angles[:,0],angles[:,1],etheta,ephi=self.geom.TangentBundle(self.rvalone[a,:,:])
		if polar:
		    ax.plot(angles[:,1],angles[:,0], color=defectsmap(a+2),linestyle='-',marker='')
		else:
		    ax1.plot(self.track2,angles[:,0], color=defectsmap(a+2),linestyle='-',marker='')
		    ax2.plot(self.track2,angles[:,1], color=defectsmap(a+2),linestyle='-',marker='')
	if polar:
	    ax.set_rmax(np.pi)
	    ax.grid(True)
	else:
	    ax1.set_xlabel('time')
	    ax2.set_xlabel('time')
	    ax1.set_ylabel('theta')
	    ax2.set_ylabel('phi')
	plt.title(title)
	return onefig
	  
      
    def trackFourDefects(self,rvalhalf0):
	# Prepare a 4xsnapx3 tracks matrix
	self.rvalhalf=np.zeros((4,len(rvalhalf0),3))
	# Stick in the first snapshot
	self.rvalhalf[:,0,:]=rvalhalf0[0]
	# Track the remaining trajectory (ignoring jumps, ignoring splits, just picking the least distance in each case)
	for n in range(1,len(rvalhalf0)):
	    # Calculate all permutations of distances to previous one and simply pick total minimum
	    idxkeep=[]
	    for b in range(4):
		dist=np.empty((4,))
		for a in range(4):
		    dist[a]=np.sqrt(np.sum((rvalhalf0[n][a,:]-self.rvalhalf[b,n-1,:])**2))
		idxsort=np.argsort(dist)
		if b>0:
		    v=0
		    hmm=np.nonzero(np.array(idxkeep)==idxsort[v])[0]
		    while len(hmm)>0:
			v+=1
			hmm=np.nonzero(np.array(idxkeep)==idxsort[v])[0]
		    idx=idxsort[v]
		else:
		    idx=idxsort[0]
		idxkeep.append(idx)
		# Stick the one we found in with label b
		self.rvalhalf[b,n,:]=rvalhalf0[n][idx,:]
	

    def trackTwoDefects(self,rvalone0):
	  # Prepare a 2xsnap#3 tracks matrix
	  self.rvalone=np.zeros((2,len(rvalone0),3))
	  # Stick in the first snapshot
	  self.rvalone[:,0,:]=rvalone0[0]
	  # Track the remaining trajectory (ignoring jumps, ignoring splits, just picking the least distance in each case)
	  for n in range(1,len(rvalone0)):
	      # Calculate all permutations of distances to previous one and simply pick total minimum
	      for b in range(2):
		  dist=np.empty((2,))
		  for a in range(2):
		      dist[a]=np.sqrt(np.sum((rvalone0[n][a,:]-self.rvalone[b,n-1,:])**2))
		  idx=np.argmin(dist)
		  #print dist
		  # Stick the one we found in with label b
		  self.rvalone[b,n,:]=rvalone0[n][idx,:]
	  

    def findroot(self,pointer,i):
	if pointer[i]<0:
	  ptri=i
	else:
	  ptri=self.findroot(pointer,pointer[i])
	return ptri
  
    # Merging very close defects
    def mergeDefects(self,rmerge,verbose=False):
	self.ndefects1=[]
	self.defects1=[]
	for u in range(len(self.ndefects)):
	    if verbose:
		print "Working on snapshot " + str(u)
	    # Dig out defects, associate position, charge and angles
	    ndef=len(self.defects[u])
	    charge0=[]
	    rval0=np.zeros((ndef,3))
	    charge0=np.zeros((ndef,))
	    for n in range(ndef):
		charge0[n]=self.defects[u][n][0]
		rval0[n,:]=self.defects[u][n][1:4]
	    # Compute pair distances and merge
	    # The old fashioned language here is because I don't want to encounter pairs twice
	    # that just makes trouble ...
	    if verbose:
		print charge0
	    close=[]
	    for n in range(ndef):
		for m in range(n+1,ndef):
		  dist=np.sqrt(np.sum((rval0[n,:]-rval0[m,:])**2))
		  if dist<rmerge:
		    close.append([n,m])
	    if close==[]:
		rval=list(rval0)
		charge=list(charge0)
	    else:
		if verbose:
		    print close
		# need to do a cluster merge on close ones (mostly 3s)
		# Newman-Ziff like ... kind of
		# Make every defect the root of a cluster of size 1
		size=[1 for k in range(ndef)]
		pointer=[-1 for k in range(ndef)]
		# Stick in the known bonds
		for pair in close:
		    if pointer[pair[1]]==-1:
			pointer[pair[1]]=pair[0]
		    else:
			pointer[pair[0]]=pointer[pair[1]]
		if verbose:
		    print pointer
		# Now recursively relabel pointers and sizes. 
		# Loop over all defects
		for k1 in range(ndef):
		    #print "Starting " + str(k1)
		    # Find the root of my defect. That's either myself, or the root of a cluster
		    k2=self.findroot(pointer,k1)
		    # If we are part of a recognizable tree, and not at its root, i.e if it's not myself
		    if k2!=k1:
		      size[k2]+=size[k1]
		      size[k1]=0
		      pointer[k1]=k2     
		if verbose:
		    print pointer
		    print size
		# reconstructing my list grouped by clusters now
		rval=[]
		charge=[]
		for k in range(ndef):
		    thiscluster = [index for index,value in enumerate(pointer) if value==k]
		    if size[k]>0:
			thiscluster.append(k)
		    if thiscluster!=[]:
			newcharge=np.sum(charge0[thiscluster])
			if newcharge!=0.0:
			    newr0=np.mean(rval0[thiscluster,:],axis=0)
			    # normalize again onto the sphere
			    newr0=self.geom.R*newr0/np.sqrt(np.sum(newr0**2))
			    rval.append(newr0)
			    charge.append(newcharge)
		if verbose:
		    print charge
	    self.ndefects1.append(len(charge))
	    self.defects1.append([])
	    for n in range(self.ndefects1[u]):
		tmp=[charge[n]]
		tmp.extend(rval[n])
		self.defects1[u].append(tmp)
		
   
    def makeHistogram(self,histbin,title,merged=False):

	# Initialize pair permutation histograms
	nbin=len(histbin)
	halfplusplus=np.zeros((nbin-1,))
	halfminmin=np.zeros((nbin-1,))
	oneone=np.zeros((nbin-1,))
	halfplusmin=np.zeros((nbin-1,))
	onehalfplus=np.zeros((nbin-1,))
	onehalfmin=np.zeros((nbin-1,))
	# Temporal loop over defect list
	if merged:
	    nsnaps=len(self.defects1)
	else:
	    nsnaps=len(self.defects)
	for u in range(0,nsnaps):
	    # Dig out defects, associate position, charge and angles
	    if merged:
		ndef=len(self.defects1[u])
	    else:
		ndef=len(self.defects[u])
	    charge=[]
	    rval=np.zeros((ndef,3))
	    if merged:
		for n in range(ndef):
		    charge.append(self.defects1[u][n][0])
		    rval[n,:]=self.defects1[u][n][1:4]
	    else:
		for n in range(ndef):
		    charge.append(self.defects[u][n][0])
		    rval[n,:]=self.defects[u][n][1:4]
	    # Sort them into defects of charge +1/2, -1/2 and 1.0
	    ishalf=[index for index,value in enumerate(charge) if value==0.5]
	    ismhalf=[index for index,value in enumerate(charge) if value==-0.5]
	    isone=[index for index,value in enumerate(charge) if value==1.0]
	    
	    # Add a histogram of pair angles (everything with everything in principle, for the three common types)
	    # Just simply invert the dot product ...)
	    if len(ishalf)>0:
		ang_halfplusplus=180.0/np.pi*np.arccos(np.einsum('ij,kj->ik',rval[ishalf,:],rval[ishalf,:])/self.geom.R**2)
		nanlist=list(np.where(np.isnan(ang_halfplusplus))[0])
		ang_halfplusplus[nanlist]=0.0
		hist, bin_edges = np.histogram(ang_halfplusplus,bins=histbin, density=True)
		if len(ang_halfplusplus)>0:
		    halfplusplus+=hist*len(ang_halfplusplus)**2 
	    if len(ismhalf)>0:
		ang_halfminmin=180.0/np.pi*np.arccos(np.einsum('ij,kj->ik',rval[ismhalf,:],rval[ismhalf,:])/self.geom.R**2)
		nanlist=list(np.where(np.isnan(ang_halfminmin))[0])
		ang_halfminmin[nanlist]=0.0
		hist, bin_edges = np.histogram(ang_halfminmin,bins=histbin, density=True)
		if len(ang_halfminmin)>0:
		    halfminmin+=hist*len(ang_halfminmin)**2 
	    if len(isone)>0:
		ang_oneone=180.0/np.pi*np.arccos(np.einsum('ij,kj->ik',rval[isone,:],rval[isone,:])/self.geom.R**2)
		nanlist=list(np.where(np.isnan(ang_oneone))[0])
		ang_oneone[nanlist]=0.0
		hist, bin_edges = np.histogram(ang_oneone,bins=histbin, density=True)
		if len(ang_oneone)>0:
		    oneone+=hist*len(ang_oneone)**2 
	    if (len(ishalf)>0) and (len(ismhalf)>0):
		ang_halfplusmin=180.0/np.pi*np.arccos(np.einsum('ij,kj->ik',rval[ishalf,:],rval[ismhalf,:])/self.geom.R**2)
		nanlist=list(np.where(np.isnan(ang_halfplusmin))[0])
		ang_halfplusmin[nanlist]=0.0
		hist, bin_edges = np.histogram(ang_halfplusmin,bins=histbin, density=True)
		if len(ang_halfplusmin)>0:
		    halfplusmin+=hist*len(ang_halfplusmin)**2 
	    if (len(isone)>0) and (len(ishalf)>0):
		ang_onehalfplus=180.0/np.pi*np.arccos(np.einsum('ij,kj->ik',rval[isone,:],rval[ishalf,:])/self.geom.R**2)
		nanlist=list(np.where(np.isnan(ang_onehalfplus))[0])
		ang_onehalfplus[nanlist]=0.0
		hist, bin_edges = np.histogram(ang_onehalfplus,bins=histbin, density=True)
		if len(ang_onehalfplus)>0:
		    onehalfplus+=hist*len(ang_onehalfplus)**2 
	    if (len(isone)>0) and (len(ismhalf)>0):
		ang_onehalfmin=180.0/np.pi*np.arccos(np.einsum('ij,kj->ik',rval[isone,:],rval[ismhalf,:])/self.geom.R**2)
		nanlist=list(np.where(np.isnan(ang_onehalfmin))[0])
		ang_onehalfmin[nanlist]=0.0
		hist, bin_edges = np.histogram(ang_onehalfmin,bins=histbin, density=True)
		if len(ang_onehalfmin)>0:
		    onehalfmin+=hist*len(ang_onehalfmin)**2 
				    
	# histogram figure
	histfig=plt.figure( figsize=(6, 4))
	# Exclude the self-bits of the histogram
	halfplusplus[0]=0
	halfminmin[0]=0
	oneone[0]=0
	halfplusmin[0]=0
	onehalfplus[0]=0
	onehalfmin[0]=0
	isdata=[index for index,value in enumerate(halfplusplus) if value>0]
	#plt.plot(histbin[isdata],halfplusplus[isdata]/nsnaps,'.-k',label='+1/2 +1/2')
	if len(isdata)>0:
	    plt.semilogy(histbin[:(nbin-1)],halfplusplus/nsnaps,'.-k',label='+1/2 +1/2')
	isdata=[index for index,value in enumerate(halfminmin) if value>0]
	#plt.plot(histbin[isdata],halfminmin[isdata]/nsnaps,'.-r',label='-1/2 -1/2')
	if len(isdata)>0:
	    plt.semilogy(histbin[:(nbin-1)],halfminmin/nsnaps,'.-',color=[1.0, 0,0],label='-1/2 -1/2')
	isdata=[index for index,value in enumerate(oneone) if value>0]
	#plt.plot(histbin[isdata],oneone[isdata],'.-g',label='+1 +1')
	if len(isdata)>0:
	    plt.semilogy(histbin[:(nbin-1)],oneone,'.-',color=[0, 1.0,0.0],label='+1 +1')
	isdata=[index for index,value in enumerate(halfplusmin) if value>0]
	#plt.plot(histbin[isdata],halfplusmin[isdata]/nsnaps,'.-',color=[0.5, 0, 0],label='+1/2 -1/2')
	if len(isdata)>0:
	    plt.semilogy(histbin[:(nbin-1)],halfplusmin/nsnaps,'.-',color=[0.5, 0, 0],label='+1/2 -1/2')
	isdata=[index for index,value in enumerate(onehalfplus) if value>0]
	#plt.plot(histbin[isdata],onehalfplus[isdata]/nsnaps,'.-',color=[0, 0.5, 0],label='+1 +1/2')
	if len(isdata)>0:
	    plt.semilogy(histbin[:(nbin-1)],onehalfplus/nsnaps,'.-',color=[0, 0.5, 0],label='+1 +1/2')
	isdata=[index for index,value in enumerate(onehalfmin) if value>0]
	#plt.plot(histbin[isdata],onehalfmin[isdata]/nsnaps,'.-',color=[0.5, 0.5, 0],label='+1 -1/2')
	if len(isdata)>0:
	    plt.plot(histbin[:(nbin-1)],onehalfmin/nsnaps,'.-',color=[0.5, 0.5, 0],label='+1 -1/2')
	plt.xlabel('angle')
	plt.ylabel('pair frequency')
	plt.legend(loc=2)
	plt.title(title)
	return histfig,halfplusplus,halfminmin,oneone,halfplusmin,onehalfplus,onehalfmin

#===============================================================================

    
    



