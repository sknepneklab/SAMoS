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
sys.path.insert(0, '/home/silke/Documents/SAMoS/source_git/apcs/utils')

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

cdict = {'red':   [(0.0,  0.05, 0.5),
				   (0.3,  1.0, 0.75),
                   (0.5,  0.75, 0.0),
                   (1.0,  0.0, 0.0)],

         'green': [(0.0,  0.0, 0.0),
				   (0.35,  0.0, 0.5),
                   (0.5, 1.0, 1.0),
                   (0.8,  0.5, 0.0),
                   (1.0,  0.0, 0.0)],

         'blue':  [(0.0,  0.0, 0.0),
                   (0.5,  0.0, 0.0),
                   (0.7, 0.5, 1.0),
                   (1.0,  0.25, 0.0)]}
         
picklefolder_silke="/home/silke/Documents/SAMoS/Nematic/nematic/twoAligning/picklefiles" 

# Hack to avoid the whole parameter read-in saga
# and use the geometry class anyways
class param:
    def __init__(self,R):
        self.r=R

testmap=LinearSegmentedColormap('test',cdict,N=4) 

# k=1 data
# Rval="5 10 15 20 25 30 40 50"
# k="1"
# Jval="0.01 0.1 1.0 10.0"
# Tau="100"
# vval="0.025 0.75 1.25"
# vval="0.01 0.05 0.1 0.25"
# step="50"
# nsteps="200"

##vval="0.1 0.5 1.5 2.0 2.5 3.0 5.0"
##step="10"
##nsteps="200"
#Rval=[30]
#vval=['0.1', '0.5', '1.5', '2.0', '2.5', '3.0', '5.0']
##vval=['0.1']
#Jval=['10.0']
#k='1'
#Tau='100'
#step='10'
#nsteps='200'


Rval=[30]
vval=['0.025', '0.05', '0.1', '0.25', '0.75', '1.25']
#vval=['0.1']
Jval=['10.0']
k='1'
Tau='100'
step='50'
nsteps='200'

#Rval=[30]
##vval=['0.075','0.1', '0.125','0.15','0.175','0.2','0.225','0.25','0.275','0.3','0.325','0.35','0.375','0.4','0.425','0.45','0.475','0.5']
#vval=['0.475','0.5']
#Jval=['1.0']
#k='1'
#Tau='100.0'
#step='1'
#nsteps='1000'


tracking=False

nj=0
for J in Jval:
    nr=0
    for R in Rval:
        nv=0
        par=param(R)
        geom=GeometrySphere(par)
        for v in vval:
            plt.figure()
            ax = plt.subplot(111, polar=True)


            filename=picklefolder_silke+"/k_" + k + "/defect_data_R" + str(R) + "_J" + J + "_v" + v+ "_Tau" + Tau + "_k" +k +"_step" +step+"_nsteps" + nsteps +".p"
            print filename
            data=pickle.load(open(filename))
            ndefects=data["numdefects_n"]
            defects=data["defects_n"]

            # Get the Euler angles for all the defects, and then plot, according to charge
            rvalhalf0=[]
            rvalone0=[]
            for u in range(100,len(defects)):
                ndef=len(defects[u])
                charge=[]
                angles=np.zeros((ndef,2))
                rval=np.zeros((ndef,3))
                for n in range(ndef):
                    charge.append(defects[u][n][0])
                    rval[n,:]=defects[u][n][1:4]
                ishalf=[index for index,value in enumerate(charge) if value==0.5]
                ismhalf=[index for index,value in enumerate(charge) if value==-0.5]
                isone=[index for index,value in enumerate(charge) if value==1.0]
                if len(ishalf)==4:
                    rvalhalf0.append(rval[ishalf,:])
                if len(isone)==2:
                    rvalone0.append(rval[isone,:])
                angles[:,0],angles[:,1],etheta,ephi=geom.TangentBundle(rval)
                ax.plot(angles[ishalf,1],angles[ishalf,0], color='k',linestyle='none',marker='.')
                ax.plot(angles[ismhalf,1],angles[ismhalf,0], color='r',linestyle='none',marker='.')
                ax.plot(angles[isone,1],angles[isone,0], color='g',linestyle='none',marker='.')
             ##data={'J':params.J,'v':params.v0,'k':params.pot_params['k'],'defects_n':defects_n_out,'defects_v':defects_v_out,'numdefects_n':numdefects_n_out,'numdefects_v':numdefects_v_out}
            if tracking:
                # Sort out the tracks if possible
                if len(rvalhalf0)>0:
                    rvalhalf=np.zeros((4,len(rvalhalf0),3))
                    rvalhalf[:,0,:]=rvalhalf0[0]
                    for n in range(1,len(rvalhalf0)):
                        # Calculate all permutations of distances to previous one and pick total minimum?
                        for b in range(4):
                            dist=np.empty((4,))
                            for a in range(4):
                                dist[a]=np.sqrt(np.sum((rvalhalf0[n][a,:]-rvalhalf[b,n-1,:])**2))
                            idx=np.argmin(dist)
                            #print dist
                            rvalhalf[b,n,:]=rvalhalf0[n][idx,:]
                    angles=np.zeros((len(rvalhalf0),2))
                    for a in range(4):
                        angles[:,0],angles[:,1],etheta,ephi=geom.TangentBundle(rvalhalf[a,:,:])
                        ax.plot(angles[:,1],angles[:,0], color=testmap(a),linestyle='-',marker='')
                        #ax.plot(angles[:,1],angles[:,0], color='k',linestyle='-',marker='')
                if len(rvalone0)>0:
                    rvalone=np.zeros((2,len(rvalone0),3))
                    rvalone[:,0,:]=rvalone0[0]
                    for n in range(1,len(rvalone0)):
                        # Calculate all permutations of distances to previous one and pick total minimum?
                        for b in range(2):
                            dist=np.empty((2,))
                            for a in range(2):
                                dist[a]=np.sqrt(np.sum((rvalone0[n][a,:]-rvalone[b,n-1,:])**2))
                            idx=np.argmin(dist)
                            #print dist
                            rvalone[b,n,:]=rvalone0[n][idx,:]
                    angles=np.zeros((len(rvalone0),2))
                    for a in range(2):
                        angles[:,0],angles[:,1],etheta,ephi=geom.TangentBundle(rvalone[a,:,:])
                        ax.plot(angles[:,1],angles[:,0], color='g',linestyle='-',marker='')
            ax.set_rmax(np.pi)
            ax.grid(True)
            #plt.xlabel('theta')
            #plt.ylabel('phi')
            plt.title('v0=' + v+ ', J=' + J + ', R=' + str(R) + ', Tau=' + Tau + 'k=' + k)
            nv+=1
        nr+=1
    nj+=1


plt.show()
    
    



