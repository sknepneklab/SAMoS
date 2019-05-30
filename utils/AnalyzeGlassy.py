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
import argparse
import pickle

from Geometry import *
from Glassy import *
from Writer import *

# WTF is wrong with my plotting??
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="input file (base name) ")
parser.add_argument("-r", "--radii", type=str, help="radii file (initial configuration) ")
parser.add_argument("-c", "--conffile", type=str, help="configuration file")
parser.add_argument("-d", "--directory", type=str, help="input directory")
parser.add_argument("-o", "--output", type=str, help="output directory")
parser.add_argument("-p", "--prefix", type=str, default="glassy",help="prefix for output file")
parser.add_argument("-s", "--skip", type=int, default=0, help="skip this many samples")
parser.add_argument("-t", "--step", type=int, default=1, help="step snapshots with this spacing in correlation function")
parser.add_argument("--tracer", action='store_true', default=False, help="Is this a configuration with tracer particles")
parser.add_argument("--drift", action='store_true', default=False, help="Should I take away the drift?")
parser.add_argument("-u", "--usetype", type=str, default="all", help="Which types of particles should I use for MSD and SelfInt?")
parser.add_argument("--getMSD", action='store_true', default=False, help="Compute mean squared displacement?")
parser.add_argument("--getSelfInt", action='store_true', default=False, help="Compute Self-Intermediate scattering function?")
parser.add_argument("--getFourPoint", action='store_true', default=False, help="Compute four point function?")
parser.add_argument("--getDynStruct", action='store_true', default=False, help="Compute the dynamic structure factor?")
parser.add_argument("--getFourier", action='store_true', default=False, help="Compute the Fourier transformed positions and velocities?")
parser.add_argument("--getVelcorr", action='store_true', default=False, help="Compute velocity correlation function")
parser.add_argument("--getNonGaussian", action='store_true', default=False, help="Compute velocity correlation function")
parser.add_argument("--plot", action='store_true', default=False, help="Plot MSD and correlations")
parser.add_argument("--ignore", action='store_true', default=False, help="Ignore complications like missing potentials for quick result (warning!)")

args = parser.parse_args()
#def __init__(self,directory,conffile,inputfile,radiusfile,skip,tracer=False,ignore=False,takeDrift=False,usetype='all'):
#simL = SimRun(confdir,conffile,prefix,radiusfile,skip,True,False,True)
sim = SimRun(args.directory,args.conffile,args.input,args.radii,args.skip,args.tracer,args.ignore,args.drift,args.usetype)
data={'input':args.input,'configuration':args.conffile}
if args.getMSD:
  #getMSD(self,verbose=True):
	tplot,msd = sim.getMSD(args.plot)
	dataMSD={'tplot':tplot,'msd':msd,}
	data.update(dataMSD)
if args.getSelfInt:
  ##This tends to be too slow except for when tracer particles are used
  #qval=np.linspace(0,np.pi,20)
  #SelfInt=np.zeros((len(qval),sim.Nsnap))
  #for q in range(len(qval)):
          #qvalintermediate=qval[q]*np.array([1,1,0])
          #print qval[q]
          #tval,SelfInt[q,:] = sim.SelfIntermediate(qvalintermediate,args.plot)
  #dataSelfInt={'qval':qval,'tval':tval,'SelfInt':SelfInt}
  #data.update(dataSelfInt)
  # Take a single slice at pi (magnitude)
  qval=np.pi
  SelfInt=np.zeros((sim.Nsnap,))
  qvalintermediate=qval*np.array([1,1,0])/np.sqrt(2)
  tval,SelfInt = sim.SelfIntermediate(qvalintermediate,args.plot)
  dataSelfInt={'qval':qval,'tval':tval,'SelfInt':SelfInt}
  data.update(dataSelfInt)
if args.getFourPoint:
  afourpoint=0.5
  qmax=np.pi
  nmax=50
  #def FourPoint(self,a,qmax=3.14,verbose=True,nmax=20):
  tvalFourPoint, FourPoint=sim.FourPoint(afourpoint,qmax,args.plot,nmax)
  dataFourPoint={'afourpoint':afourpoint,'qmax':qmax,'nmax':nmax,'tvalFourPoint':tvalFourPoint,'FourPoint':FourPoint}
  data.update(dataFourPoint)
if args.getDynStruct:  
  nmax=50
  qmax=np.pi
  omegamax=0.1
  # def getDynStruct(self,qmax,omegamax,verbose=True,nmax=50):
  omega,qrad,DynStruct=sim.getDynStruct(qmax,omegamax,args.plot,nmax)
  dataDynStruct={'omegamax':omegamax,'qrad':qrad,'DynStruct':DynStruct}
  data.update(dataDynStruct)
if args.getFourier:
  Sqvel=np.zeros((107,))
  Sqrad=np.zeros((107,))
  qmaxFourier=4.0
  npts=sim.Nsnap/args.step
  #plt.figure()
  for u in range(0,sim.Nsnap,args.step):
    #qradv,velrad,Sqvel=sim.FourierTransVel(u,qmaxFourier,args.plot)
    #Sqvel+=Sqvel
    qrad2,posrad=sim.FourierTrans(u,qmaxFourier,args.plot)
    Sqrad+=posrad
    #plt.plot(qrad2,posrad)
    #plt.text(qrad2[100],posrad[100],str(u))
  Sqvel/=npts
  Sqrad/=npts
  if args.plot:
    plt.figure()
    plt.plot(qrad2,Sqrad)
    plt.xlabel('q')
    plt.ylabel('S(q)')
    plt.title('Positions - after averaging')
  #dataFourier={'npts':npts,'qmaxFourier':qmaxFourier,'qrad2':qrad2,'qradv':qradv,'Sqrad':Sqrad,'Sqvel':Sqvel}
  dataFourier={'npts':npts,'qmaxFourier':qmaxFourier,'qrad2':qrad2,'Sqrad':Sqrad}
  data.update(dataFourier)
if args.getVelcorr:
  dx=0.1
  xmax=20
  nbins=int(xmax/dx)
  velcorr=np.zeros((nbins,))
  npts=sim.Nsnap/args.step
  for u in range(0,sim.Nsnap,args.step):
    #def getVelcorrSingle(self,whichframe,dx,xmax,verbose=True):
    bins,velcorr0=sim.getVelcorrSingle(u,0.1,20,args.plot)
    velcorr+=velcorr0
  velcorr/=npts
	dataVelcorr={'dx':dx,'xmax':xmax,'nbins':nbins,'bins':bins,'velcorr':velcorr}
	data.update(dataVelcorr)	
if args.getNonGaussian:
    #getMSD(self,verbose=True):
    tplot,msd, kurtosis, nongaussian=sim.getNonGaussian(args.plot)
	dataNonGauss={'tplot':tplot,'msd':msd,'kurtosis':kurtosis,'nongaussian':nongaussian}
	data.update(dataNonGauss)


# Save the output files
outglassy=args.output + args.prefix +'.p'
pickle.dump(data,open(outglassy,'wb'))
plt.show()
if args.plot:
	plt.show()
	

	
	
