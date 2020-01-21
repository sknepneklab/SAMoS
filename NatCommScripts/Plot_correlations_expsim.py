# *****************************************************************************
# *
# *  This Python script is a part of tha analysis of the data published in 
# *  the paper: "Universal motion patterns in confluent cell monolayers"
# *  by Silke Henkes, Kaja Kostanjevec, J. Martin Collinson, Rastko Sknepnek, 
# *  and Eric Bertin, Jounral name, vol, page (2019).
# *
# *  Please refer to the document Computational_summary.pdf for a detailed
# *  description of the tasks performed by this script.
# * 
# *****************************************************************************

import random
import sys, os, glob
import pickle as pickle
import copy as cp
import numpy as np
import scipy as sp
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as ptch
import matplotlib.lines as lne
from matplotlib.colors import LinearSegmentedColormap
import argparse
from matplotlib import rc

rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)

matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['axes.linewidth'] = 2
matplotlib.rcParams['xtick.major.size'] = 16
matplotlib.rcParams['xtick.minor.size'] = 0
matplotlib.rcParams['ytick.major.size'] = 16
matplotlib.rcParams['ytick.minor.size'] = 0
matplotlib.rcParams['font.size']=24.0
matplotlib.rcParams['legend.fontsize']=18.0

cdict = {'red':   [(0.0,  0.0, 0.5),
				           (0.35,  1.0, 0.75),
                   (0.45,  0.75, 0.0),
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
                   



# analytical correlation function with a single fit parameter mutauzeta (for experiments)
def Corrfun(q,v0,mutauzeta):
    return 0.5*v0**2/(1+mutauzeta*q**2)

# analytical velocity scaling with fit parameters mutauzeta and a2, the area cutoff (for experiments)
def Vscaling(qmax,qmin,mutauzeta,a2):
    return np.sqrt(a2/(8*np.pi)*1.0/mutauzeta*np.log((1.0+mutauzeta*qmax**2)/(1.0+mutauzeta*qmin**2)))

# analytical velocity correlation function with all its parameters (for simulations)
def Corrfun2(q,v0,mu,zeta,tau):
    return 0.5*v0**2/(1+mu*tau/zeta*q**2)

# analytical velocity scaling with all its parameters and a2, the area cutoff (for simulations)
def Vscaling2(qmax,qmin,mu,zeta,tau,a2):
    return np.sqrt(a2/(8*np.pi)*1.0/(mu*tau/zeta)*np.log((1.0+mu*tau/zeta*qmax**2)/(1.0+mu*tau/zeta*qmin**2)))

# velocity autocorrelation analytical predictions with all the parameters
def Vauto_integrand(t,q,mu,zeta,tau):
    return q/(4*np.pi*tau)*(mu*q**2*np.exp(-mu/zeta*q**2*t)-zeta/tau*np.exp(-t/tau))/(mu**2*q**4-(zeta/tau)**2)


if 'OUTDIR' in os.environ:
  fileout = os.environ['OUTDIR']+'/'
else:
  print "Warning!"
  print "Output directory has not been set. Using current directory."
  print "You can set output directory by setting the shell environment variable OUTDIR."
  print "In bash, e.g., you can type export OUTDIR=$HOME/path/to/output"
  print
  fileout = './'

# final parameters of the matching simulations
nu = '0.8'
k = '55'
v0 = '90'
phi = '0.95'
rdeath = '0.01'
mult = '0.05'
# area per particle
a2 = 380
# we decided to divide it by pi, taking particle radius as the fundamental length scale
a2 = a2/np.pi
N = 1434
        
m = 0

plt.figure()

### Soft particles without division
outpickle = fileout + "Velcorr_" + "phi" + phi + "_nu" + nu +"_v0" + v0 + "_k" + k + "_expsim.p"
data = pickle.load(open(outpickle, "rb"))
locals().update(data)
plt.loglog(qrad,Fouriervel/np.mean(vav)**2,'.-',color='r',label='soft')

### Soft particles with division
outpickle = fileout + "Velcorr_" + "rdeath" + rdeath + "_nu" + nu +"_v0" + v0 + "_ratio" + mult + "_k" + k + "_expsim.p"
data = pickle.load(open(outpickle, "rb"))
locals().update(data)
plt.loglog(qrad,Fouriervel/np.mean(vav)**2,'.-',color='g',label='dividing')

# and the expected analytical line for soft particles
tau = 2.0/float(nu)
zeta = 1.0
mu = 0.51*11**2*float(k)
B = 1.684*11**2*float(k)
vscal = Vscaling2(0.25,0,B+mu,zeta,tau,a2)+Vscaling2(0.25,0,mu,zeta,tau,a2)
v02 = vscal**(-1)
qrad2 = np.logspace(-2.5,-0.5,100)
Gofq_analytics_real = Corrfun2(qrad2,v02,B+mu,zeta,tau)+Corrfun2(qrad2,v02,mu,zeta,tau)
plt.plot(qrad2,Gofq_analytics_real,'--',color='r')

## Vertex model
phi = '3.6'
outpickle = fileout + "Velcorr_" + "phi" + phi + "_nu" + nu +"_v0" + v0 + "_k" + k + "vertex_expsim.p"
data = pickle.load(open(outpickle, "rb"))
locals().update(data)
plt.loglog(qrad,Fouriervel/np.mean(vav)**2,'.-',color='b',label='vertex')
# and the expected analytical line ...
tau = 2.0/float(nu)
zeta = 1.0
mu = 0.5*11**2*float(k)
B = 7.0*11**2*float(k)
vscal = Vscaling2(0.25,0,B+mu,zeta,tau,a2)+Vscaling2(0.25,0,mu,zeta,tau,a2)
print vscal
v02 = vscal**(-1)
print v02
qrad2 = np.logspace(-2.5,-0.5,100)
Gofq_analytics_real = Corrfun2(qrad2,v02,B+mu,zeta,tau)+Corrfun2(qrad2,v02,mu,zeta,tau)
plt.plot(qrad2,Gofq_analytics_real,'--',color='b')

# expected ratio of bulk and shear modulus 
rat = 4.3
muBtauzeta = 10000
# This is in fact the ratio v0/vav
# Should come out of the velocity scaling
vscal = Vscaling(0.25,0,muBtauzeta,a2)+Vscaling(0.25,0,rat*muBtauzeta,a2)
#print vscal
# but doesn't after putting the a^2 in there
v0 = 1./vscal
print 'v0 = ', v0
#Gofq_analytics=Corrfun(qrad,v0,muBtauzeta)+Corrfun(qrad,v0,mutauzeta)
qrad2 = np.logspace(-2.5,-0.5,100)
Gofq_analytics = Corrfun(qrad2,v0,muBtauzeta) + Corrfun(qrad2,v0,rat*muBtauzeta)
plt.plot(qrad2,Gofq_analytics,'k-',label='exp. fit')


plt.xlabel('q')
plt.ylabel('Fourier')
plt.title('Simulations')
plt.legend()
plt.xlim(0.003,0.4)
plt.ylim(0.02,300)

######=================================================
# Correlation functions from experiment,
nx = 54
ny = 40
N = nx*ny

topdir = "./"
expnamelist = ['7']
testmap = LinearSegmentedColormap('test',cdict,N=len(expnamelist)+1) 

en = 0
plt.figure()
v0estimate = np.zeros((len(expnamelist),))
for expname in expnamelist:
    filename = topdir + "PIVresults_exp7.p"
    data = pickle.load(open(filename))
    qrad = data["qrad"]
    Sqrad = data["Sqrad"]
    velav = data["velav"]
    vav = np.mean(velav)
    #print vav
    v0estimate[en] = np.sqrt(Sqrad[0]*N)
    plt.loglog(qrad[:len(Sqrad)],N*Sqrad,'o',color = testmap(en),label = expname)
    en += 1
#0.5*v0**2/(1+mu*tau/zeta*q**2)
# assume the bulk modulus and the shear moduli scale as in our simulations, i.e. B = 7 and mu = 0.4 in reduced units
rat = 4.3
muBtauzeta = 10000
# This is in fact the ratio v0/vav
# Should come out of the velocity scaling
vscal = Vscaling(0.25,0,muBtauzeta,a2)+Vscaling(0.25,0,rat*muBtauzeta,a2)

v0 = 1.0/vscal
print 'v0 = ', v0
qrad2 = np.logspace(-2.5,-0.5,100)
Gofq_analytics = Corrfun(qrad2,v0,muBtauzeta) + Corrfun(qrad2,v0,rat*muBtauzeta)
plt.plot(qrad2, Gofq_analytics, 'k-', label='10000')
# comparison lines
muBtauzeta = 5000
vscal = Vscaling(0.25,0,muBtauzeta,a2)+Vscaling(0.25,0,rat*muBtauzeta,a2)
v0 = vscal**(-1)
print v0
Gofq_analytics = Corrfun(qrad2,v0,muBtauzeta)+Corrfun(qrad2,v0,rat*muBtauzeta)
plt.plot(qrad2,Gofq_analytics,'r--',label='5000')
muBtauzeta = 20000
vscal = Vscaling(0.25,0,muBtauzeta,a2)+Vscaling(0.25,0,rat*muBtauzeta,a2)
v0 = 1.0/vscal
print 'v0 = ', v0
Gofq_analytics = Corrfun(qrad2,v0,muBtauzeta)+Corrfun(qrad2,v0,rat*muBtauzeta)
plt.plot(qrad2, Gofq_analytics, 'g--',label='20000')



plt.xlabel('q')
plt.ylabel('Fourier')
plt.title('Experiment')
plt.legend(loc = 3)
plt.xlim(0.003,0.4)
plt.ylim(0.02,300)


plt.show()
		







