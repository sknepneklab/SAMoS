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
import numpy as np
import matplotlib.pyplot as plt
# Samos analysis libraries
try:
	from Geometry import *
	from Configuration import *
	from Hessian2d import *
	from Glassy import *
except:
	print "Could not load SAMoS analysis modules."
	print "Please make sure that the shell variable PYTHONPATH"
	print "includes the location of the SAMoS analysis modules."
	print "These modules are typically located in the analysis"
	print "directory inside you SAMoS installation directory."
	sys.exit(1)

if 'BASEDIR' in os.environ:
  filebase = os.environ["BASEDIR"]+'/'
else:
  print "Warning!"
  print "Base directory has not been set. Using current directory."
  print "You can set base directory by setting the shell environment variable BASEDIR."
  print "In bash, e.g., you can type export BASEDIR=$HOME/path/to/base"
  print
  filebase = './'
if 'OUTDIR' in os.environ:
  fileout = os.environ['OUTDIR']+'/'
else:
  print "Warning!"
  print "Output directory has not been set. Using current directory."
  print "You can set output directory by setting the shell environment variable OUTDIR."
  print "In bash, e.g., you can type export OUTDIR=$HOME/path/to/output"
  print
  fileout = './'

# form of the analytical correlation function
def Corrfun(q,v0,mu,zeta,tau):
	return 0.5*v0**2/(1+mu*tau/zeta*q**2)



tempval  = ['0.01']
nuval    = ['0.01']
taus     = np.array([200])
zeta     = 1.0

ntemp = len(tempval)
nnu   = len(nuval)

# careful counting of the number of q values that fit into this box
# from output of qrad (either of velocity Fourier transforms)
# need to know this now to save output as numpy array
nqout = 224

for temp in tempval:
	ncount=0
	for nu in nuval:
		# get the directory where the data is located
		confdir = filebase +  "vertex/data_temp_" + temp +"/data_nu_" + nu + "_plane/"
		# get the intitial configuration
		conffile="vertex_temp_" + temp + "_nu_" + nu + ".conf"
		# prefix of the file names, and location of the particle radii (required for internal compatibility, not used for vertex model)
		prefix = 'cell_0'
		radiusfile = 'vertex_input.dat'
		# Read in all of the simulated data. Skip is normally 250, but second half of simulation here
		# Usetype = 1 - use only the interior of the system, exclude the boundary particles from the analysis
		sim = SimRun(confdir,conffile,prefix,radiusfile,100,False,True,True,usetype=1)
		
		# Glassy dynamics quantities
		# compute the mean square displacement
		tplot,msd = sim.getMSD(False)
		# compute the self-intermediate scattering function (at a number of q; we use the last point)
		qval = np.linspace(0,2*np.pi,40)
		SelfInt = np.zeros((len(qval),sim.Nsnap))
		for q in range(len(qval)):
			qvalintermediate = qval[q]*np.array([1,1,0])
			tval,SelfInt[q,:] = sim.SelfIntermediate(qvalintermediate,False)
                
                # Do the Fourier transform of the simulated velocities
		print "Starting comparison to actual simulation data velocity correlations"
		Fouriervel = np.zeros((nqout,))
		u = 0
		for k in range(0,sim.Nsnap,1):
      if k % 10 == 0:
				print k
			qrad,valrad,Sqrad = sim.FourierTransVel(k,5.0,False,True,100)
			Fouriervel += Sqrad
			u += 1
		Fouriervel /= u
		
		# use the values of the bulk and shear modulus from Sussmann et al
		B  = 7.0
		mu = 0.5
		print 'mu = ', mu
		print 'B = ', B
		tau = 2.0/float(nu)
		v = np.sqrt(float(nu)*float(temp))
		# Analytic prediction (equation 10 of the paper)
		Gofq_analytics = Corrfun(qrad,v,B+mu,zeta,taus[ncount])+Corrfun(qrad,v,mu,zeta,taus[ncount])
		
		# plotting both correlation functions (i.e. one of the temperature points of Figure 3c)
		plt.figure()
		plt.loglog(qrad,Fouriervel/v**2,'o',color='k',label='simulation')
		plt.loglog(qrad,Gofq_analytics/v**2,'.--',color='b',label='analytics')
		plt.xlabel('q')
		plt.ylabel('Fourier/v2')
		plt.title('Fourier correlations, tau = ' + str(tau))
		
		data = {'tempval': tempval, 
		        'nuval': nuval,
						'Nsnap': sim.Nsnap, 
						'tplot': tplot, 
						'msd': msd, 
						'qval': qval, 
						'SelfInt': SelfInt, 
						'Gofq_analytics': Gofq_analytics, 
						'qrad': qrad,
						'Fouriervel': Fouriervel}

		outpickle=fileout +"Vertex_temp_" + temp + "_nu_" + nu + ".p"
		pickle.dump(data,open(outpickle,'wb'))
	ncount+=1


plt.show()
		







