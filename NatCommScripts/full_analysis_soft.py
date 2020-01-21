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
	return 0.5*v0**2/(1 + mu*tau/zeta*q**2)


tempval = ['0.01']
nuval   = ['0.01']
taus    = np.array([200])
zeta    = 1.0

ntemp = len(tempval)
nnu   = len(nuval)

# careful counting of the number of q values that fit into this box
# from output of qrad (either of velocity Fourier transforms)
# need to know this now to save output as numpy array
nqout = 224

step = 1
for temp in tempval:
	ncount = 0
	for nu in nuval:
		print "=================================================================================================================="
		print "Starting Work on simulation temp=" + temp + ", noise=" + nu 
		# Get the correct equilibrated positions
		equilpositions = filebase + "soft/data_temp_" + temp +"/data_nu_" + nu + "_plane/equilibration_0000100000.dat"
		# Get the parameters for that one (the initial ones)
		conffile="plane_temp_" + temp + "_nu_" + nu + ".conf"
		param = Param(filebase + "soft/data_temp_" + temp +"/data_nu_" + nu + "_plane/"+conffile)
		# Construct the configuration based on these positions
		folder = filebase +  "soft/data_temp_" + temp +"/data_nu_" + nu + "_plane/"
		prefix="plane_temp_"+temp+"_nu_"+nu+"_00"
		# radii are / were not contained in all of the output text files. Retrieve them from the input configuration
		radiusfile="plane_temp_"+temp+"_nu_"+nu+".txt"
		
		print "Creating Hessian and computing modes"
		# Use equilibrated positions to create configuration
		equilpos = Configuration(param,equilpositions)
		# feed this into the Hessian routine 
		hess=Hessian2d(equilpos)
		# and create the Hessian matrix
		hess.makeMatrix()
		
		# Compute the Fourier space orientationally averaged linear response, longitudinal and transverse
		print "Computing Fourier space response for moduli"
		qrad2, longitudinal2, transverse2=hess.getModuli(1.5, False)
		
		# now diagonalise the matrix to obtain the eigenvalues and eigenvectors (normal modes)
		hess.getModes()
		
		# compute the Fourier transform of the modes up to qmax= 5.0 (about an inverse particle radius)
		print "Computing Fourier mode projections"
		Fourier_qproj = np.zeros((2*hess.N+1, nqout))
		u = 0
		# replace 100 by 1 to properly average over all modes: The full projection will take several hours
		for k in range(0,2*hess.N,100):
			if k % 10 == 0:
				print k
			qrad,Fourier_qproj[u,:] = hess.ModesFourier(k,5.0,False)
			u += 1
		
		# Let's compute the Fourier correlation function G(q) out of this like in our draft
		# Full numerical version
		print "Computing velocity correlations from normal modes"
		tau = 2.0/float(nu)
		v = np.sqrt(float(nu)*float(temp))
		Gofq = np.zeros((len(qrad),))
		n = 0
		for q in qrad:
			# Equation 6 of the paper
			Gofq[n] = np.sum(v**2/(2.0*(1.0+hess.eigval*tau/zeta))*(Fourier_qproj[:len(hess.eigval),n]))
			n += 1
		
		# Combine this with the Fourier transform of the actual data
		print "Starting comparison to actual simulation data velocity correlations"
		# for test case: skip first half, usually this is set to 250
		sim = SimRun(folder,conffile,prefix,radiusfile,100,False,False,True)
		
		# First compute the glassy analysis data
		# Mean square displacement
		tplot,msd = sim.getMSD(False)
		# Self-intermediate scattering function (at a number of q; we use the last point for Figure 2)
		qval = np.linspace(0,2*np.pi,40)
		SelfInt = np.zeros((len(qval),sim.Nsnap))
		for q in range(len(qval)):
			qvalintermediate = qval[q]*np.array([1,1,0])
			print qval[q]
			tval,SelfInt[q,:] = sim.SelfIntermediate(qvalintermediate,False)
			
		# compute the correlation function and average it over snapshots
		Fouriervel = np.zeros((nqout,))
		u = 0
		for k in range(0,sim.Nsnap,1):
			if k % 20 == 0:
				print k
			qrad,valrad,Sqrad = sim.FourierTransVel(k,5.0,False)
			Fouriervel += Sqrad
			u += 1
		Fouriervel /= u
			
		# get bulk and shear moduli from fit of longitudinal and transverse response to q^2
		nout = 15
		lfit = np.polyfit(qrad2[:nout],longitudinal2[:nout],1)
		tfit = np.polyfit(qrad2[:nout],transverse2[:nout],1)
		# mu is prefactor of transverse
		mu = tfit[0]
		# prefactor of longitudinal is B + mu
		B =lfit[0] - mu
		print 'mu = ', mu
		print 'B = ', B
		# Analytic prediction (equation 10 of the paper)
		Gofq_analytics = Corrfun(qrad,v,B+mu,zeta,taus[ncount])+Corrfun(qrad,v,mu,zeta,taus[ncount])
		
		# plotting all three correlation functions (i.e. one of the temperature points of Figure 3a)
		plt.figure()
		plt.loglog(qrad,Gofq/v**2,'.-',color='r', label = 'modes')
		plt.loglog(qrad,Fouriervel/v**2,'o',color='k',label='simulation')
		plt.loglog(qrad,Gofq_analytics/v**2,'.--',color='b',label='analytics')
		plt.xlabel('q')
		plt.ylabel('Fourier/v2')
		plt.title('Fourier correlations, tau = ' + str(tau))
		
		
		# saving pickle file for final plots
		data={'tempval': tempval, 
          'nuval': nuval,
          'Nsnap': sim.Nsnap,
          'tplot': tplot,
          'msd': msd,
          'qval': qval,
          'SelfInt': SelfInt, 
          'qrad': qrad,
          'Gofq': Gofq,
          'Gofq_analytics': Gofq_analytics,
          'Fouriervel': Fouriervel,
          'B': B,
          'mu': mu}
		outpickle=fileout +"Soft_temp_" + temp + "_nu_" + nu + ".p"
		pickle.dump(data,open(outpickle,'wb'))			
	ncount+=1		
		
plt.show()
		











