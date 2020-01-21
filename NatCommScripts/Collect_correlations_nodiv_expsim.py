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
import sys
import os
import glob
import pickle as pickle
import copy as cp
import numpy as np
import scipy as sp
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as ptch
import matplotlib.lines as lne
from matplotlib.colors import LinearSegmentedColormap
try:
	from Glassy import *
except:
	print "Could not load SAMoS analysis modules."
	print "Please make sure that the shell variable PYTHONPATH"
	print "includes the location of the SAMoS analysis modules."
	print "These modules are typically located in the analysis"
	print "directory inside you SAMoS installation directory."
	sys.exit(1)

matplotlib.rcParams['text.usetex']      = 'false'
matplotlib.rcParams['lines.linewidth']  = 2
matplotlib.rcParams['axes.linewidth']   = 2
matplotlib.rcParams['xtick.major.size'] = 8
matplotlib.rcParams['ytick.major.size'] = 8
matplotlib.rcParams['font.size']        = 16.0
matplotlib.rcParams['legend.fontsize']  = 14.0

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


nuval   = ['0.8']
kval    = ['55']
v0val   = ['90']
phival  = ['0.95']


plotting = True
N = 1434

step = 1
for phi in phival:
  for v0 in v0val:
    for nu in nuval:
      for k in kval:
        confdir = filebase + "nodiv/phi_" + phi + "/nu_" + \
          nu + "/v0_" + v0 + "/k_" + k + "/conf3/"
        radiusfile = "plane.txt"
        prefix = "cornea_cells_00"
        conffile = 'epithelial_nodiv_expsim.conf'
        sim = SimRun(confdir, conffile, prefix, radiusfile, 50, False, False, True)

        qrad, valrad, Sqrad = sim.FourierTransVel(0, 0.25, False)
        Fouriervel = np.zeros((len(qrad),))
        u = 0
        for l in range(0, sim.Nsnap, 10):
          if l % 10 == 0:
            print l
          qrad, valrad, Sqrad = sim.FourierTransVel(l, 0.25, False)
          Fouriervel += Sqrad
          u += 1
        Fouriervel /= u
        if plotting:
          plt.figure()
          plt.loglog(qrad, Fouriervel, '.-')
          plt.xlabel('q')
          plt.ylabel('<|v(q)|^2>')
          plt.title('Computational Fourier velocity, v0 = ' + v0)

        tval_auto, velauto, v2av = sim.getVelAuto(False)
        tval_msd, msd = sim.getMSD(False)
        qval = [0.3, 0.3, 0.0]/np.sqrt(2)
        tval_self, SelfInt = sim.SelfIntermediate(qval, False)

        bins = np.linspace(0, 5, 251)
        bins2 = np.linspace(-5, 5, 501)
        vav, vdist, vdist2 = sim.getVelDist(bins, bins2)

        print "Saving data to pickle file"
        outpickle = fileout + "Velcorr_" + "phi" + phi + \
            "_nu" + nu + "_v0" + v0 + "_k" + k + "_expsim.p"
        print outpickle
        data = {"nu": nu, 
                "v0": v0, 
                "k": k, 
                "qrad": qrad, 
                "Fouriervel": Fouriervel, 
                "msd": msd, 
                "SelfInt": SelfInt, 
                "tval_msd": tval_msd, 
                "tval_self": tval_self,
                "bins": bins, 
                "bins2": bins2, 
                "phi": phi, 
                "tval_auto": tval_auto, 
                "velauto": velauto, 
                "v2av": v2av, 
                "vav": vav, 
                "vdist": vdist, 
                "vdist2": vdist2}
        pickle.dump(data, open(outpickle, 'wb'))

plt.show()
