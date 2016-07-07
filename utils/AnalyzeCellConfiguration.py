# * *************************************************************
# *  
# *   Soft Active Mater on Surfaces (SAMoS)
# *   
# *   Author: Rastko Sknepnek
# *  
# *   Division of Physics
# *   School of Engineering, Physics and Mathematics
# *   University of Dundee
# *   
# *   (c) 2013, 2014
# * 
# *   School of Science and Engineering
# *   School of Life Sciences 
# *   University of Dundee
# * 
# *   (c) 2015
# * 
# *   Author: Silke Henkes
# * 
# *   Department of Physics 
# *   Institute for Complex Systems and Mathematical Biology
# *   University of Aberdeen  
# * 
# *   (c) 2014, 2015
# *  
# *   This program cannot be used, copied, or modified without
# *   explicit written permission of the authors.
# * 
# * ***************************************************************

import sys
import argparse
import pickle

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as ptch
import matplotlib.lines as lne
from matplotlib.colors import LinearSegmentedColormap


matplotlib.rcParams['text.usetex'] = 'false'
matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['axes.linewidth'] = 2
matplotlib.rcParams['xtick.major.size'] = 8
matplotlib.rcParams['ytick.major.size'] = 8
matplotlib.rcParams['font.size']=16.0
matplotlib.rcParams['legend.fontsize']=14.0

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
	 

from Geometry import *
from CellConfiguration import *
from Writer import *

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="input file (base name)")
parser.add_argument("-c", "--conffile", type=str, help="configuration file")
parser.add_argument("-d", "--directory", type=str, help="input directory")
parser.add_argument("-o", "--output", type=str, help="output directory")
parser.add_argument("-s", "--skip", type=int, default=0, help="skip this many samples")
parser.add_argument("-z", "--nz", type=int, default=10, help="number of (integer) z bins")
parser.add_argument("-p", "--nprat", type=int, default=100, help="number of p ratio")
parser.add_argument("-m", "--mask", type=int, default=1000, help="radius of the area within which to compute statistics")
parser.add_argument("--prefix",type=str,default='CellStats',help="prefix of pickle output")
parser.add_argument("--verbose",action='store_true', default=False, help="Plot output")
args = parser.parse_args()


params = Param(args.directory+args.conffile)
print params

files = sorted(glob(args.directory + args.input+'*.dat'))[args.skip:]

if len(files) == 0:
  files = sorted(glob(args.directory + args.input+'*.dat.gz'))[args.skip:]

vel2av=np.zeros((len(files),))
f2av=np.zeros((len(files),))
pratav=np.zeros((len(files),))
zav=np.zeros((len(files),))
bfrac=np.zeros((len(files),))
Ninside=np.zeros((len(files),))

ratbin=np.linspace(3.5,4.5,args.nprat+1)
conbin=np.linspace(0,10,args.nz+1)
print conbin
pratdist=np.zeros((len(files),args.nprat))
zdist=np.zeros((len(files),args.nz))

u=0
for f in files:
	print f
	# Ignore here is to simply calculate interactions of multiple types of particles (they have the same potential)
	conf = CellConfiguration(params,f,True)
	#plt.show()#
	vel2av[u], f2av[u], pratav[u],pratdist[u,:],zav[u],zdist[u,:],bfrac[u],Ninside[u]=conf.getStatsCells(ratbin,conbin,args.mask)
	print "Mean square velocity: " + str(vel2av[u])
	print "Mean square force: " + str(f2av[u])
	print "Mean perimeter ratio: " + str(pratav[u])
	print "Contact number: " + str(zav[u])
	u+=1
data={}
#try:
	#data={'v':params.v0,'kappa':params.pot_params['kappa'],'gamma':params.pot_params['gamma'],'lambdaval':params.pot_params['lambdaval'],'population':params.population,'pop_params':params.pop_params}
#except:
	#data={'v':params.v0,'kappa':params.pot_params['kappa'],'gamma':params.pot_params['gamma'],'lambdaval':params.pot_params['lambdaval']}
dataS={'vel2av':vel2av,'f2av':f2av,'pratav':pratav,'ratbin':ratbin,'pratdist':pratdist,'zbin':conbin,'zdist':zdist,'zav':zav,'bfrac':bfrac,'Ninside':Ninside,'mask':args.mask}
data.update(dataS)
if args.prefix=='CellStats':
	outpickle=args.directory+'CellStats.p'
else:
	outpickle=args.directory+args.prefix+'_CellsStats.p'
print outpickle
pickle.dump(data,open(outpickle,'wb'))

if args.verbose:
	pratplot=np.mean(pratdist,axis=0)
	plt.figure()
	plt.plot(ratbin[:-1],pratplot,'.-')
	plt.xlabel('shape parameter')
	plt.ylabel('Probability')
	
	zplot=np.mean(zdist,axis=0)
	plt.figure()
	plt.plot(conbin[:-1],zplot,'.-')
	plt.xlabel('side number')
	plt.ylabel('Probability')
	
	plt.show()
	
	
