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
parser.add_argument("-f", "--inface", type=str, help="input file (base name - faces)")
parser.add_argument("-c", "--conffile", type=str, help="configuration file")
parser.add_argument("-d", "--directory", type=str, help="input directory")
parser.add_argument("-o", "--output", type=str, help="output directory")
parser.add_argument("-s", "--skip", type=int, default=0, help="skip this many samples")
parser.add_argument("-z", "--nz", type=int, default=10, help="number of (integer) z bins")
parser.add_argument("-p", "--nprat", type=int, default=200, help="number of p ratio")
parser.add_argument("-m", "--mask", type=int, default=1000, help="radius of the area within which to compute statistics")
parser.add_argument("--prefix",type=str,default='CellStats',help="prefix of pickle output")
parser.add_argument("--verbose",action='store_true', default=False, help="Plot output")
args = parser.parse_args()


params = Param(args.directory+args.conffile)
print params

files_cells = sorted(glob(args.directory + args.input+'*.dat'))[args.skip:]
files_faces = sorted(glob(args.directory + args.inface+'*.fc'))[args.skip:]

if len(files_cells) == 0:
  files_cells = sorted(glob(args.directory + args.input+'*.dat.gz'))[args.skip:]
  files_faces = sorted(glob(args.directory + args.inface+'*.fc.gz'))[args.skip:]

vel2av = np.zeros((len(files_cells),))
f2av = np.zeros((len(files_cells),))
areav = np.zeros((len(files_cells),))
pratav = np.zeros((len(files_cells),))
zav = np.zeros((len(files_cells),))
borderlen = np.zeros((len(files_cells),))
bfrac = np.zeros((len(files_cells),))
Ninside = np.zeros((len(files_cells),))

areabin = np.linspace(1,5,args.nprat+1)
ratbin = np.linspace(3,5,args.nprat+1)
conbin = np.linspace(0,10,args.nz+1)

areadist = np.zeros((len(files_cells),args.nprat))
pratdist = np.zeros((len(files_cells),args.nprat))
zdist = np.zeros((len(files_cells),args.nz))

u = 0
for u in range(len(files_cells)):
  fcells=files_cells[u]
  ffaces=files_faces[u]
	# Ignore here is to simply calculate interactions of multiple types of particles (they have the same potential)
	conf = CellConfiguration(params,fcells,ffaces,True,False)
	#plt.show()#
	vel2av[u], f2av[u], areav[u], areadist[u,:], pratav[u],pratdist[u,:],zav[u],zdist[u,:],borderlen[u],bfrac[u],Ninside[u]=conf.getStatsCells(areabin,ratbin,conbin,args.mask)
	print "Mean square velocity: " + str(vel2av[u])
	print "Mean square force: " + str(f2av[u])
	print "Mean area: " + str(areav[u])
	print "Total perimeter length: " + str(borderlen[u])
	print "Mean perimeter ratio: " + str(pratav[u])
	print "Contact number: " + str(zav[u])
	u+=1
data = {}

dataS = {'vel2av':vel2av,'f2av':f2av,'areav':areav,'areadist':areadist,'pratav':pratav,'ratbin':ratbin,'pratdist':pratdist,'zbin':conbin,'zdist':zdist,'zav':zav,'borderlen':borderlen,'bfrac':bfrac,'Ninside':Ninside,'mask':args.mask}
data.update(dataS)
if args.prefix == 'CellStats':
	outpickle = args.output+'CellStats.p'
else:
	outpickle = args.output+args.prefix+'_CellsStats.p'
print 'Using pickle file : ', outpickle
pickle.dump(data,open(outpickle,'wb'))

if args.verbose:
  areaplot=np.mean(areadist,axis=0)
	plt.figure()
	plt.plot(areabin[:-1],areaplot,'.-')
	plt.xlabel('cell area')
	plt.ylabel('Probability')
	
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
	
	plt.figure()
	plt.plot(borderlen,'.-')
	plt.xlabel('time')
	plt.ylabel('Border length')
	
	plt.show()
	
	
