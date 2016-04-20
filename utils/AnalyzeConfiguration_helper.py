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

from Geometry import *
from Configuration import *
from Tesselation import *
from Defects import *
from Writer import *

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="input file (base name)")
parser.add_argument("-c", "--conffile", type=str, help="configuration file")
parser.add_argument("-d", "--directory", type=str, help="input directory")
parser.add_argument("-o", "--output", type=str, help="output directory")
parser.add_argument("-s", "--skip", type=int, default=0, help="skip this many samples")
parser.add_argument("--writeD",action='store_true', default=False, help="Output defects")
parser.add_argument("--getStatsBasic",action='store_true',default=False, help="Output basic stats (v2av, packing fraction, energy, pressure)")

args = parser.parse_args()


params = Param(args.directory+args.conffile)
files = sorted(glob(args.directory + args.input+'*.dat'))[args.skip:]
nfiles=len(files)
print '----------------------------------------------------------------'
print nfiles
# Get the relevant runtime parameters
print params.dump['freq']
# Not this one - need to fix, expects a NVE as minimisation step. Need it to take *any* first run as a minimisatin step
#print params.nsteps
print params.dt
# Also not necessary
#print params.int_params
DataRun={'dt':params.dt, 'freqsave':int(params.dump['freq']),'nfiles':nfiles}
print DataRun

if args.writeD:	
	outpickle=args.directory+'defect_data.p'	
else:
	outpickle=args.directory+'configuration_data.p'
print outpickle
# First read the dictionary out of the file
data=pickle.load(open(outpickle, "rb"))
# add our new data to it
data.update(DataRun)
# and save back
pickle.dump(data,open(outpickle,'wb'))
print data

	
	
