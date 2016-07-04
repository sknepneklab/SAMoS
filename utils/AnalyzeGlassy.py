# * *************************************************************
# *  
# *   Soft Active Mater on Surfaces (SAMoS)
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
from Glassy import *
from Writer import *

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="input file (base name) ")
parser.add_argument("-r", "--radii", type=str, help="radii file (initial configuration) ")
parser.add_argument("-c", "--conffile", type=str, help="configuration file")
parser.add_argument("-d", "--directory", type=str, help="input directory")
parser.add_argument("-o", "--output", type=str, help="output directory")
parser.add_argument("-s", "--skip", type=int, default=0, help="skip this many samples")
parser.add_argument("--tracer", action='store_true', default=False, help="Is this a configuration with tracer particles")
parser.add_argument("--getMSD", action='store_true', default=False, help="Compute mean squared displacement")
parser.add_argument("--getVelcorr", action='store_true', default=False, help="Compute velocity correlation function")
parser.add_argument("--plot", action='store_true', default=False, help="Plot MSD and correlations")
parser.add_argument("--ignore", action='store_true', default=False, help="Ignore complications for quick result (warning!)")

args = parser.parse_args()
sim = SimRun(args.directory,args.conffile,args.input,args.radii,args.skip,args.ignore,args.tracer,args.plot)
output=False
if args.getMSD:
	tplot,msd = sim.getMSD()
	data={'tplot':tplot,'msd':msd,'configuration':args.conffile}
	output=True
if args.getVelcorr:
	bins,velcorr,fig=sim.getVelcorr(0.5)
	print bins
	print velcorr
	data={'bins':bins,'velcorr':velcorr,'configuration':args.conffile}
	output=True
if output:
	outglassy=args.output + '/glassy.p'
	pickle.dump(data,open(outglassy,'wb'))
if args.plot:
	plt.show()

	
	
