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
from Writer import *
from Interaction import *
from Hessian import *

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="input file (base name)")
parser.add_argument("-c", "--conffile", type=str, help="configuration file")
parser.add_argument("-d", "--directory", type=str, help="input directory")
parser.add_argument("-o", "--output", type=str, help="output directory")
parser.add_argument("-n", "--number", type=int, default=1, help="number of file to work on")

args = parser.parse_args()


params = Param(args.directory+args.conffile)
files = sorted(glob(args.directory + args.input+'*.dat'))
nfiles=len(files)
print '----------------------------------------------------------------'
print nfiles
print args.number
if (args.number<= nfiles):
	conf = Configuration(params,files[args.number])
	hess = Hessian(conf,False)
	hess.makeMatrix(True,10.0)
	hess.getModes()
	hess.plotModes(8.0,100)
else:
	print "Error: this configuration doesn't exist!"
	
plt.show()


	
	
