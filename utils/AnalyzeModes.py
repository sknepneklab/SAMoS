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


	
	
