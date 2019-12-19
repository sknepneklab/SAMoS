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

sys.path.insert(1,'/home/sh18581/Documents/Cornea/Analysis/SAMoS/utils/')

from Geometry import *
from Cornea import *
from Writer import *

# WTF is wrong with my plotting??
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()
parser.add_argument("-c", "--conffile", type=str, help="configuration file")
parser.add_argument("-d", "--directory", type=str, help="input directory")
parser.add_argument("-o", "--output", type=str, help="output directory")
parser.add_argument("-p", "--prefix", type=str, default="cornea",help="prefix for output file")
parser.add_argument("-s", "--skip", type=int, default=3000, help="skip this many samples")
parser.add_argument("-m", "--howmany", type=int, default=30, help="read this many samples")
parser.add_argument("-t", "--step", type=int, default=1, help="step snapshots with this spacing in flow field")
parser.add_argument("-u", "--maxtype", type=int, default=3, help="Up to what maximum type should I process data?")
#parser.add_argument("--getMSD", action='store_true', default=False, help="Compute mean squared displacement?")
parser.add_argument("--plot", action='store_true', default=False, help="Plot MSD and correlations")
parser.add_argument("--ignore", action='store_true', default=False, help="Ignore complications like missing potentials for quick result (warning!)")

args = parser.parse_args()
#def __init__(self,directory,conffile,skip,howmany,ignore=True,maxtype=3):
mycornea = Cornea(args.directory,args.conffile,args.skip,args.howmany,args.ignore,args.maxtype)
data={'configuration':args.conffile}
#if args.getMSD:
  ##getMSD(self,verbose=True):
	#tplot,msd = sim.getMSD(args.plot)
	#dataMSD={'tplot':tplot,'msd':msd,}
	#data.update(dataMSD)
mycornea.getFlowField(1,args.step)


# Save the output files
outcornea=args.output + args.prefix +'.p'
pickle.dump(data,open(outcornea,'wb'))
plt.show()

	

	
	
