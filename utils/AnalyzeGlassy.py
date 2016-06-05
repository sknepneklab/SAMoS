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
#parser.add_argument("--nematic", action='store_true', default=False, help="Shift n vectors such that particle is in the middle of director.")
#parser.add_argument("--closeHoles", action='store_true', default=False, help="Closes the holes in the tesselation to help tracking of defects (recommended for low density and/or nematic)")
#parser.add_argument("--makeEdges",action='store_true', default=False, help="Make edges to the tesselation along borders")
#parser.add_argument("--writeP",action='store_true', default=False, help="Output particle positions velocities directors.")
#parser.add_argument("--writeT",action='store_true', default=False, help="Output tesselation")
#parser.add_argument("--writeD",action='store_true', default=False, help="Output defects")
args = parser.parse_args()
#sim = SimRun(True,args.directory,args.conffile,args.input,args.skip,True)
sim = SimRun(args.directory,args.conffile,args.input,args.radii,args.skip,True,True)
tplot,msd = sim.getMSD()
bins,velcorr,fig=sim.getVelcorr(0.5)
data={'bins':bins,'velcorr':velcorr,'tplot':tplot,'msd':msd,'configuration':args.conffile}
outglassy=args.output + '/glassy.p'
outfig=args.output + '/velcorr.pdf'
pickle.dump(data,open(outglassy,'wb'))
plt.savefig(outfig)
plt.show()


	
	
