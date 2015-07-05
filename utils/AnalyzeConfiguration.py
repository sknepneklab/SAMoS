#    
#    Author: Silke Henkes
#   
#    ICSMB, Department of Physics
#    University of Aberdeen
#    
#    (c)  2015
#    
#    This program cannot be used, copied, or modified without
#    explicit permission of the author.
#  
# ################################################################

import sys
import argparse

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
parser.add_argument("--nematic", action='store_true', default=False, help="Shift n vectors such that particle is in the middle of director.")
parser.add_argument("--writeP",action='store_true', default=False, help="Output particle positions velocities directors.")
parser.add_argument("--writeT",action='store_true', default=False, help="Output tesselation")
parser.add_argument("--writeD",action='store_true', default=False, help="Output defects")
args = parser.parse_args()


params = Param(args.directory+args.conffile)

files = sorted(glob(args.directory + args.input+'*.dat'))[args.skip:]
    
u=0
for f in files:
	conf = Configuration(params,f)
	writeme = Writer(args.nematic)
	if args.writeP:
		outparticles = args.output + '/frame' + str(u) + '_particles.vtp'
		print outparticles
		writeme.writeConfigurationVTK(conf,outparticles)
	#plt.show()
	if args.writeD or args.writeT:
		conf.getTangentBundle()
		tess = Tesselation(conf)
		print "initialized tesselation"
		LoopList,Ival,Jval = tess.findLoop()
		print "found loops"
		#print LoopList
		if args.writeD:
			#print "Still to be done ..."
			outdefects = args.output + '/frame' + str(u) + '_defects.vtp'	
			print outdefects
			defects = Defects(tess,conf)
			if args.nematic:
				defects_n, defects_v,numdefect_n,numdefect_v=defects.getDefects('nematic')
			else:
				defects_n, defects_v,numdefect_n,numdefect_v=defects.getDefects('polar')
			print "found defects"
			#defects.PlotDefects()
			writeme.writeDefects(defects_n, defects_v,numdefect_n,numdefect_v,outdefects)
		if args.writeT:
			outpatches = args.output + '/frame' + str(u) + '_patches.vtp'
			print outpatches
			#tess.makeEdges(0.85)   
			tess.OrderPatches()
			print "ordered patches"
			writeme.writePatches(tess,outpatches)
	
	u+=1
	
	
