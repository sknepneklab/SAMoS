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
parser.add_argument("--nematic", action='store_true', default=False, help="Shift n vectors such that particle is in the middle of director.")
parser.add_argument("--closeHoles", action='store_true', default=False, help="Closes the holes in the tesselation to help tracking of defects (recommended for low density and/or nematic)")
parser.add_argument("--makeEdges",action='store_true', default=False, help="Make edges to the tesselation along borders")
parser.add_argument("--writeP",action='store_true', default=False, help="Output particle positions velocities directors.")
parser.add_argument("--writeT",action='store_true', default=False, help="Output tesselation")
parser.add_argument("--writeD",action='store_true', default=False, help="Output defects")
args = parser.parse_args()


params = Param(args.directory+args.conffile)

files = sorted(glob(args.directory + args.input+'*.dat'))[args.skip:]

defects_n_out=[[] for u in range(len(files))]
defects_v_out=[[] for u in range(len(files))]
numdefects_n_out=np.zeros(len(files))
numdefects_v_out=np.zeros(len(files))

u=0
for f in files:
	print f
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
		LoopList,Ival,Jval = tess.findLoop(args.closeHoles)
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
			defects_n_out[u]=defects_n
			defects_v_out[u]=defects_v
			numdefects_n_out[u]=numdefect_n
			numdefects_v_out[u]=numdefect_v
			#defects.PlotDefects()
			writeme.writeDefects(defects_n, defects_v,numdefect_n,numdefect_v,outdefects)
		if args.writeT:
			outpatches = args.output + '/frame' + str(u) + '_patches.vtp'
			print outpatches
			if args.makeEdges:
				tess.makeEdges(0.85)   
			tess.OrderPatches()
			print "ordered patches"
			writeme.writePatches(tess,outpatches)
	
	u+=1
data={'J':params.J,'v':params.v0,'k':params.pot_params['k'],'defects_n':defects_n_out,'defects_v':defects_v_out,'numdefects_n':numdefects_n_out,'numdefects_v':numdefects_v_out}
outpickle=args.directory+'defect_data.p'
pickle.dump(data,open(outpickle,'wb'))
	
	
