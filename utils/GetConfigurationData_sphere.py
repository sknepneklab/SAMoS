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
import pickle as pickle
#sys.path.append('/home/silke/Documents/CurrentProjects/Rastko/source_git/apcs/utils/') 

from Geometry import *
from Configuration import *
from Tesselation import *
from Defects import *
from Writer import *

parser = argparse.ArgumentParser()
parser.add_argument("-b", "--basefolder", type=str, help="Base folder with all data files")
parser.add_argument("-s", "--skip", type=int, default=0, help="skip this many samples")
parser.add_argument("--writeP", action='store_true', default=True, help="write particle positions velocities directors vtp file")
parser.add_argument("--writeT", action='store_true', default=True, help="write Tesselation vtp file")
parser.add_argument("--writeD", action='store_true', default=True, help="write Defects vtp file")
parser.add_argument("--nematic", action='store_true', default=True, help="Assume nematic.")
parser.add_argument("--closeHoles", action='store_true', default=True, help="Closes the holes in the tesselation to help tracking of defects (recommended for low density and/or nematic)")
parser.add_argument("--makeEdges",action='store_true', default=False, help="Make edges to the tesselation along borders")
parser.add_argument("--name_1", type=str,default='J', help="Name of the first quantity to loop over")
parser.add_argument("--name_2", type=str,default='v0', help="Name of the second quantity to loop over")
parser.add_argument("--val_1", type=str, action='append', help="Values of quantity one to loop over")
parser.add_argument("--val_2", type=str, action='append', help="Values of quantity two to loop over")
args = parser.parse_args()



basefolder=args.basefolder

Jlist=args.val_1
vlist=args.val_2

writeP=args.writeP
writeT=args.writeT
writeD=args.writeD
nematic=args.nematic
skip=args.skip

defects_n_out=[]
defects_v_out=[]
numdefects_n_out=[]
numdefects_v_out=[]
for J in Jlist:
	for v in vlist:
		directory=basefolder+'/' + args.name_1 +'_' + J + '/'+args.name_2 + '_' + v +'/'
		conffile='nematic_'+args.name_1 + '_' + J + '_'+args.name_2+'_' + v +'.conf'
		finput='nematic_'+args.name_1+'_'+J+'_'+args.name_2+'_'+v+'_0'
		print finput
		params = Param(directory+conffile)
		files = sorted(glob(directory + finput+'*.dat'))[skip:]
		u=0
		for f in files:
			conf = Configuration(params,f)
			writeme = Writer(nematic)
			if writeP:
				outparticles = directory + 'frame' + str(u) + '_particles.vtp'
				print outparticles
				writeme.writeConfigurationVTK(conf,outparticles)
			#plt.show()
			conf.getTangentBundle()
			tess = Tesselation(conf)
			print "initialized tesselation"
			LoopList,Ival,Jval = tess.findLoop(args.closeHoles)
			print "found loops"
			#print LoopList
			if writeD:
				#print "Still to be done ..."
				outdefects = directory + '/frame' + str(u) + '_defects.vtp'	
				print outdefects
				defects = Defects(tess,conf)
				if nematic:
					defects_n, defects_v,numdefect_n,numdefect_v=defects.getDefects('nematic')
				else:
					defects_n, defects_v,numdefect_n,numdefect_v=defects.getDefects('polar')
				defects_n_out.append(defects_n)
				defects_v_out.append(defects_v)
				numdefects_n_out.append(numdefect_n)
				numdefects_v_out.append(numdefect_v)
				print "found defects"
				#defects.PlotDefects()
				writeme.writeDefects(defects_n, defects_v,numdefect_n,numdefect_v,outdefects)
			if writeT:
				outpatches = directory + '/frame' + str(u) + '_patches.vtp'
				print outpatches
				if args.makeEdges:
					tess.makeEdges(0.85)   
				tess.OrderPatches()
				print "ordered patches"
				writeme.writePatches(tess,outpatches)
			
			u+=1
		data={args.name_1:J,args.name_2:v,'k':params.pot_params['k'],'defects_n':defects_n_out,'defects_v':defects_v_out,'numdefects_n':numdefects_n_out}
		outpickle=directory+'defect_data.p'
		pickle.dump(data,open(outpickle,'wb'))
			
		







	
	
