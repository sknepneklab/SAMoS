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
import pickle as pickle
sys.path.append('/home/silke/Documents/CurrentProjects/Rastko/source_git/apcs/utils/') 

from Geometry import *
from Configuration import *
from Tesselation import *
from Defects import *
from Writer import *

basefolder='/home/silke/Documents/CurrentProjects/Rastko/nematic/NEWRUNS/2015-07-04/'
Jlist=['0.1','0.5','1.0','5.0']
vlist=['0.5','1.0','1.5','2.0','2.5','5.0']

writeP=True
writeT=True
writeD=True
outD=True
nematic=True
skip=0

defects_n_out=[]
defects_v_out=[]
numdefects_n_out=[]
numdefects_v_out=[]
for J in Jlist:
	for v in vlist:
		directory=basefolder+'J_' + J + '/v0_' + v +'/'
		conffile='nematic_J_' + J+ '_v0_' + v +'.conf'
		finput='nematic_J_'+J+'_v0_'+v+'_0'
		params = Param(directory+conffile)
		files = sorted(glob(directory + finput+'*.dat'))[skip:]
		u=0
		for f in files:
			conf = Configuration(params,f)
			writeme = Writer(nematic)
			if writeP:
				outparticles = directory + '/frame' + str(u) + '_particles.vtp'
				print outparticles
				writeme.writeConfigurationVTK(conf,outparticles)
			#plt.show()
			conf.getTangentBundle()
			tess = Tesselation(conf)
			print "initialized tesselation"
			LoopList,Ival,Jval = tess.findLoop()
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
				#tess.makeEdges(0.85)   
				tess.OrderPatches()
				print "ordered patches"
				writeme.writePatches(tess,outpatches)
			
			u+=1
		data={'J':J,'v':v,'k':params.pot_params['k'],'defects_n':defects_n_out,'defects_v':defects_v_out,'numdefects_n':numdefects_n_out}
		outpickle=directory+'defect_data.p'
		pickle.dump(data,open(outpickle,'wb'))
			
		







	
	
