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
parser.add_argument("--nematic", action='store_true', default=False, help="Shift n vectors such that particle is in the middle of director.")
parser.add_argument("--contractile", action='store_true',default=False, help="Adds contractile stresses to calculation")
parser.add_argument("-a", "--alpha", type=float, default=0.0, help="Prefactor of the contractile term")
parser.add_argument("--closeHoles", action='store_true', default=False, help="Closes the holes in the tesselation to help tracking of defects (recommended for low density and/or nematic)")
parser.add_argument("--makeEdges",action='store_true', default=False, help="Make edges to the tesselation along borders")
parser.add_argument("--writeP",action='store_true', default=False, help="Output particle positions velocities directors.")
parser.add_argument("--writeT",action='store_true', default=False, help="Output tesselation")
parser.add_argument("--writeD",action='store_true', default=False, help="Output defects")
parser.add_argument("--getStatsBasic",action='store_true',default=False, help="Output basic stats (v2av, packing fraction, energy, pressure)")
args = parser.parse_args()


params = Param(args.directory+args.conffile)

files = sorted(glob(args.directory + args.input+'*.dat'))[args.skip:]

if len(files) == 0:
  files = sorted(glob(args.directory + args.input+'*.dat.gz'))[args.skip:]

if args.writeD:
	defects_n_out=[[] for u in range(len(files))]
	defects_v_out=[[] for u in range(len(files))]
	numdefects_n_out=np.zeros(len(files))
	numdefects_v_out=np.zeros(len(files))
	
if args.getStatsBasic:
	vel2av=np.zeros(len(files))
	phival=np.zeros(len(files))
	pressav=np.zeros(len(files))
	energy=np.zeros(len(files))

u=0
for f in files:
	print f
	conf = Configuration(params,f)
	if args.contractile:
		writeme = Writer(args.nematic,args.alpha)
	else:
		writeme = Writer(args.nematic)
	if args.writeP:
		outparticles = args.output + '/frame%06d_particles.vtp' % u  # + str(u) + '_particles.vtp'
		print outparticles
		writeme.writeConfigurationVTK(conf,outparticles)
	#plt.show()
	if args.getStatsBasic:
		vel2av[u], phival[u],pressav[u],energy[u]= conf.getStatsBasic()
		print "Mean square velocity: " + str(vel2av[u])
		print "Packing fraction: " + str(phival[u])
		print "Mean pressure: " + str(pressav[u])
		print "Energy: " + str(energy[u])
	if args.writeD or args.writeT:
		conf.getTangentBundle()
		tess = Tesselation(conf)
		print "initialized tesselation"
		LoopList,Ival,Jval = tess.findLoop(args.closeHoles)
		print "found loops"
		#print LoopList
		if args.writeD:
			#print "Still to be done ..."
			outdefects = args.output + '/frame%06d_defects.vtp' % u # + str(u) + '_defects.vtp'	
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
			outpatches = args.output + '/frame%06d_patches.vtp' % u #+ str(u) + '_patches.vtp'
			print outpatches
			if args.makeEdges:
				tess.makeEdges(0.85)   
			tess.OrderPatches()
			print "ordered patches"
			if args.contractile:
				writeme.writePatches(tess,outpatches,True)
			else:
				writeme.writePatches(tess,outpatches,True)
	
	u+=1
data={'J':params.J,'v':params.v0,'k':params.pot_params['k'],'pot_params':params.pot_params,'population':params.population,'pop_params':params.pop_params}
if args.writeD:
	dataD={'defects_n':defects_n_out,'defects_v':defects_v_out,'numdefects_n':numdefects_n_out,'numdefects_v':numdefects_v_out}
	data.update(data2)
if args.getStatsBasic:
	dataS={'vel2av':vel2av,'phival':phival,'pressav':pressav,'energy':energy}
	data.update(dataS)
if args.writeD:	
	outpickle=args.directory+'defect_data.p'
else:
	outpickle=args.directory+'configuration_data.p'
print outpickle
pickle.dump(data,open(outpickle,'wb'))
	
	
