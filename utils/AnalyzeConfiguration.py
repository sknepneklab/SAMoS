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
from Writer import *

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="input file (base name)")
parser.add_argument("-c", "--conffile", type=str, help="configuration file")
parser.add_argument("-d", "--directory", type=str, help="input directory")
parser.add_argument("-o", "--output", type=str, help="output directory")
parser.add_argument("-s", "--skip", type=int, default=0, help="skip this many samples")
parser.add_argument("--nematic", action='store_true', default=False, help="Shift n vectors such that particle is in the middle of director.")
parser.add_argument("-P", "--writeP",action='store_true', default=True, help="Output particle positions velocities directors.")
parser.add_argument("-T", "--writeT",action='store_true', default=True, help="Output tesselation")
parser.add_argument("-D", "--writeD",action='store_true', default=False, help="Output defects")
args = parser.parse_args()


basefolder = '/home/silke/Documents/CurrentProjects/Rastko/Cells/Inke/'
conffile = 'goblet.conf'


filename = 'test3/goblet_test_0001100000.dat'
outname = 'testing'

params = Param(args.directory+conffile)

files = sorted(glob(args.directory + args.input+'*.dat'))[args.skip:]
    
u=0
for f in files:
	conf = Configuration(params,f)
	writeme = Writer()
	if args.writeT:
		outparticles = args.output + '/frame' + str(u) + '_particles.vtp'
		print outparticles
		writeme.writeConfigurationVTK(conf,outparticles)
	#plt.show()
	conf.getTangentBundle()
	tess = Tesselation(conf)
	#print "initialized tesselation"
	LoopList,Ival,Jval = tess.findLoop()
	#print "found loops"
	tess.OrderPatches()
	#print "ordered patches"
	if args.writeT:
		outpatches = args.output + '/frame' + str(u) + '_patches.vtp'
		print outpatches
		writeme.writePatches(tess,outpatches)
	if args.writeD:
		print "Still to be done ..."
		outdefects = basefolder + '/' + outname + '_defects.vtp'	
	u+=1
	
	
