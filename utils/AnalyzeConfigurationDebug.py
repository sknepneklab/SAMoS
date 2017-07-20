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
from Tesselation import *
from Defects import *
from Writer import *
from ConfigurationNoParameters import *

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="input file (full name!)")
parser.add_argument("-d", "--directory", type=str, help="input directory")
parser.add_argument("-o", "--output", type=str, help="output directory")
parser.add_argument("-g", "--geometry", type=str, default = "plane_periodic", help="geometry of the system")
parser.add_argument("--lx", type =float, help = "Lx")
parser.add_argument("--ly", type =float, help = "Ly")
parser.add_argument("--lz", type =float, help = "Lz")
parser.add_argument("--delaunay", action='store_true',default=True, help="Use Delaunay triangulation for tesselation.")
parser.add_argument("--nematic", action='store_true', default=False, help="Shift n vectors such that particle is in the middle of director.")
parser.add_argument("--closeHoles", action='store_true', default=False, help="Closes the holes in the tesselation to help tracking of defects (recommended for low density and/or nematic)")
parser.add_argument("--makeEdges",action='store_true', default=False, help="Make edges to the tesselation along borders")
parser.add_argument("-m", "--mult",type = float, default=1.0, help="initial Multiplier for tesselation neighbour radius")
parser.add_argument("--writeP",action='store_true', default=False, help="Output particle positions velocities directors.")
parser.add_argument("--writeT",action='store_true', default=False, help="Output tesselation")
parser.add_argument("--writeD",action='store_true', default=False, help="Output defects")
parser.add_argument("--getStatsBasic",action='store_true',default=False, help="Output basic stats (v2av, packing fraction)")
parser.add_argument("--prefix",type=str,default='frame',help="prefix of vtp output files and pickle output")
args = parser.parse_args()



filename=args.directory  + '/' + args.input
print filename
box = [args.lx, args.ly, args.lz]
conf = Configuration(filename,args.geometry,box)
writeme = Writer(args.nematic)
if args.writeP:
        outparticles = args.output + '/'+ args.prefix + '_particles.vtp' 
        print outparticles
        writeme.writeConfigurationVTK(conf,outparticles)
#plt.show()
if args.getStatsBasic:
        vel2av, phival,ndensity,zav=conf.getStatsBasic()
        #vel2av[u], phival[u],pressav[u],energy[u]= conf.getStatsBasic()
        print "Mean square velocity: " + str(vel2av)
        print "Packing fraction: " + str(phival)
        print "Number density: " + str(ndensity)
        print "Contact number: " + str(zav)
if args.writeD or args.writeT:
        conf.getTangentBundle()
        tess = Tesselation(conf)
        print "initialized tesselation"
        #LoopList,Ival,Jval = tess.findLoop(args.closeHoles,args.mult,1.1)
        LoopList,Ival,Jval = tess.findLoopDelaunay(True)
        print "found loops"
        #print LoopList
        if args.writeD:
                #print "Still to be done ..."
                outdefects = args.output + '/' + args.prefix + '_defects.vtp'
                print outdefects
                defects = Defects(tess,conf)
                # Look for nematic defects in the director field, but do not look for velocity defects (since it's a mess)
                if args.nematic:
                        defects_n, defects_v,numdefect_n,numdefect_v=defects.getDefects('nematic',False)
                else:
                        defects_n, defects_v,numdefect_n,numdefect_v=defects.getDefects('polar')
                print "found defects"
                writeme.writeDefects(defects_n, defects_v,numdefect_n,numdefect_v,outdefects)
        if args.writeT:
                outpatches = args.output + '/' + args.prefix + '_patches.vtp'
                print outpatches
                if args.makeEdges:
                        tess.makeEdges(3.0)   
                tess.OrderPatches()
                print "ordered patches"
                writeme.writePatches(tess,outpatches,True)
                        

	
