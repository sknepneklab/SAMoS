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


from Geometry import *
from Configuration import *
from Tesselation import *
from Writer import *


basefolder = '/home/silke/Documents/CurrentProjects/Rastko/Cells/Inke/'
conffile = 'goblet.conf'
filename = 'test3/goblet_test_0001100000.dat'
outname = 'testing'

params = Param(basefolder+conffile)
conf = Configuration(params,basefolder+filename)
plt.show()
conf.getTangentBundle()
tess = Tesselation(conf)
print "initialized tesselation"
LoopList,Ival,Jval = tess.findLoop()
print "found loops"
tess.OrderPatches()
print "ordered patches"
writeme = Writer()
outparticles = basefolder + '/' + outname + '_particles.vtp'
outdefects = basefolder + '/' + outname + '_defects.vtp'
outpatches = basefolder + '/' + outname + '_patches.vtp'
writeme.writePatches(tess,outpatches)
writeme.writeConfigurationVTK(conf,outparticles)
