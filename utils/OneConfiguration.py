
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
