
from Geometry import *
from Configuration import *
from Tesselation import *
from Writer import *


basefolder = '/home/sknepnek/DATA/Active/polar/sphere/2014-04-15/J_1.00/v0_5.0/'
conffile = 'sphere_J_1.00_v0_5.0.conf'
filename = 'sphere_J_1.00_v0_5.0_0007450000.dat'
outname = 'out.vtp'

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
writeme.writePatches(tess,outname)
