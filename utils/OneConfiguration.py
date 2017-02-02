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
