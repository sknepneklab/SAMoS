#!/usr/bin/env python

import numpy as np
import ioutils as io
#from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib import pyplot as plt
import time
import glob

# is building a class on top of matplotlib which manages the graph formatting 
#overengineering the problem

# i.e.
#class Extfigure(object):
    #def __init__(self, ddump):
        #self.ddump = ddump

# we are using numpys .npz now
def plot_radial(fils):
    fig = plt.figure()
    axi = fig.add_subplot(111)
    #axi.show()
    axi = plt
    for fi in fils:
        data = np.load(fi)
        xs = data['rspace']
        assert len(xs) > 1
        xdiff = xs[1] - xs[0]
        x = xs[:-1] + xdiff/2
        rad = data['radial_pressure']
        axi.plot(x, rad, marker='o')
        #axi.draw()
        axi.show()
        time.sleep(0.05)
        
def readdotplot(fi):
        ddump = io.readdump(fi, firstc=float)
        return ddump

def plot_pressures(ddump):
    xtitle, xls =  ddump.items()[0]
    fig = mpl.figure()
    axi = fig.add_subplot(111)
    axi.xlabel = xtitle
    for ytitle, yls in ddump.items()[1:]:
        axi.plot(xls, yls, label=ytitle)
    axi.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    axi.show()

if __name__=='__main__':

    import sys
    import argparse

    parser = argparse.ArgumentParser()
 
    parser.add_argument("-i", "--input", type=str, default='*.npz', help="Input dat file")
    #parser.add_argument("-d", "--dir", type=str, default='/home/dan/tmp/', help="Output directory")
    args = parser.parse_args()

    fnls = sorted(glob.glob(args.input))
    plot_radial(fnls)
    
    #plot_pressures(ddump)
