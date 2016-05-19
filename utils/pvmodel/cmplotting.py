#!/usr/bin/env python

import numpy as np
import ioutils as io
#from matplotlib import pyplot as plt
import matplotlib as mpl

# is building a class on top of matplotlib which manages the graph formatting 
#overengineering the problem

# i.e.
#class Extfigure(object):
    #def __init__(self, ddump):
        #self.ddump = ddump

def plot_pressures(ddump):
    xtitle, xls =  ddump.items()[0]
    fig = mpl.figure()
    plt = fig.add_subplot(111)
    plt.xlabel = xtitle
    for ytitle, yls in ddump.items()[1:]:
        plt.plot(xls, yls, label=ytitle)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.show()

if __name__=='__main__':

    import sys
    args = sys.argv[1:]
    filein = args[0]
    print 'plotting from file ', filein
    with open(filein, 'r') as fi:
        ddump = io.readdump(fi, firstc=float)
        plot_pressures(ddump)
