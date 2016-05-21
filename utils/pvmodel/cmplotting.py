#!/usr/bin/env python

import numpy as np
import ioutils as io
#from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib import pyplot as plt
import time
import glob
import os.path as path

# is building a class on top of matplotlib which manages the graph formatting 
#overengineering the problem?

# i.e.
#class Extfigure(object):
    #def __init__(self, ddump):
        #self.ddump = ddump

# we are using numpys .npz now
def plot_radial(*fils):
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

# we are using numpys .npz now
def plot_radial_all(fglob, *fdirs):
    fig = plt.figure()
    #axi = fig.add_subplot(111)
    npzset = []
    for fd in fdirs:
        npzset.append([])
    npzset.append( sorted(glob.glob(path.join(fd, fglob))) )

    axi = plt
    for i, fi in enumerate(npzset[0]):
        data = np.load(fi)
        xs = data['rspace']
        assert len(xs) > 1
        xdiff = xs[1] - xs[0]
        x = xs[:-1] + xdiff/2
        rad = data['radial_pressure']
        axi.plot(x, rad, marker='o')
        #axi.draw()
        time.sleep(0.05)


def justifynum(num):
    return '{:0>10}'.format(num)
#cmplotting.py compare_radial 4000 stress_st . ../../rA_2.5/pressure_test/ ../../vpotential_only/pressure_test/
def compare_radial(num, fname, *fdirs):
    plt.clf()
    fname = '_'.join([fname, justifynum(num)])  + '.npz'
    paths = [path.join(fd, fname) for fd in fdirs]
    plt.xlabel('radial distance')
    plt.ylabel('2 * pressure')
    for fi in paths:
        data = np.load(fi)
        xs = data['rspace']
        assert len(xs) > 1
        xdiff = xs[1] - xs[0]
        x = xs[:-1] + xdiff/2
        rad = data['radial_pressure']
        plt.plot(x, rad, marker='o', label=fi)
    #http://matplotlib.org/users/legend_guide.html
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)
    plt.show()

#def compare(

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

    # This is a cute command line function that takes the first argument to be the name
    #  of a function defined locally and any other arguments as arguments to the function
    args = sys.argv[1:]
    nargs = len(args)
    if nargs is 0:
        print 'python constructmesh <function> [arguments]'
        sys.exit()

    fname = args[0]
    print 'using function ', fname
    fargs = args[1:]
    try:
        f_using = locals()[fname]
    except KeyError:
        print 'I don\'t have a function \'%s\'', fname
        raise

    if len(fargs) is 0:
        print 'No arguments given to constructmesh.%s' % fname

        if hasattr(f_using, 'defaults'):
            fargs = f_using.defaults
            print 'Using defaults ', f_using.defaults
        else:
            fargs = []

    print fargs
    #print map(eval, fargs)
    #ff = map(eval, fargs)
    ff = fargs
    outd = f_using(*ff)



    #parser = argparse.ArgumentParser()
    #parser.add_argument("-i", "--input", type=str, default='*.npz', help="Input dat file")
    ##parser.add_argument("-d", "--dir", type=str, default='/home/dan/tmp/', help="Output directory")
    #args = parser.parse_args()

    #fnls = sorted(glob.glob(args.input))
    #plot_radial(fnls)
    
    #plot_pressures(ddump)
