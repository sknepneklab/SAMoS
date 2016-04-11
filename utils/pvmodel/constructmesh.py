#!/usr/bin/env python

# Construct test cases for debugging and later construct simulation initial conditions

import writemesh as wr

import numpy as np
import sys

from collections import OrderedDict



# Triangulation needs positions and associated areas 
def makeout(rvals, parea):
    outd = OrderedDict()
    keys = ['id', 'x', 'y', 'z', 'area']
    for k in keys:
        outd[k] = []

    for i, r in enumerate(rvals):
        x, y, z = r
        outd['x'].append( x )
        outd['y'].append( y )
        outd['z'].append( z ) 
        outd['id'].append( i )
        outd['area'].append( parea[i] )
    return outd

# Single cell, n sides


def single(nsides, origin=[0.,0.,0.]):
    parea = 5.
    rvals = np.array(origin).reshape((1,3))
    rad = np.linspace(-np.pi, np.pi, nsides, endpoint=False)
    nb = len(rad)
    x = np.cos(rad)
    y = np.sin(rad)
    z = np.zeros(nb)
    area = np.full(nb+1, parea)
    rbound = np.column_stack([x, y, z])
    print rvals
    print rbound
    rvals = np.vstack([rvals, rbound])
    #print np.column_stack([x, y, z])
    print rvals

    outd = makeout(rvals, area)
    return outd
single.defaults = [6]
single.call = ['single( nsides )']

    #wr.datdump

if __name__=='__main__':
    import argparse

    args = sys.argv[1:]
    nargs = len(args)
    if nargs is 0:
        print 'python constructmesh <function> [*arguments]'
        sys.exit()

    fname = args[0]
    fargs = args[1:]
    try:
        f_using = locals()[fname]
    except KeyError:
        print 'I don\'t have a function \'%s\'', fname
        raise

    if len(fargs) is 0:
        print 'No arguments given to constructmesh.%s' % fname
        fargs = f_using.defaults
        print 'Using defaults ', f_using.defaults

    print fargs
    print map(eval, fargs)
    ff = map(eval, fargs)
    outd = f_using(*ff)

    with open('test.dat', 'w') as fo:
        print outd
        wr.datdump(outd, fo)


