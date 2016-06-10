#!/usr/bin/env python

# Construct test cases for debugging and later construct simulation initial conditions

import writemesh as wr

import numpy as np
from numpy.linalg import norm
import sys
import ioutils as io

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

npnormal= np.array([0.,0.,1.])
def fillout(outd, normal=npnormal):
    ld= len(outd.values()[0])
    outd['nvz'] = zeros = np.full(ld, 1.)
    for hh in ['vx', 'vy', 'vz', 'nvx', 'nvy']:
        outd[hh] = np.full(ld, 0.)
    return outd

import collections
#class Pointset(collections.abc.Set):

def _extend_ptset(pts, npts):
    rpts = []
    for npt in npts:
        isnew = True
        for pt in pts:
            if np.allclose(np.array(pt), np.array(npt), rtol=0.1):
                isnew = False
                break
        if isnew:
            rpts.append(npt)
    return pts.extend(rpts)

# minimal hexagonal configuration
def areahex(area=3.0):
    # first find the radius length of a hexagon with the correct area
    l = np.sqrt((2*area)/(3*np.sqrt(3)))
    print l/2.
    # distance between hex centers
    lc = np.cos(np.pi/6) * 2 * l
    print lc
    rads = np.linspace(0, 2*np.pi, 6, False)
    def hexes(rad):
        return [np.cos(rad), np.sin(rad), 0.]
    hpts = lc*np.array(map(hexes, rads))
    pts = list(hpts)
    crads = lc*np.array(map(hexes, rads+(np.pi/6)))
    for hcenter in hpts:
        _extend_ptset(pts, [hcenter+hpt for hpt in hpts])
    print len(pts)
    outd = makeout(pts, [area]*len(pts))
    outd =fillout(outd)
    return outd

from random import random, uniform
def random_direction():
    rtheta = random()* 2*np.pi
    return np.cos(rtheta), np.sin(rtheta)
from read_data import ReadData
# Generic function for operating on the data input files as I read and write them using 
#  separate functions
def rwoperate(ifile, operate):
    rdat = ReadData(ifile)
    keys = rdat.keys
    lrdat = len(rdat.data[0])
    outd = OrderedDict()
    for k in keys:
        outd[k] = []
    for i, _ in enumerate(rdat.data[0]):
        for kn, kv in keys.items():
            outd[kn].append( rdat.data[kv][i] )
        if operate:
            operate(outd, i)
    return outd

# read the dat files and then adjust the preferred areas to a uniform random distrubution
def randomise_area(ifile, minmax=[3.0, 3.5]):
    mn,mx = minmax
    def rarea(outd, i):
        # set the area
        uarea = uniform(mn, mx) 
        #print 'setting area ', uarea
        outd['area'][i] = uarea
    outd = rwoperate(ifile, rarea)

    with open('areaset.input', 'w') as fo:
        io.dump(outd, fo, htag='keys:')

def shift_area(ifile, shift):
    def fshift(outd, i):
        outd['area'][i] += shift
    outd= rwoperate(ifile, fshift)
    with open('shifted.input', 'w') as fo:
        io.dump(outd, fo, htag='keys:')

# hardcoded
def addtype(ifile):
    outd = rwoperate(ifile, None)
    outd['type'] = list(np.zeros(len(outd['id'])))
    print outd['type'][0]
    for i, idd in enumerate(outd['id']):
        if idd < 1000:
            outd['type'][i] = 1
        else:
            outd['type'][i] = 2

    outstr = ' %d\t ' + '%f\t '*(len(outd)-2) + '%d\t ' + '\n'
    with open('typed.input', 'w') as fo:
        io.dump(outd, fo, htag='keys:', outstr=outstr)

# cut a piece of an existing input dat file
#constructmesh.py circleslice epithelial_equilini.dat 
def circleslice(ifile, radius=5.5, boundary=False,  center=[0.,0.,0.]):
    rdat = ReadData(ifile)
    keys = rdat.keys
    x = np.array(rdat.data[keys['x']])
    y = np.array(rdat.data[keys['y']])
    z = np.array(rdat.data[keys['z']])
    rlists = [rdat.data[k] for k in keys.values()]
    outd = OrderedDict()
    for k in keys:
        outd[k] = []
    for i, xval in enumerate(x):
        cd= norm(center - np.array([x[i],y[i],z[i]]))
        if cd < radius:
            # keep this line
            for kn, kv in keys.items():
                outd[kn].append( rdat.data[kv][i] )
    # fix the ids
    ids = range(len(outd.values()[0]))
    lz = ids[-1] +1
    outd['id'] = ids
    # set all the velocities to zero
    for kn in ['vx', 'vy']:
        outd[kn] = list(np.zeros(lz))

    print [len(v) for v in outd.values()]
    #io.stddict(outd)

    # add circular boundary of cells
    if boundary:
        rr = radius + 0.5
        bcelldensity = 0.9
        nbcells = round( (2* np.pi * rr) * bcelldensity)
        thetas = np.linspace(0, 2*np.pi, nbcells, False)
        print 'number of boundary particles', len(thetas)

        xb = rr *np.cos(thetas)
        yb = rr *np.sin(thetas)
        area = rdat.data[keys['area']][0]
        for istepped, xval in enumerate(xb):
            i = lz + istepped
            outd['id'].append(i)
            outd['radius'].append(1.)
            outd['x'].append(xb[istepped])
            outd['y'].append(yb[istepped])
            outd['nvz'].append(1.)
            outd['area'].append(area)
            # make a random direction
            nx, ny = random_direction()
            outd['nx'].append(nx)
            outd['ny'].append(ny)
            for kn in ['z', 'vx', 'vy', 'vz', 'nz', 'nvx', 'nvy']:
                outd[kn].append(0.)
    print map(len, outd.values())
    print outd['id']

    # dump
    with open('slice.input', 'w') as fo:
        io.dump(outd, fo)

# Single cell, n sides
def single(nsides, prefarea=5., origin=[0.,0.,0.]):
    parea = prefarea
    rvals = np.array(origin).reshape((1,3))
    rad = np.linspace(-np.pi, np.pi, nsides, endpoint=False)
    nb = len(rad)
    x = np.cos(rad)
    y = np.sin(rad)
    z = np.zeros(nb)
    area = np.full(nb+1, parea)
    rbound = np.column_stack([x, y, z])
    rvals = np.vstack([rvals, rbound])
    #print np.column_stack([x, y, z])

    outd = makeout(rvals, area)
    outd = fillout(outd)
    return outd
single.defaults = ['6']
single.call = ['single( nsides )']


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

        if hasattr(f_using, 'defaults'):
            fargs = f_using.defaults
            print 'Using defaults ', f_using.defaults
        else:
            fargs = []

    ff = fargs
    print fargs
    for i, f in enumerate(ff):
        try:
            ff[i] = eval(f)
        except:
            # This string has a '.' and look like an object
            #ff[i] = f
            pass
    print ff

    outd = f_using(*ff)
    if outd:
        with open('test.dat', 'w') as fo:
            print outd
            io.datdump(outd, fo)


