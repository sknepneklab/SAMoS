#!/usr/bin/env python

# Construct test cases for debugging and later construct simulation initial conditions

from cellmesh import *
import ioutils as io

import numpy as np
import sys, os
import math as m
import random

from collections import OrderedDict
from matplotlib import pyplot as plt

# Define analysis routines here

# We need one for catching the pressure of some random cells
def choices(pv, nr=5):
    # nr is the number of random points to choose
    bulk = pv.get_mesh_bulk(pv.tri)
    #boundary = pv.mesh.get_mesh_boundary(pv.mesh)
    bulk_choice = np.random.choice(bulk, size=nr, replace=False)
    #boundary_choice = np.random.choice(boundary, size=nr, replace=False)
    return bulk_choice


# template for Tracker objects
class Tracker(object):
    def __init__(self,fileo):
        # setup
        self.outf = open(fileo, 'w')
        self.track_id = 2
        self.outd= OrderedDict()
    
    def update(self):
        self.track_id += 1
    
    def write_whole(self):
        io.dump(self.outd, self.outf)

    def cleanup(self):
        self.outf.close()

class Pressure_t(Tracker):
    def __init__(self, choice, fileo, xvar='time'):
        super(Pressure_t, self).__init__(fileo)
        print 'printing to ', fileo
        self.choice = choice
        self.outd[xvar]= []
        for cid in choice:
            self.outd[cid]= []
        # header
        self.outf.write('# ' + xvar + ' ' + '  '.join(map(str, list(choice))) + '\n')
        self.line = '{}  ' + '  '.join(['{}']*len(choice))
        self.xvar = xvar

    def update(self, pv, xxv):
        super(Pressure_t, self).update()
        pr = [pv.pressure[ch] for ch in self.choice]
        #nones = [i for i, x in enumerate(pr) if x == None]
        #for none in reverse(nones):
            #self.choice.pop(none)
        for i, cid in enumerate(self.choice):
            self.outd[self.xvar].append(xxv)
            self.outd[cid].append(pr[i])
        outv = [xxv] + list(pr)
        outl = self.line.format(*outv)
        self.outf.write(outl + '\n')

def hash_function(f, linspace, ne):
    l = linspace[-1] - linspace[0]
    hf = map(f, linspace)
    def nf(r):
        i = int(round((ne-1) * r/l))
        return hf[i]
    return nf

def quartic_wl(wl):
    def quartic(r):
        return 5/(np.pi*wl**2) * (1 + 3*(r/wl))*(1-(r/wl))**3 if r < wl else 0. 
    return quartic

# define a different function to do 
ne = 1001 # The number of times the smoothing function is evaluated
set_choice = [844, 707, 960, 606, 958, 801, 864, 859, 870, 302, 500]
wrange = np.arange(0.1, 4, 0.1)
def range_wl(args, pv, wls, ellipses=True):
    #choice = choices(pv, nr=8)
    choice = args.selection
    print 'chioce',  choice
    fileo=  os.path.join(args.dir, 'pressure_wl.plot') 
    ptrack = Pressure_t(choice, fileo, xvar='wl')
    pv._stress_setup()
    for wl in wls:
        print 'using wl ', wl
        om = quartic_wl(wl)
        lspace = np.linspace(0., wl, ne)
        #omega = hash_function(om, lspace, ne)
        omega = om
        omega.wl = wl
        pv._set_wl(omega)
        pv.stress_on_centres(omega, clist=choice, kinetic = False) 
        
        ptrack.update(pv, wl)

        # change the numbering in this line if you change wrange spacing
        if ellipses:
            wlnum = ('%d' % int(100*wl)).zfill(3)
            stress_outname = 'hardy_stress_' + args.outnum + '_' + wlnum + '.vtp'
            sout = path.join(args.dir, stress_outname)
            wr.write_stress_ellipses(pv, sout, pv.stress)
    ptrack.cleanup()

if __name__=='__main__':

    epidat = '/home/dan/cells/run/soft_rpatch/epithelial_equilini.dat'

    import argparse

    import writemesh as wr
    import os.path as path

    parser = argparse.ArgumentParser()
 
    parser.add_argument("-i", "--input", type=str, default='cell*.dat', help="Input dat file")
    parser.add_argument("-d", "--dir", type=str, default='/home/dan/tmp/', help="Output directory")
    parser.add_argument("-w", action='store_true')
    parser.add_argument("-wl", type=float, default=1., help='smoothing length')
    parser.add_argument("-s", action='store_true')
    #parser.add_argument("-o", "--output", type=str, default=epidat, help="Input dat file")
    args = parser.parse_args()

    if not os.path.exists(args.dir):
        os.mkdir(args.dir)

    inp_file = args.input
    import glob
    infiles = sorted(glob.glob(inp_file))
    if infiles == []:
        print 'Did not find any files using %s', inp_file
        sys.exit()

    # Cell parameters
    k = 1.
    gamma = 0.
    L = 0.
    wl = args.wl

    trackers = []
    def initial_setup(args, pv):
        fileo=  os.path.join(args.dir, 'pressure_t.plot') 
        ptrack = Pressure_t(set_choice, fileo)
        return ptrack

    first = True
    for fin in infiles:

        
        inp_dir, base_file = path.split(fin)
        base_name, ext = path.splitext(base_file)
        # Naming convention
        outnum= base_name.split('_')[-1]
        args.outnum = outnum
        print 'working on file number ', outnum

        rdat = ReadData(fin)
        facefile = 'faces_' + outnum + '.fc'
        simp, _ = io.readfc(facefile)
        pv = PVmesh.datbuild(rdat, simp)
        if first:
            ptrack = initial_setup(args, pv)
            first = False
            
        if args.s:
            args.selection=set_choice
            print 'using choice'
            print set_choice
        else:
            args.selection=pv.tri_bulk


        # Handle K, and Gamma
        nf = pv.tri.n_vertices()
        # Could easily read these from a .conf file
        K = np.full(nf, k)
        Gamma = np.full(nf, gamma)
        cl = pv._construct_cl_dict(L)
        pv.set_constants(K, Gamma, cl)

        pv.calculate_energy()
        #pv.calculate_forces()

        if args.w:
            wls = wrange
            print 'Going to calculate stress for a range of smoothing lengths'
            print wls
            range_wl(args, pv, wls)
        else:
            print 'Just calculate for one value of smoothing length'
            range_wl(args, pv, [wl], ellipses=True)

        ptrack.update(pv, outnum)
        outdir = args.dir

        cell_outname = 'cell_dual_' + outnum+ '.vtp'
        mout = path.join(outdir, cell_outname)
        wr.writemeshenergy(pv, mout)

        tri_outname = 'cell_' + outnum + '.vtp'
        tout = path.join(outdir, tri_outname)
        wr.writetriforce(pv, tout)

        if False:
            naive_stress = 'naive_stress_' + outnum + '.vtp'
            sout = path.join(outdir, naive_stress)
            wr.write_stress_ellipses(pv, sout, pv.n_stress)

        if pv.stress:
            stress_outname = 'hardy_stress_' + outnum + '.vtp'
            sout = path.join(outdir, stress_outname)
            wr.write_stress_ellipses(pv, sout, pv.stress)
            #stddict(pv.stress)


