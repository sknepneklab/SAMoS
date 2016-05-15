#!/usr/bin/env python

# Construct test cases for debugging and later construct simulation initial conditions

from cellmesh import *
import ioutils as io

import numpy as np
import sys, os
import random

from collections import OrderedDict

# Define analysis routines here

# We need one for catching the pressure of some random cells
def choices(pv, nr=5):
    # nr is the number of random points to choose
    bulk = pv.get_mesh_bulk(pv.tri)
    #boundary = pv.mesh.get_mesh_boundary(pv.mesh)
    bulk_choice = np.random.choice(bulk, size=nr, replace=False)
    #boundary_choice = np.random.choice(boundary, size=nr, replace=False)
    return bulk_choice

# template
class Tracker(object):
    def __init__(self,fileo):
        # setup
        self.outf = open(fileo, 'w')
        self.track_id = 0
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


# define a different function to do 
def range_wl(args, pv, wls):
    choice = choices(pv, nr=8)
    fileo=  os.path.join(args.dir, 'pressure_wl.plot')
    ptrack = Pressure_t(choice, fileo)
    pv._stress_setup(wls[0])
    for wl in wls:
        print 'using wl ', wl
        pv._set_wl(wl)
        pv.stress_on_centres(wl, clist=choice)
        ptrack.update(pv, wl)
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
    wl = 1.

    trackers = []
    def initial_setup(pv):
        choice = choices(pv)
        ptrack = Pressure_t(choice)
        trackers.append(ptrack)

    first = True
    for fin in infiles:

        
        inp_dir, base_file = path.split(fin)
        base_name, ext = path.splitext(base_file)
        # Naming convention
        outnum= base_name.split('_')[-1]
        print 'working on file number ', outnum

        rdat = ReadData(fin)
        facefile = 'faces_' + outnum + '.fc'
        simp, _ = io.readfc(facefile)
        pv = PVmesh.datbuild(rdat, simp)
        #if first:
            #initial_setup(pv)
            #first = False

        # Handle K, and Gamma
        nf = pv.tri.n_vertices()
        # Could easily read these from a .conf file
        K = np.full(nf, k)
        Gamma = np.full(nf, gamma)
        cl = pv._construct_cl_dict(L)
        pv.set_constants(K, Gamma, cl)

        pv.calculate_energy()
        #pv.calculate_forces()

        #pv.stress_on_centres(wl)
        if args.w:
            print 'w'
            wls = np.arange(0.1, 5, 0.1)
            print wls
            range_wl(args, pv, wls)

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


