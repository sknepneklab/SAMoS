#!/usr/bin/env python

# Construct test cases for debugging and later construct simulation initial conditions

from cellmesh import *

import numpy as np
import sys, os


if __name__=='__main__':

    epidat = '/home/dan/cells/run/soft_rpatch/epithelial_equilini.dat'

    import argparse

    from writemesh import *
    import os.path as path

    parser = argparse.ArgumentParser()
 
    parser.add_argument("-i", "--input", type=str, default=epidat, help="Input dat file")
    parser.add_argument("-d", "--dir", type=str, default='/home/dan/tmp/', help="Output directory")
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

    for fin in infiles:

        
        inp_dir, base_file = path.split(fin)
        base_name, ext = path.splitext(base_file)
        # Naming convention
        outnum= base_name.split('_')[-1]
        print 'working on file number ', outnum

        rdat = ReadData(fin)
        pv = PVmesh.datbuild(rdat)

        # Handle K, and Gamma
        nf = pv.tri.n_vertices()
        # Could easily read these from a .conf file
        k = 1.
        gamma = 0.
        L = 0.
        K = np.full(nf, k)
        Gamma = np.full(nf, gamma)
        cl = pv._construct_cl_dict(L)
        pv.set_constants(K, Gamma, cl)

        pv.calculate_energy()
        print 'Calculating forces and stresses'
        pv.calculate_forces()

        outdir = args.dir

        cell_outname = 'cell_dual_' + outnum+ '.vtp'
        mout = path.join(outdir, cell_outname)
        writemeshenergy(pv, mout)

        tri_outname = 'cell_' + outnum + '.vtp'
        tout = path.join(outdir, tri_outname)
        writetriforce(pv, tout)


        # Dump the force and energy
        #fef = 'force_energy.dat'
        #fefo = path.join(outdir, fef)
        #pv.out_force_energy(fefo)



