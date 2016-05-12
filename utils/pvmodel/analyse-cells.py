#!/usr/bin/env python

# Construct test cases for debugging and later construct simulation initial conditions

from cellmesh import *

import numpy as np
import sys, os


if __name__=='__main__':

    epidat = '/home/dan/cells/run/soft_rpatch/epithelial_equilini.dat'

    import argparse

    import writemesh as wr
    import os.path as path

    parser = argparse.ArgumentParser()
 
    parser.add_argument("-i", "--input", type=str, default='cell*.dat', help="Input dat file")
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
        #pv.calculate_forces()

        wl = 1.
        pv._stress_setup()
        pv.stress_on_centres(wl)

        outdir = args.dir

        if pv.forces:
            cell_outname = 'cell_dual_' + outnum+ '.vtp'
            mout = path.join(outdir, cell_outname)
            wr.writemeshenergy(pv, mout)

            tri_outname = 'cell_' + outnum + '.vtp'
            tout = path.join(outdir, tri_outname)
            wr.writetriforce(pv, tout)
            
            naive_stress = 'naive_stress_' + outnum + '.vtp'
            sout = path.join(outdir, naive_stress)
            wr.write_stress_ellipses(pv, sout, pv.n_stress)

        if pv.stress:
            stress_outname = 'hardy_stress_' + outnum + '.vtp'
            sout = path.join(outdir, stress_outname)
            wr.write_stress_ellipses(pv, sout, pv.stress)
            #stddict(pv.stress)


