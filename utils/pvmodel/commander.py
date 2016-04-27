#!/usr/bin/env python

# Construct test cases for debugging and later construct simulation initial conditions

from cellmesh import *

import numpy as np
import sys


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
    inp_file = args.input
    inp_dir, base_file = path.split(inp_file)
    base_name, ext = path.splitext(base_file)

    rdat = ReadData(inp_file)
    PV = PVmesh.datbuild(rdat)

    # Handle K, and Gamma
    nf = PV.tri.n_vertices()
    # Could easily read these from a .conf file
    k = 1.
    gamma = 0.
    K = np.full(nf, k)
    Gamma = np.full(nf, gamma)
    PV.set_constants(K, Gamma)

    PV.calculate_energy()
    PV.calculate_forces()
    PV.calculate_stress()

    outdir = args.dir

    mout = path.join(outdir, 'cellmesh.vtp')
    #print 'saving ', mout
    writemeshenergy(PV, mout)

    tout = path.join(outdir, 'trimesh.vtp')
    #print 'saving ', tout
    writetriforce(PV, tout)

    # Dump the force and energy
    fef = 'force_energy.dat'
    fefo = path.join(outdir, fef)
    PV.out_force_energy(fefo)



