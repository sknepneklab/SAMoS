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

import argparse

import writemesh as wr
import os.path as path


# Define analysis routines here

# We need one for catching the pressure of some random cells
def random_choice(pv, nr=5):
    # nr is the number of random points to choose
    bulk = pv.get_mesh_bulk(pv.tri)
    bulk_choice = np.random.choice(bulk, size=nr, replace=False)
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
        for i, cid in enumerate(self.choice):
            self.outd[self.xvar].append(xxv)
            self.outd[cid].append(pr[i])
        outv = [xxv] + list(pr)
        outl = self.line.format(*outv)
        self.outf.write(outl + '\n')

# Preferably would like to write out data to a structured format at each step
# h5py and hdf5 are powerful tools for this but I also want a simple solution
# numpy.savez is the simple solution for now.
# Can't seem to decide if I want the number of data points over which the radial
# averaging is calculated to change between time steps
# for now. let it change. store a whole set of files.
class STracker(Tracker):
    def __init__(self, fileo, xvar='time'):
        # track radially averaged profiles of pressure, pressure variance, stress components ...
        self.xvar = xvar
        self.fileo = fileo
        self.qdict = OrderedDict()

    def update(self, pv, xxv):
        fileo = self.fileo
        fileo = '_'.join([fileo, xxv])
        print 'saving radial data to ', fileo+'.npz'

        qdict = self.qdict
        qdict['rspace'] = pv.rspace 
        for rname, rav in zip(pv.rnames, pv.ravg):
            qdict[rname] = rav
        #qdict['radial_pressure'] = pv.radial_pressure
        if len(qdict) != 0:
            np.savez(fileo, **qdict)

def hash_function(f, linspace):
    l = linspace[-1] - linspace[0]
    ne = len(linspace)
    hf = map(f, linspace)
    def nf(r):
        i = int(round((ne-1) * r/l))
        return hf[i]
    return nf

def quartic_wl(wl):
    def quartic(r):
        return 5/(np.pi*wl**2) * (1 + 3*(r/wl))*(1-(r/wl))**3 if r < wl else 0. 
    return quartic

import glob
# The template senario
class Senario(object):
    def __init__(self, args):
        # Use the args object from argsparse and add to it as necessary to keep
        #  track of all the parameters
        self.args= args

        inp_file = args.input
        self.infiles = sorted(glob.glob(inp_file))
        if self.infiles == []:
            print 'Did not find any files using %s', inp_file
            sys.exit()

        self.fid = 0
        self.numf = len(self.infiles)
        self.verb = True
        self.alog = 'analysis.log'

    def __iter__(self):
        return self

    def next(self, step=1):
        if self.fid < self.numf:
            self.dataf = self.infiles[self.fid]
            self.outnum = self._outnum() # must call this function
            if self.verb:
                print 'working on file number ', self.outnum
            try:
                self._operate()
            except:
                self._finishup()
                raise
            self.fid += step
        else:
            self._finishup()
            raise StopIteration()

    def _outnum(self):
        inp_dir, base_file = path.split(self.dataf)
        self.args.inp_dir = inp_dir
        base_name, ext = path.splitext(base_file)
        outnum= base_name.split('_')[-1]
        self.args.outnum = outnum
        return outnum

    def _operate(self):
        pass # the code to execure on each data file

    def _finishup(self):
        pass # the cleanup code that needs to be executed

    def save_parameters(self, paramsf='out.parameters'):
        with open(paramsf, 'w') as fp:
            if self.verb:
                'recording the analysis parameters to ', paramsf
            io.dumpargs(self.args.__dict__, fp)


set_choice = [844, 707, 960, 606, 958, 801, 864, 859, 870, 302, 500]
class Stress_Senario(Senario):
    def __init__(self, args):
        super(Stress_Senario, self).__init__(args)
        print 'initializing the run generator'
        wl= args.wl

        self.fileo=  os.path.join(args.dir, 'stress_st') 

        om = quartic_wl(wl)
        #hspace = np.linspace(0, wl, 10001, endpoint=True)
        #om = hash_function(om, hspace)
        om.wl = wl

        self.omega = om
        self.st = STracker(self.fileo)

    def _operate(self):
        args = self.args
        outnum = self._outnum()
        omega = self.omega
        k, gamma, L  = args.k, args.gamma, args.L

        rdat = ReadData(self.dataf)
        facefile = path.join(args.inp_dir, 'faces_' + outnum + '.fc')
        simp, _ = io.readfc(facefile)
        pv = PVmesh.datbuild(rdat, simp)
            
        # Handle K, and Gamma
        nf = pv.tri.n_vertices()
        # Could easily read these from a .conf file
        K = np.full(nf, k)
        Gamma = np.full(nf, gamma)
        cl = pv._construct_cl_dict(L)
        pv.set_constants(K, Gamma, cl)

        pv.calculate_energy()
        #pv.calculate_forces()
        #tmp test
        cell_outname = 'cell_dual_' + outnum+ '.vtp'
        mout = path.join(args.dir, cell_outname)
        wr.writemeshenergy(pv, mout)

        tri_outname = 'cell_' + outnum + '.vtp'
        tout = path.join(args.dir, tri_outname)
        wr.writetriforce(pv, tout)


        if self.verb: print 'Just calculate for one value of smoothing length'
        pv._stress_setup()
        pv._set_wl(omega)
        clist = set_choice if args.s else None
        pv.stress_on_centres(omega, clist=clist, kinetic = False) 
        #range_wl(args, pv, [wl], ellipses=True)
        pv.radial()
        self.st.update(pv, outnum)

        #track.update(pv, outnum)
        outdir = args.dir

        cell_outname = 'cell_dual_' + outnum+ '.vtp'
        mout = path.join(outdir, cell_outname)
        wr.writemeshenergy(pv, mout)

        tri_outname = 'cell_' + outnum + '.vtp'
        tout = path.join(outdir, tri_outname)
        wr.writetriforce(pv, tout)

        stress_outname = 'hardy_stress_' + outnum + '.vtp'
        sout = path.join(outdir, stress_outname)
        wr.write_stress_ellipses(pv, sout, pv.stress)

# define a different function to do 
ne = 1001 # The number of times the smoothing function is evaluated
wrange = np.arange(0.1, 4, 0.1)
def range_wl(args, pv, wls, ellipses=True):
    #choice = choices(pv, nr=8)
    choice = args.selection
    #print 'choice',  choice
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

    stressrun = Stress_Senario(args)
    stressrun.save_parameters()
    for _ in stressrun: pass


#    sys.exit()

    #inp_file = args.input
    #import glob
    #infiles = sorted(glob.glob(inp_file))
    #if infiles == []:
        #print 'Did not find any files using %s', inp_file
        #sys.exit()

    ## Cell parameters
    #k = 1.
    #gamma = 0.
    #L = 0.
    #wl = args.wl

    #trackers = []
    #def initial_setup(args, pv):
        #fileo=  os.path.join(args.dir, 'pressure_t.plot') 
        #ptrack = Pressure_t(set_choice, fileo)
        #return ptrack

    #first = True
    #for fin in infiles:

        
        #inp_dir, base_file = path.split(fin)
        #base_name, ext = path.splitext(base_file)
        ## Naming convention
        #outnum= base_name.split('_')[-1]
        #args.outnum = outnum
        #print 'working on file number ', outnum

        #rdat = ReadData(fin)
        #facefile = 'faces_' + outnum + '.fc'
        #simp, _ = io.readfc(facefile)
        #pv = PVmesh.datbuild(rdat, simp)
        #if first:
            #ptrack = initial_setup(args, pv)
            #first = False
            
        #if args.s:
            #args.selection=set_choice
            #print 'using choice'
            #print set_choice
        #else:
            #args.selection=pv.tri_bulk


        ## Handle K, and Gamma
        #nf = pv.tri.n_vertices()
        ## Could easily read these from a .conf file
        #K = np.full(nf, k)
        #Gamma = np.full(nf, gamma)
        #cl = pv._construct_cl_dict(L)
        #pv.set_constants(K, Gamma, cl)

        #pv.calculate_energy()
        ##pv.calculate_forces()

        #if args.w:
            #wls = wrange
            #print 'Going to calculate stress for a range of smoothing lengths'
            #print wls
            #range_wl(args, pv, wls)
        #else:
            #print 'Just calculate for one value of smoothing length'
            #range_wl(args, pv, [wl], ellipses=True)

        #outdir = args.dir

        #cell_outname = 'cell_dual_' + outnum+ '.vtp'
        #mout = path.join(outdir, cell_outname)
        #wr.writemeshenergy(pv, mout)

        #tri_outname = 'cell_' + outnum + '.vtp'
        #tout = path.join(outdir, tri_outname)
        #wr.writetriforce(pv, tout)

        #if False:
            #naive_stress = 'naive_stress_' + outnum + '.vtp'
            #sout = path.join(outdir, naive_stress)
            #wr.write_stress_ellipses(pv, sout, pv.n_stress)

        #if pv.stress:
            #stress_outname = 'hardy_stress_' + outnum + '.vtp'
            #sout = path.join(outdir, stress_outname)
            #wr.write_stress_ellipses(pv, sout, pv.stress)
            ##stddict(pv.stress)


