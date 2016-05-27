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
        # prepend to 'rspace'
        for stress in pv.stresses.values():
            space_name = '_'.join([stress.name, 'rspace'])
            qdict[space_name] = stress.rspace 
            for rname, rav in stress.ravg.items():
                qdict[rname] = rav

        if len(qdict) != 0:
            print 'saving'
            print fileo
            #print qdict

        np.savez(fileo, **qdict)

import pickle
# Seems like pickling the data might be better
# could pickle pv.stresses
# what other data is useful to pickle?
class Pickler(Tracker):
    def __init__(self, fileo='stresses'):
        self.fileo = fileo

    def update(self, pv, xxv):
        fileo = self.fileo
        fileo = '_'.join([fileo, xxv]) + '.pkl'
        print 'saving stress data to ', fileo

        with open(fileo, 'wb') as fo:
            pickle.dump(pv.stresses, fo)
        


def hash_function(f, linspace):
    l = linspace[-1] - linspace[0]
    ne = len(linspace)
    hf = map(f, linspace)
    def nf(r):
        i = int(round((ne-1) * r/l))
        return hf[i]
    return nf

from scipy import integrate
def quartic_wl(wl):
    print 'generating smoothing function with wl = ', wl
    def quartic(r):
        return 5/(np.pi*wl**2) * (1 + 3*(r/wl))*(1-(r/wl))**3 if r < wl else 0. 
        #return 5/(4* np.pi*wl) * (1 + 3*(r/wl))*(1-(r/wl))**3 if r < wl else 0. 
    space = np.linspace(0, wl, 100, True)
    qq = np.array(map(quartic, space))
    ii = integrate.simps(qq* space*2*np.pi, space)
    print 'integrated, ii = ', ii
    #io.plotrange(quartic, 0, wl)
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

    def next(self):
        step = self.args.step
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

    def _name_vtp(self, name, i=None):
        args= self.args
        outname = name + args.outnum
        if i:
             outname += '_' + str(i) 
        outname += '.vtp'
        return path.join(args.dir, outname)


set_choice = [844, 707, 960, 606, 958, 801, 864, 859, 870, 302, 500]
class Stress_Senario(Senario):
    def __init__(self, args):
        super(Stress_Senario, self).__init__(args)
        print 'initializing the run generator'
        wl= args.wl

        self.fileo=  os.path.join(args.dir, 'stress_st') 

        self.omega= None
        om = quartic_wl(wl)
        self.omega = om

        self.st = STracker(self.fileo)
        wlf = os.path.join(args.dir, 'stress_wl')
        self.wlst = STracker(wlf, xvar='wl')
        self.pkr = Pickler(os.path.join(args.dir, 'stresses'))

    def _read_pv(self):
        rdat = ReadData(self.dataf)
        facefile = path.join(args.inp_dir, 'faces_' + args.outnum + '.fc')
        simp, _ = io.readfc(facefile)
        pv = PVmesh.datbuild(rdat, simp)
        #pv = PVmesh.datbuild(rdat)
        return pv

    def _operate(self):
        args = self.args
        outnum = self._outnum()
        omega = self.omega

        pv = self._read_pv()
        self._handle_constants(pv)
        pv.calculate_energy()
        pv.calculate_forces(exclude_boundary=args.exclude)
        # ready to output simple stress immediately
        nsout = self._name_vtp('simple_stress_')
        wr.write_stress_ellipses(pv, nsout, pv.stresses['simple'])
        if args.simple: return

        if self.verb: print 'Just use one value of smoothing length'
        pv._stress_setup()
        self._iterate_wl(pv)
        pv._set_wl(args.wl)
        clist = set_choice if args.s else None

        pv.stress_on_centres(omega, clist=clist, hardy=args.hardy)
        self.st.update(pv, outnum)
        self.pkr.update(pv, outnum)
        vsout = self._name_vtp('virial_stress_')
        wr.write_stress_ellipses(pv, vsout, pv.stresses['virial'])
        
        mout = self._name_vtp('cell_dual_')
        wr.writemeshenergy(pv, mout)

        tout = self._name_vtp('cell_')
        wr.writetriforce(pv, tout)

        if args.hardy:
            sout = self._name_vtp('hardy_stress_')
            wr.write_stress_ellipses(pv, sout, pv.stresses['hardy'])

    def _handle_constants(self, pv):
        # Handle K, and Gamma
        k, gamma, L  = args.k, args.gamma, args.L
        nf = pv.tri.n_vertices()
        # Could easily read these from a .conf file
        K = np.full(nf, k)
        Gamma = np.full(nf, gamma)
        cl = pv._construct_cl_dict(L)
        pv.set_constants(K, Gamma, cl)


    def _iterate_wl(self,pv):
        args = self.args
        if not args.w:
            return
        start, stop, step = eval(args.w)
        wls = range(start, stop, step)
        for wl in wls:
            pv._set_wl(wl)
            pv.stress_on_centres(self.omega)
            #pv.radial()
            wlnum = ('%d' % int(1000*wl)).zfill(4)
            self.wlst.update(pv, wlnum)

            wlout = self._name_vtp('hardy_stress_', wlnum)
            print 'the name of one stress output ', wlout
            wr.write_stress_ellipses(pv, wlout, pv.stress)

        # Fow now just stop here
        #sys.exit()



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
    parser.add_argument("-wl", type=float, default=1., help='smoothing length')
    parser.add_argument("-s", action='store_true')
    parser.add_argument("-f", '--fast', action='store_true')
    parser.add_argument("-k", type=float, default=1., help='Area constant')
    parser.add_argument("-g", '--gamma', type=float, default=0., help='Perimeter constant')
    parser.add_argument("-L", type=float, default=0., help='Contact length constant')
    parser.add_argument("-w", type=str, default='', 
            help='A list containing start finish and step for wl values')
    parser.add_argument('-step', type=int, default=1, help='Step through time')
    # flags
    parser.add_argument("--simple", action='store_true', 
            help='Only perform simple stress calculation')
    parser.add_argument("--hardy", action='store_true', 
            help='Turn on Hardy stress calculation')
    parser.add_argument("--exclude", action='store_true', 
            help='Mirrors exclude_boundary flag in samso')
    args = parser.parse_args()

    if not os.path.exists(args.dir):
        os.mkdir(args.dir)

    print 'vertex model parameters'
    print 'k ', args.k
    print 'gamma, ',  args.gamma
    print 'L ',args.L
    if args.w:
        print 'Calculating stress for wl = ', range(*eval(args.w))

    stressrun = Stress_Senario(args)
    stressrun.save_parameters()
    for _ in stressrun: pass

