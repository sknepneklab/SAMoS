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

divider = '###########  ############'

# Define analysis routines here

# We need one for catching the pressure of some random cells
def random_choice(pv, nr=5):
    # nr is the number of random points to choose
    bulk = pv.get_mesh_bulk(pv.tri)
    bulk_choice = np.random.choice(bulk, size=nr, replace=False)
    return bulk_choice

verb = False

from scipy import integrate
def uniformkernel(wl):
    def unif(r):
        return 1/(np.pi * wl**2)
    return unif

def quartic_wl(wl):
    if verb: print 'generating smoothing function with wl = ', wl
    def quartic(r):
        return 5/(np.pi*wl**2) * (1 + 3*(r/wl))*(1-(r/wl))**3 if r < wl else 0. 
        #return 5/(4* np.pi*wl) * (1 + 3*(r/wl))*(1-(r/wl))**3 if r < wl else 0. 
    space = np.linspace(0, wl, 100, True)
    qq = np.array(map(quartic, space))
    ii = integrate.simps(qq* space*2*np.pi, space)
    # assert ii is close to 1
    #io.plotrange(quartic, 0, wl)
    return quartic

def hash_function(f, linspace):
    l = linspace[-1] - linspace[0]
    ne = len(linspace)
    hf = map(f, linspace)
    def nf(r):
        i = int(round((ne-1) * r/l))
        return hf[i]
    return nf


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
        self.alog = 'analysis.log'

    def __iter__(self):
        return self

    def next(self):
        step = self.args.step
        if self.fid < self.numf:
            self.dataf = self.infiles[self.fid]
            self.outnum = self._outnum() # must call this function
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
            if verb:
                'recording the analysis parameters to ', paramsf
            io.dumpargs(self.args.__dict__, fp)

    def _name_vtp(self, name, i=None):
        args= self.args
        outname = name + args.outnum
        if i:
             outname += '_' + str(i) 
        outname += '.vtp'
        return path.join(args.dir, outname)


# template for Tracker objects
class Tracker(object):
    def __init__(self,fileo):
        # setup
        self.fileo =fileo
        self.outd= OrderedDict()
    
    def update(self, pv, xxv):
        fileo = self.fileo
        foutname = '_'.join([fileo, xxv]) + '.pkl'
        if verb: print 'saving radial data to ', foutname
        self.foutname = foutname
        self._update(pv)

    # overwrite this to change behaviour
    def _update(self, pv):
        pass
    
    def write_whole(self):
        io.dump(self.outd, self.outf)

    def cleanup(self):
        self.outf.close()

import pickle
# Can wrap the whole code in a try: finally: block.
# This will ensure that relevent quantities are saved
class FPickler(Tracker):
    def _update(self, pv):
        forces = pv.tri.forces
        bulkforces = [forces[i] for i in pv.tri.bulk]
        avg_force = np.mean(map(norm, bulkforces))
        self.outd['avg_force'] = avg_force
        vforces = [pv.mesh.vertex_forces[i] for i in mesh.bulk]
        self.outd['avg_vertex_force'] = np.mean(map(norm, vforces))
        with open(self.foutname, 'wb') as fo:
            pickle.dump(self.outd, fo)

# Seems like pickling the data might be better
# could pickle pv.stresses
# what other data is useful to pickle?
class SPickler(Tracker):
    def __init__(self, fileo='stresses'):
        self.fileo = fileo

    def _update(self, pv):
        with open(self.foutname, 'wb') as fo:
            pickle.dump(pv.stresses, fo)


set_choice = [844, 707, 960, 606, 958, 801, 864, 859, 870, 302, 500]
class Stress_Senario(Senario):
    def __init__(self, args):
        super(Stress_Senario, self).__init__(args)
        print 'initializing the run generator'
        wl= args.wl

        self.fileo=  os.path.join(args.dir, 'stress_st') 

        self.omega= None
        if args.kernel == 'quartic':
            om = quartic_wl(wl)
        elif args.kernel == 'uniform':
            om = uniformkernel(wl)
        self.omega = om

        wlf = os.path.join(args.dir, 'stress_wl')
        self.glo = FPickler(os.path.join(args.dir, 'avg'))
        self.pkr = SPickler(os.path.join(args.dir, 'stresses'))

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
        print divider
        print 'areametric', pv.areametric()
        print 'primmetric', pv.perimmetric()

        if verb: print 'Just use one value of smoothing length'
        pv._stress_setup(centroid=args.centroid)
        #self._iterate_wl(pv)
        pv._set_wl(args.wl)
        clist = set_choice if args.s else None

        pv.makecellparts()
        pv.stress_on_centres(omega, clist=clist, hardy=args.hardy)
        print 'virial avg_pressure', pv.stresses['virial'].avg_pressure
        #print pv.stresses['virial'].pressure.values()[:100]
        #pv.stress_on_vertices(omega)

        forces = pv.tri.forces
        bulkforces = [forces[i] for i in pv.tri.bulk]
        avg_force = np.mean(map(norm, bulkforces))
        print 'avg_force', avg_force
        pv.stresses['avg_force'] = avg_force

        #print pv.stresses['virial'].pressure
        #print 'hardy_vertex pressure', pv.stresses['hardy_vertices'].avg_pressure
        # Polygons mesh. Just for debugging
        #wr.writepoints(pv.mesh.centroids.values(), self._name_vtp('centroids_'))
        #wr.writemesh(pv.polygons, self._name_vtp('polygons_'))

        if args.hardy:
            pv.stresses['virial'].compare(pv.stresses['hardy'])

        # pickle useful datas
        self.pkr.update(pv, outnum)

        ### vtp output 
        stressnames = ['virial', 'hardy']
        for sname in stressnames:
            if sname in pv.stresses:
                vsout = self._name_vtp(sname+'_stress_')
                wr.write_stress_ellipses(pv, vsout, pv.stresses[sname], 
                        scale=args.scale)
        #vsout = self._name_vtp('hardy_vertex_stress_')
        #wr.write_stress_ellipses(pv, vsout, pv.stresses['hardy_vertices'],usecentres=False)
        
        mout = self._name_vtp('cell_dual_')
        wr.writemeshenergy(pv, mout)

        tout = self._name_vtp('cell_')
        wr.writetriforce(pv, tout)

    def _handle_constants(self, pv):
        # Handle K, and Gamma
        k, gamma, L  = args.k, args.gamma, args.L
        nf = pv.tri.n_vertices()
        # Could easily read these from a .conf file
        K = np.full(nf, k)
        Gamma = np.full(nf, gamma)
        cl = pv.mesh._construct_cl_dict(L)
        pv.set_constants(K, Gamma, cl)


    # out of order
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

if __name__=='__main__':

    epidat = '/home/dan/cells/run/soft_rpatch/epithelial_equilini.dat'

    parser = argparse.ArgumentParser()
 
    parser.add_argument("-i", "--input", type=str, default='cell*.dat', help="Input dat file")
    parser.add_argument("-c", "--conf", type=str, default='*.conf', help="Input dat file")
    parser.add_argument("-d", "--dir", type=str, default='tmp/', help="Output directory")
    parser.add_argument("-wl", type=float, default=1., help='smoothing length')
    parser.add_argument("-s", action='store_true')
    parser.add_argument("-k", type=float, default=1., help='Area constant')
    parser.add_argument("-g", '--gamma', type=float, default=0., help='Perimeter constant')
    parser.add_argument("-L",  type=float, default=0., help='Contact length constant')
    parser.add_argument("-w", type=str, default='', 
            help='A list containing start finish and step for wl values')
    parser.add_argument('-step', type=int, default=1, help='Step through time')
    parser.add_argument('-kf', "--kernel", type=str, default='quartic',
            help = 'the type of smoothing function to use')
    parser.add_argument('--scale', type=float, default=1., 
        help='scale stress ellipses')

    valid_smoothing_functions = ['uniform', 'quartic']
    # flags
    parser.add_argument("--centroid", action='store_true', 
            help='Use the polygonal centroid as the center for stress calcuation')
    parser.add_argument("--simple", action='store_true', 
            help='Only perform simple stress calculation')
    parser.add_argument("--hardy", action='store_true', 
            help='Turn on Hardy stress calculation')
    parser.add_argument("--exclude", action='store_true', 
            help='Mirrors exclude_boundary flag in samso')
    #parser.add_argument("--norm", action='store_true', 
            #help='Noramlize stress ellipses')


    args = parser.parse_args()

    if not os.path.exists(args.dir):
        os.mkdir(args.dir)

    # find that configuration file
    if '*' in args.conf:
        # probably gave a regular expression to find the configuration file
        confs = glob.glob(args.conf)
        print 'found configuration file(s)', confs
        args.conf = confs[0]

    from read_conf import ReadConf
    rc = ReadConf(args.conf)
    myparams = ['L', 'k', 'gamma']
    samparams = ['lambda', 'K', 'gamma']
    pmap = dict(zip(samparams, myparams))
    try:
        attrs = rc.key_words['pair_potential'][0].attributes
        for attr in attrs:
            if attr.name in samparams:
                setattr(args, pmap[attr.name], float(attr.val))
            if attr.name == 'exclude_boundary':
                args.exclude = True
    except KeyError:
        raise

    print 'vertex model parameters'
    print 'k ', args.k
    print 'gamma ',  args.gamma
    print 'L ',args.L
    print 'Excluding the boundary:', args.exclude
    print 'Using the cell centroids for stress calculation:', args.centroid
    if args.w:
        print 'Calculating stress for wl = ', range(*eval(args.w))

    stressrun = Stress_Senario(args)
    stressrun.save_parameters()
    for _ in stressrun: pass


