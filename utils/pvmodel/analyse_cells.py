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

import glob

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


def f_outnum(dataf):
        inp_dir, base_file = path.split(dataf)
        base_name, ext = path.splitext(base_file)
        outnum= base_name.split('_')[-1]
        return outnum 

# The template senario
class Senario(object):
    def __init__(self, args):
        # Use the args object from argsparse and add to it as necessary to keep
        #  track of all the parameters
        self.args= args
        # Declare tracking variables
        self.outnum=None
        self.inp_dir=None
        self.dataf=None
        self.fid = 0
        # Glob and slice dataset
        self.infiles =self._process_input(args)

        self.numf = len(self.infiles)
        self.alog = 'analysis.log'

    def _process_input(self, args):
        inp_file = args.input
        infiles = sorted(glob.glob(inp_file))
        if infiles == []:
            print 'Did not find any files using %s', inp_file
            sys.exit()

        if not (args.start or args.stop):
            return infiles
        # Use the start and stop numbers to cut down the input files ot a managable number
        onums = map(lambda s: int(f_outnum(s)), infiles)
        starti = onums.index(args.start) if args.start else 0
        stopi = onums.index(args.stop)+1 if args.stop else len(onums)
        infiles = infiles[starti:stopi]
        return infiles


    def __iter__(self):
        return self

    def next(self):
        step = self.args.step
        if self.fid < self.numf:
            self.dataf = self.infiles[self.fid]
            self.inp_dir, base_file = path.split(self.dataf)
            self.outnum = f_outnum(self.dataf)
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
        outname = name + self.outnum
        if i:
             outname += '_' + str(i) 
        outname += '.vtp'
        return path.join(args.dir, outname)

    def _read_pv(self):
        rdat = ReadData(self.dataf)
        facefile = path.join(self.inp_dir, 'faces_' + self.outnum + '.fc')
        simp, _ = io.readfc(facefile)
        pv = PVmesh.datbuild(rdat, simp)
        #pv = PVmesh.datbuild(rdat)
        return pv

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

# Object representing a T1 transition
class Tone(object):

    def __init__(self, vids, tstamp):
        # the current state of transition
        self.vids  = vids 
        # the flipped state
        self.xvids = None
        # The number of times flipped
        self.nx = 0
        # the step number when we started tracking
        self.tstamp = tstamp

    def flipped(self):
        return self.nx % 2 == 1

    def flip(self, vids):
        self.xvids = self.vids
        self.vids = vids
        self.nx +=  1

# invert a dictionary
def invert(dct):
    return dict([(v,k) for k,v in dct.items()])

# a collect functions which operate 
#import collections.abc.MutableSequence as MutableSequence
import collections
import bisect
class Tlist(collections.MutableMapping):

    tcut = 0.1
    ttrack = 0.8

    def __init__(self, lst):
        self.lst = lst
        self.nt = 0

    def __getitem__(self, i):
        return self.lst[i]

    def __setitem__(self, i, v):
        self.lst[i] = v

    #def insert(self, i, v):
        #self.lst[i] = v

    def __iter__(self):
        return self

    def __len__(self):
        return len(self.lst)

    def __delitem__(self, i):
        del self.lst[i]

    def append(self, tone):
        assert self.keyinc not in self.lst
        self.lst[self.keyinc] = tone
        self.byvids[tone.vids] = self.keyinc
        self.keyinc += 1

    def remove(self, vpair):
        i = self.byvids[vpair]
        xvids = self.lst[i].xvids
        if xvids and xvids in self.byvids:
            del self.byvids[xvids]
        del self.lst[i]
        del self.byvids[vpair]

    # helper methods for the populate and update functions
    @staticmethod
    def _find_short(pv):
        lls = sorted(pv.mesh.lle.items(), key= lambda t: t[1])
        eids, llen = zip(*lls)
        short = eids[:bisect.bisect_left(llen, Tlist.ttrack)]
        #short = eids[:bisect.bisect_left(llen, 3.1)]

        rshort = []
        dual_to_tri = invert(pv.tri.to_mesh_edge)
        #assert set(dual_to_tri.keys()) == set(pv.tri.to_mesh_edge.values())
        # Need to remove boundary trimesh edges
        for se in short:
            if se in dual_to_tri:
                rshort.append(se)
            else:
                #print 'boundary', se
                pass #boundary edge

        pv.mesh.dual_to_tri = dual_to_tri
        return rshort
    @staticmethod
    def _maketovpair(pv):
        #dual_to_tri = invert(pv.tri.to_mesh_edge)
        dual_to_tri = pv.mesh.dual_to_tri
        return dict([(eid, pv.tri.getvpair(teid)) for eid, teid in dual_to_tri.items()])

    # constructor which takes pv object
    @classmethod
    def populate(cls, pv, outnum):

        short = cls._find_short(pv)
        # dual to tri edges
        tovpair = cls._maketovpair(pv)
        lst = dict(enumerate([Tone(tovpair[eid], outnum) for eid in short]))
        byvids = dict([(tone.vids, i) for i, tone in lst.items()])

        print 'started tracking {} edges'.format(len(lst))
        tt = cls(lst)

        tt.last_tovpair = tovpair
        tt.byvids = byvids
        tt.keyinc = len(lst)
        return tt

    def get_aux(self, tri, vpair):
        vid = tuple(vpair)
        heh = tri.halfedge_handle(tri.vedge[vid])
        oheh = tri.opposite_halfedge_handle(heh)
        ft =tri.face_handle(heh)
        fto = tri.face_handle(oheh)
        adjv = [v.idx() for v in tri.fv(ft)]
        adjvo = [v.idx() for v in tri.fv(fto)]
        tquad = set(adjv).union(set(adjvo))
        auxvids = frozenset(set(adjv).union(set(adjvo)) - vpair)
        return auxvids

    # now we need methods to update this list whenever edges change length or flip
    def update(self, pv, outnum):
        # setup
        pv.tri.vedgemap()
        tri = pv.tri; mesh = pv.mesh

        # assume that any edge that is created on this step is below track
        # Just check all edges every time to start with
        short = Tlist._find_short(pv)

        # Important! remeber that the mesh ids may not be consistent
        #  We are sure that the cell centre ids are consistent, that is all.

        tovpair = Tlist._maketovpair(pv)

        vidsls = [tovpair[eid] for eid in short]
        # look for new pairs of vertices
        newvids = set(tovpair.values()) - set(self.last_tovpair.values())
        # find which vid pairs that the new pairs correspond to (flips have occured)
        print 'new', len(newvids)
        unaccounted = []
        for vpair in newvids:
            # Either auxvids is being tracked or 
            auxvids = self.get_aux(tri, vpair)
            if vpair in self.byvids:
                continue
            if auxvids not in self.byvids:
                unaccounted.append(vpair)
                continue
            tt = self.lst[self.byvids[auxvids]]
            # Keeping both references to the transition object
            self.byvids[vpair] = self.byvids[auxvids]
            tt.flip(vpair)
        
        ### deal with the unaccounted edges
        ppairs = []
        pairwith = newvids - set(unaccounted)
        #print pairwith
        #print unaccounted
        for vpair in unaccounted:
            # pair them up by finding the common vertices
            va, vb = vpair
            pairw = None
            for pair in pairwith:
                if va in pair or vb in pair:
                    pairw = pair
                    break
            if not pairw:
                print 'warning, was looking for a pair of edges, only found 1'
                print 'Need to track more egdes'
                print vpair
            ppairs.append((vpair, pairw))

        # Now we have the edge pairs, going to attempt to flip both of them and
        # see which pair of flips recovers a pair of edges which are being tracked

        # we want to find the flipped edge of the second pair given that the first
        # edge is flipped
        
        # going to assume that the connectivity of the triangulation is such that
        #  there is a loop of 5 vertices with two central edges.
        for vpair, pairw in ppairs:
            vd, vb = vpair # this the edge we want to flip
            ive, ivb = pairw # the edge we want to ignore (pre-flip)
            # we can obtain all the vertices in the loop by tracing the three
            #  faces associated with the two internal edges
            hehv_id= tri.vedge[(vd,vb)]
            hehi_id= tri.vedge[(ive, ivb)]
            hehv =tri.halfedge_handle(hehv_id)
            hehi =tri.halfedge_handle(hehi_id)
            # find the third face
            ovid = tri.opposite_halfedge_handle(hehv)
            oiid = tri.opposite_halfedge_handle(hehi)
            hhfaces = [hehv, hehi, ovid, oiid] # should be one duplicate face
            loop = set()
            for heh in hhfaces:
                fh = tri.face_handle(heh)
                verts = [v.idx() for v in tri.fv(fh)]
                loop = loop.union( set(verts) )
            assert len(loop) == 5
            # now identify the remaining edge by process of elimination
            original = loop - set(pairw)
            original = original- set(vpair)
            original = frozenset(original)
            # finally found the missing transition so count it and remove it
            assert original in self.byvids
            self.remove(original)
            self.nt += 1


        # resolve edges which aren't short anymore
        for tone in self.lst.values():
            #tmp 
            try:
                hehid = tri.vedge[tuple(tone.vids)]
            except KeyError:
                print 'dropped a tone object'
                self.remove(tone.vids)
                continue
            el = tri.ll[hehid]
            if el > Tlist.tcut:
                if tone.xvids:
                    tri.ll[hehid]
                    # count this transition
                    self.nt += 1
                self.remove(tone.vids)

        #Find which of the entries in short are newvids and which arent
        shortpairs = [tovpair[eid] for eid in short]
        newshort = set(shortpairs) - set(newvids)
        #find the ones we aren't tracking yet
        untracked= newshort - set(self.byvids.keys())
        print 'transitions', self.nt
        print 'untracked', len(untracked)
        print 'unaccounted', len(unaccounted)
        # track any short edges
        for i, vpair in enumerate(untracked):
            self.append(Tone(vpair, outnum))

        # finishup
        self.last_tovpair = tovpair



class Transitions(Senario):
    def __init__(self, args):
        super(Transitions, self).__init__(args)
        # we need to start keeping track of all the pv objects we read in
        self.helist= OrderedDict()
        self.prevon = None
        self.totalt1 = 0

        # set the length under which we start to track edges
        self.ltracking = 0.1
        # main object for dealing with transistions
        self.tlist = None

    def _operate(self):
        args = self.args
        outnum = self.outnum

        # basic setup
        rdat = ReadData(self.dataf)
        facefile = path.join(self.inp_dir, 'faces_' + self.outnum + '.fc')
        simp, _ = io.readfc(facefile)
        pv = PVmesh.datbuild(rdat, simp)
        if self.fid == 0:
            self.tlist = Tlist.populate(pv, outnum)
        else:
            self.tlist.update(pv, outnum)

            ## two halfedges for every edge
            #nt1 = len(gained)/2
            #print 'no. of t1 transitions'
            #print nt1
            #self.totalt1 += nt1

    def _finishup(self):
        print 'finished'
        #print 'total t1 transitions', self.tlist.nt

    def _halfedges(self, simp):
        lamhedges = lambda si: zip(si, np.roll(si, -1))
        hearr = np.array(map(lamhedges, simp))
        # all halfedges
        hearr= hearr.reshape(-1, hearr.shape[-1])
        return set(map(tuple, hearr))

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

    def _operate(self):
        args = self.args
        outnum = self.outnum
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

def process_commands():
    parser = argparse.ArgumentParser()
 
    parser.add_argument("-i", "--input", type=str, default='cell*.dat', help="Input dat file")
    parser.add_argument("-c", "--conf", type=str, default='*.conf', help="Input dat file")
    parser.add_argument("-d", "--dir", type=str, default='tmp/', help="Output directory")
    parser.add_argument("-wl", type=float, default=1., help='smoothing length')
    parser.add_argument("-s", action='store_true')
    parser.add_argument("-k", type=float, default=1., help='Area constant')
    parser.add_argument("-g", '--gamma', type=float, default=0., help='Perimeter constant')
    parser.add_argument("-L",  type=float, default=0., help='Contact length constant')
    parser.add_argument("-a0",  type=float, default=3., help='Give the target area explicitely')
    parser.add_argument("-w", type=str, default='', 
            help='A list containing start finish and step for wl values')
    parser.add_argument('-step', type=int, default=1, help='Step through time')
    parser.add_argument('-start', type=int, default=None, help='Start number')
    parser.add_argument('-stop', type=int, default=None, help='Stop number')
    parser.add_argument('-kf', "--kernel", type=str, default='quartic',
            help = 'the type of smoothing function to use')
    parser.add_argument('--scale', type=float, default=1., 
        help='scale stress ellipses')

    parser.add_argument('-lt', '--ltransition', type=float, default=0.02, 
            help = 'Threshold length for a t1 transition')

    valid_smoothing_functions = ['uniform', 'quartic']
    # flags
    parser.add_argument("--centroid", action='store_true', 
            help='Use the polygonal centroid as the center for stress calcuation')
    parser.add_argument("--simple", action='store_true', 
            help='Only perform simple stress calculation')
    parser.add_argument("--hardy", action='store_true', 
            help='Turn on Hardy stress calculation')
    parser.add_argument("--include", action='store_true', default=False,
            help='Mirrors exclude_boundary flag in samso')
    # flags for wildly changing the type of analysis
    parser.add_argument("--trans", action='store_true',
            help='Identify transitions.')
    parser.add_argument("--all", type=str, default=None,
            help='Run some analysis on every sub-directory')

    args = parser.parse_args()
    
    # cleanup
    args.exclude = not args.include 


    return args 

def configuration(args):
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
    
    #target perimeter
    tperim = -args.L/args.gamma
    print 'target perimeter', tperim
    sindex= tperim / np.sqrt(args.a0)
    print 'shape index', sindex
    args.sindex = sindex

    return args


def print_parameters():
    print 'vertex model parameters'
    print 'k ', args.k
    print 'gamma ',  args.gamma
    print 'L ',args.L
    print 'Excluding the boundary:', args.exclude
    if args.w:
        print 'Calculating stress for wl = ', range(*eval(args.w))
    print

if __name__=='__main__':

    epidat = '/home/dan/cells/run/soft_rpatch/epithelial_equilini.dat'

    args = process_commands()

    #stressrun.save_parameters()
    if not args.all:
        
        if not os.path.exists(args.dir):
            os.mkdir(args.dir)
        configuration(args)
        print_parameters()
        if args.trans:
            stressrun = Transitions(args)
        else:
            stressrun = Stress_Senario(args)
        for _ in stressrun: pass
    else:
        # Want to be able to iterate over all folders in a directory, run analyse and collect data
        print 'ntt'
        ntt = OrderedDict()
        import opdir as op
        def macro():
            configuration(args)
            print_parameters()
            stressrun = Transitions(args)
            for _ in stressrun: pass
            ntrans = stressrun.tlist.nt
            ntt[args.sindex] = ntrans
        op.diriterate(macro, args.all)
        print ntt

        import cmplotting as cm
        cm.nttplot(ntt, log=True)




    



