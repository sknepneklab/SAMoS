# Module. Senario object is the framework which enables operating on folders of .dat data
import ioutils as io
import sys, os
import numpy as np
from numpy import linalg as lg
import glob
import math 

from collections import OrderedDict
import os.path as path
from matplotlib import pyplot as plt

from tmp_read_data import ReadData

from cellmesh import *
import cmplotting as cm
import constructmesh as ct

# helper methods
divider = '###########  ############'

def f_outnum(dataf):
        inp_dir, base_file = path.split(dataf)
        base_name, ext = path.splitext(base_file)
        outnum= base_name.split('_')[-1]
        return outnum 

verb = False

# Hardy stress things, forget them
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


# Classes for outputing the data.

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

# could pickle pv.stresses
# what other data is useful to pickle?
class SPickler(Tracker):
    def __init__(self, fileo='stresses'):
        self.fileo = fileo

    def _update(self, pv):
        with open(self.foutname, 'wb') as fo:
            pickle.dump(pv.stresses, fo)

# Define analysis routines here

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

        runfiles = []
        if args.last or args.first:
            if args.last:
                runfiles.append(infiles[-1])
            if args.first:
                runfiles.append(infiles[0])
            return runfiles

        if not (args.start or args.stop):
            return infiles
        else:
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
        pass # the code to execute on each data file

    def _finishup(self):
        pass # the cleanup code that needs to be executed
        # actually going to os.rmdir the tmp/ and plots/data if they are empty
        exdirs = ['tmp/', 'plots/data/', 'plots/']
        for ex in exdirs:
            if not os.path.exists(ex): continue

            if not os.listdir(ex):
                os.rmdir(ex)
                print 'successfully cleaned up {}'.format(ex)

    def save_parameters(self, paramsf='out.parameters'):
        with open(paramsf, 'w') as fp:
            if verb:
                'recording the analysis parameters to ', paramsf
            io.dumpargs(self.args.__dict__, fp)

    # helper functions 

    def _name_vtp(self, name, i=None):
        args= self.args
        outname = name + self.outnum
        if i:
             outname += '_' + str(i) 
        outname += '.vtp'
        return path.join(args.dir, outname)

    def _read_pv(self, env=False):
        rdat = ReadData(self.dataf, isenv=env)
        if self.args.faces:
            facefile = self.args.faces
        else:
            facefile = path.join(self.inp_dir, 'faces_' + self.outnum + '.fc')
        if self.args.triangulate:
            pv = PVmesh.datbuild(rdat)
        else:
            simp, _ = io.readfc(facefile)
            pv = PVmesh.datbuild(rdat, simp)
        self.rdat = rdat
        return pv

    def _vtp_output(self, pv):

        stressnames = ['virial', 'hardy']
        for sname in stressnames:
            if sname in pv.stresses:
                vsout = self._name_vtp(sname+'_stress_')
                wr.write_stress_ellipses(pv, vsout, pv.stresses[sname].stress, 
                        scale=self.args.scale)
    
        mout = self._name_vtp('cell_dual_')
        wr.writemeshenergy(pv, mout)

        tout = self._name_vtp('cell_')
        wr.writetriforce(pv, tout)

    def _handle_constants(self, pv):
        args = self.args
        # Handle K, and Gamma
        k, gamma, L  = args.k, args.gamma, args.L
        nf = pv.tri.n_vertices()

        K = np.full(nf, k)
        Gamma = np.full(nf, gamma)

        cl = pv._construct_cl_dict(args.lpairs)
        pv.set_constants(K, Gamma, cl)

# Bare minimum mesh senario.
# i.e. want to get the area distributions of the cell types 

class Property(Senario):
    def __init__(self, args):
        super(Property, self).__init__(args)
        if not os.path.exists(args.dir):
            os.mkdir(args.dir)

    def _operate(self):
        args = self.args

        pv = self._read_pv()
        # sets K, Gamma, Lambda 
        if hasattr(args, 'k'):
            self._handle_constants(pv)

        # Tell pv to organise the types into lists
        #pv._type_lists()

        if args.areas:
        # Get the typed area distribution object
            areatp = pv.type_area()
            self.areas_main(areatp)
            print pv.mesh.prims

        if args.texture:
            pv.makecellparts()
            vsout = self._name_vtp('structure_')
            wr.write_stress_ellipses(pv, vsout, pv.structure, normalise=True)
        
        if args.top:
            self.topology(pv)
        
        if args.boundary:
            self.boundary_main(pv)

        if False:
            # calculate the difference between real areas and critical division areas
            max_area = 2.8  # read from configuration file
            rsp, rad = self.radial_area(pv, max_area)
            outd, ax = cm.arearadial(rsp, rad)

            name = os.path.join('plots/data/','arearadial_'+str(self.outnum)+'.dat' )
            with open(name, 'w') as fo:
                io.dump(outd, fo)

    def radial_area(self, pv, max_area, nstat=4000.):
        nbins = int(math.ceil(len(pv.tri.rcmd)/nstat))
        if nbins < 4: nbins = 4
        rcmd = pv.tri.rcmd

        # subtract real and div_max_area
        bulk = pv.tri.bulk
        realareas = [pv.mesh.areas[pv.tri.to_mesh_face[i]] for i in bulk]

        adiff = np.subtract(realareas, np.full(len(bulk), max_area))
        adiff[adiff < 0] = 0 #replace negative values with 0
        print adiff
        adiff = dict(zip(bulk, adiff))

        # now radially average
        mrc = np.max(rcmd)
        rspace = np.linspace(0, mrc, nbins+1, endpoint=True)
        npd = np.digitize(rcmd, rspace, right=True)-1
        rsum, hcount = np.zeros(nbins), np.zeros(nbins)
        for i, bn in zip(bulk, npd):
            rsum[bn] += adiff[i]
            hcount[bn] += 1
        radiff = np.true_divide(rsum, hcount)
        self.rspace = rspace
        self.radiff = radiff
        return rspace, radiff


    def areas_main(self, areatp):
        for tp, areas in areatp.items():

            print 'mean areas, type = {}'.format(tp)
            print np.mean(areas), len(areas)
            plt.clf()
            plt.hist(areas, bins=30, alpha=0.6)
            plt.title('area distribution. type={}'.format(tp))
            plt.show()

    def topology(self, pv):
        # note that the number of sides of boundary cells depends on the number of boundary
        #  particles. Could exclude the outer ring of cells from being analysed here.
        tri = pv.tri;
        nnlist = []
        for i in tri.bulk:
            # loop over all the particles with boundary flag = 0
            # count their neighbours
            vh = tri.vertex_handle(i)
            nnlist.append(len(list(tri.vv(vh))))
        cm.nntopology(nnlist, self.outnum)

    def boundary_main(self, pv):
        # Extract the boundary connectivity information for this timestep
        # Useful for restarting a simulation from a .dat output
        tri = pv.tri
        bconnect = tri.get_boundary_connectivity()
        # write out that connectivity file
        out = path.splitext(self.dataf)[0] + '.boundary'
        with open(out, 'w') as fo:
            fo.write('#\n')
            for idx, (i, j) in enumerate(bconnect):
                i += pv.startid
                j += pv.startid
                fo.write('{} {} {}\n'.format(idx, i, j))

        # why did we ever do this?
        # Also set the boundary flags in the input file
        #Ci = ct.CellInput(self.dataf)
        #boundaries= tri.get_boundaries()
        #Ci.boundary_flags(boundaries)
        #Ci.dump()
        



# old at this point, be careful
# Senario for calculating stress and all related quantities while we are at it
# Should use this senario when implementing strain and Graner metrics
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

        self.ref_index = None
        if args.ref:
            nums = map(f_outnum, self.infiles)
            # string comparision 
            self.ref_index = nums.index(args.ref)
        else:
            print 'Warning: reference timestep not set for calculating Strain'
            self.ref_index = 0
        self.ref_structure = None


        self.pkr = SPickler(os.path.join(args.dir, 'stresses'))

    def _operate(self):
        args = self.args
        outnum = self.outnum
        omega = self.omega

        pv = self._read_pv()
        self._handle_constants(pv)
        #pv.calculate_energy()
        pv.calculate_forces()

        pv.makecellparts()
        pv.on_centres(args.adj)

        if self.fid == self.ref_index:
            self.ref_structure = {}


        if args.test:
            self.test_convergence(pv)

        # for real 
        print 'stress averaging adj = ', args.adj
        adjlist = [0, args.adj]
        self.calculate_U(pv, adjlist)

        # pickle useful datas
        self.pkr.update(pv, outnum)

        ### vtp output 
        self._vtp_output(pv)

    def calculate_U(self, pv, adjlist):
        self.adjavg_structure(pv, adjlist)
        #pv.structure_xx


    def test_convergence(self, pv):
        adjlist = [0, 1, 5, 10]
        self.adjavg_structure(pv, adjlist)
        cm.dev_texture(pv.xx_shear)

    def adjavg_structure(self, pv, adjlist=[0,1,5]):
        structure = pv.structure
        structure_xx = {}
        xx_trace = OrderedDict()
        xx_shear = OrderedDict()

        
        # generic method for decomposing tensors. should be moved somewhere.
        def decompose(ovtens):
            #take (.., 3, 3) array of tensors and find principle directions and eigenvalues
            ovtens = ovtens[:,:2,:2] # assume we need to go to two dimensions
            evals, evecs = lg.eig(ovtens)
            def orderedst(eva_evec): 
                eva, evec = eva_evec
                # for one set of eigenvalues/vectors order them and find components
                if abs(eva[1]) > abs(eva[0]):
                    eva[0], eva[1] = eva[1], eva[0]
                    evec[0], evec[1] = evec[1], evec[0]
                shear = (eva[0] - eva[1])/2
                trace = (eva[0] + eva[1])/2
                return eva, evec, trace, shear
            ll = map(orderedst, zip(evals, evecs))
            # transpose list of list
            rll = map(list, zip(*ll)) 
            # cleanup
            return map(np.array, rll)
            
        for adjn in adjlist:
            adjstr = pv.structure_on_centres(adjn)
            structure_xx[adjn] = adjstr
            #evecs, evals = lg.eig(np.array(adjstr.values())[:,:2,:2])
            evals, evecs, trace, shear = decompose(np.array(adjstr.values()))

            xx_trace[adjn] = trace 
            xx_shear[adjn] = shear

        pv.xx_trace = xx_trace
        pv.xx_shear = xx_shear

    #def avg_dev_texture(self, xx_shear):
        #x = xx_shear.keys()
        #ymean = map(np.mean, xx_shear.values())
        #yupper = map(max, xx_shear.values())
        #ylower = map(min, xx_shear.values())
        #plt.plot(x, ymean, marker='o')
        #plt.plot(x, yupper, linestyle='--', color='g', marker='o')
        #plt.plot(x, ylower, linestyle='--', color='g', marker='o')
        #plt.show()

    def old_testing_code(self):
        pass
        #bulkforces = [forces[i] for i in pv.tri.bulk]
        #avg_force = np.mean(map(norm, bulkforces))
        #print 'avg_force', avg_force
        #pv.stresses['avg_force'] = avg_force

        # Polygons mesh. Just for debugging
        #wr.writepoints(pv.mesh.centroids.values(), self._name_vtp('centroids_'))
        #wr.writemesh(pv.polygons, self._name_vtp('polygons_'))



# Senario for reading the .dat files and returning some simple data
#  For example the number of ghost particles
# Use if building the meshes is not necessary
class Basic(Senario):
    def __init__(self, args):
        super(Basic, self).__init__(args)
        self.timeline = []
        self.numghosts = []
        self.n_cells = []
        self.avgarea = []
        self.avgshape = []

    def _operate(self):
        args = self.args
        self.timeline.append(self.outnum)
        rdat = ReadData(self.dataf)
        keys = rdat.keys
        assert 'boundary' in keys

        boundary = np.array(rdat.data[keys['boundary']])
        a0 = np.array(rdat.data[keys['area']])
        realareas = np.array(rdat.data[keys['cell_area']])
        realperim = np.array(rdat.data[keys['cell_perim']])
        realshape = np.true_divide(realperim,np.sqrt(realareas))
        self.avgshape.append( np.mean(realshape[np.where(boundary==0)[0]]) )

        nghosts = len(np.where(boundary==1)[0])
        self.numghosts.append(nghosts)
        self.n_cells.append(len(np.where(boundary==0)[0]))
        self.avgarea.append( np.mean(realareas[np.where(boundary==0)[0]]) )



# Senario for just reading the .dat files and calculating MSD curves
class MSD(Senario):
    def __init__(self, args):
        super(MSD, self).__init__(args)
        # going to deal with the whole dataset as a block 
        self.timeline = []
        # A list of dictionaries containing data blocks of position data for each type
        # [{<type>:<posdata>}]
        self.posarrlist = []
        self.types = None

        #I want the number of cells at each step
        self.n_cells = []
        
    def _operate(self):
        args = self.args
        outnum = self.outnum

        rdat = ReadData(self.dataf)
        keys = rdat.keys
        x = np.array(rdat.data[keys['x']])
        y = np.array(rdat.data[keys['y']])
        z = np.array(rdat.data[keys['z']])
        assert 'type' in keys, 'Need type data for msd calculation'
        alltp = np.array(rdat.data[keys['type']], dtype=int)
        allboundary= np.array(rdat.data[keys['boundary']], dtype=int)
        # break the data based on the type of particle
        # throw away the boundary type, it doesn't make sense to calculate MSD for them
        # We assume that particles are not dividing for now
         
        tpids = {}
        for i, tp in enumerate(alltp):
            if allboundary[i] == 1:
                continue # Ignore boundary particles, no type 1 particles remain
            if args.type and tp != args.type:
                # specifying a particular type
                continue
            try:
                tpids[tp].append(i)
            except KeyError:
                tpids[tp] = [i]

        # Now divide up the block using numpy indexing
        posarr = np.column_stack([x,y,z])
        postp = {}
        for tp, posids in tpids.items():
            postp[tp] = posarr[posids]

        self.timeline.append(outnum)
        self.posarrlist.append(postp)
        if not self.types:
            self.types  = postp.keys()
        else:
            assert self.types == postp.keys(), \
                'All the same types should exist throughout the simulation'

        # cell growth, assume no environment
        bd = np.array(rdat.data[keys['boundary']])
        self.n_cells.append( len(np.where(bd==0)[0]) )

    def retrieve(self):
        return self.timeline, self.posarrlist

    # Operate on the data block generated from MSD senario
    # Call this only after finishing the iteration
    def allmsd(self):
        def fmsd(disp):
            return np.mean(np.sum(disp*disp, 1))
        # construct the 3d block: time, particle, xyz coord
        # assume particles are not dividing
        blocks = dict(zip(self.types, [[] for _ in range(len(self.types))]))
        print blocks
        # has to be done for each type
        for postp in self.posarrlist:
            for tp, posarr in postp.items():
                blocks[tp].append(posarr)

        # {tp:msd}
        self.msd = dict(zip(self.types, [[] for _ in range(len(self.types))]))            
        # 
        for tp, block in blocks.items():
            # Finally construct the numpy block
            allblock = np.rollaxis(np.dstack(block), 2, 0)
            taulist = range(len(self.timeline)/2)
            for tau in taulist:
                msdval = []
                for i in range(tau, len(self.timeline)-tau):
                   disp =  allblock[i-tau,:,:] - allblock[i+tau,:,:]
                   msdval.append( fmsd(disp) )
                self.msd[tp].append(np.mean(msdval))

        #return self.timeline[1:], self.msd
        return taulist, self.msd
    
    

### Some work towards tracing t1 transtions -- never very effective -- better to start over
# The Tone and Tlist classes for tracking transitions
from transitions import *

# associated with the transtions.py module 
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
        
        self.pkr = SPickler(os.path.join(args.dir, 'stresses'))

    def _operate(self):
        args = self.args
        outnum = self.outnum

        # basic setup
        rdat = ReadData(self.dataf)
        facefile = path.join(self.inp_dir, 'faces_' + self.outnum + '.fc')
        simp, _ = io.readfc(facefile)
        pv = PVmesh.datbuild(rdat, simp)
        self._handle_constants(pv)

        # Calculate cell-level virial stress
        pv.calculate_forces(exclude_boundary=args.exclude)
        pv.makecellparts()
        pv.on_centres(adj=args.adj)

        # Save the pv object in an OrderedDict
        self.helist[outnum] = pv

        # track the transitions

        #if self.fid == 0:
            #self.tlist = Tlist.populate(pv, outnum)
        #else:
            #self.tlist.update(pv, outnum)

        # Analysis output
        #self.pkr.update(pv, outnum)

        # Basic vtp output 
        #self._vtp_output(pv)

        # Output structure  tensor
        vsout = self._name_vtp('structure_')
        wr.write_stress_ellipses(pv, vsout, pv.structure, normalise=True)

    def _finishup(self):
        #print 'finished'
        if self.tlist and hasattr(self.tlist, 'nt'):
            print 'total t1 transitions', self.tlist.nt
        else:
            print 'Did not produce a tlist object'

    def _halfedges(self, simp):
        lamhedges = lambda si: zip(si, np.roll(si, -1))
        hearr = np.array(map(lamhedges, simp))
        # all halfedges
        hearr= hearr.reshape(-1, hearr.shape[-1])
        return set(map(tuple, hearr))


