#!/usr/bin/env python

# Create an object to represent the input file, then call functions to manipulate the input 
# and finally have a consistent way of writing out the file

# Construct test cases for debugging and later construct simulation initial conditions

#import writemesh as wr

import os
import numpy as np
from numpy.linalg import norm
import numpy.linalg as lg
import sys
import random as rn
from random import random, uniform
import ioutils as io

from os import path
from collections import OrderedDict


from read_data import ReadData
# This function has been hijacked to convert ReadData objects into an ordered dictionary because that's what I like to work with
def rwoperate(ifile):
    # is this ReadData really giving list instead of numpy array?
    rdat = ReadData(ifile)
    keys = rdat.keys
    lrdat = len(rdat.data[0])
    outd = OrderedDict()
    #for k in keys:
        #outd[k] = []
    #for i, _ in enumerate(rdat.data[0]):
        #for kn, kv in keys.items():
            #outd[kn].append( rdat.data[kv][i] )
    for ki, k in enumerate(keys):
        outd[k] = rdat.data[ki]
    if 'id' in outd:
        outd['id'] = map(int, outd['id'])
    if 'type' in outd:
        outd['type'] = map(int, outd['type'])
    if 'boundary' in outd:
        outd['boundary']= map(int, outd['boundary'])
    if 'in_tissue' in outd:
        outd['in_tissue'] = map(int, outd['in_tissue'])
    
    # the rest are hopefully floats at the moment
    return outd

# Samos utils
from particle import Particle

# lattice code - too lazy to write it myself as usual
#http://stackoverflow.com/questions/6141955/efficiently-generate-a-lattice-of-points-in-python
# this creates a square lattice and shears into into the desired shape
# image shape is the dimensions of the lattice
# the lattice is built up by multiplying and adding the lattice vectors
def generate_lattice(image_shape, lattice_vectors) :
    center_pix = np.array([0.,0.])
    # Get the lower limit on the cell size.
    dx_cell = max(abs(lattice_vectors[0][0]), abs(lattice_vectors[1][0]))
    dy_cell = max(abs(lattice_vectors[0][1]), abs(lattice_vectors[1][1]))
    # Get an over estimate of how many cells across and up.
    nx = int(1.5 * image_shape[0]//dx_cell)
    ny = image_shape[1]//dy_cell
    # Generate a square lattice, with too many points.
    # Here I generate a factor of 4 more points than I need, which ensures 
    # coverage for highly sheared lattices.  If your lattice is not highly
    # sheared, than you can generate fewer points.
    x_sq = np.arange(-nx, nx, dtype=float)
    y_sq = np.arange(-ny, nx, dtype=float)
    x_sq.shape = x_sq.shape + (1,)
    y_sq.shape = (1,) + y_sq.shape
    # Now shear the whole thing using the lattice vectors
    x_lattice = lattice_vectors[0][0]*x_sq + lattice_vectors[1][0]*y_sq
    y_lattice = lattice_vectors[0][1]*x_sq + lattice_vectors[1][1]*y_sq
    # Trim to fit in box.
    mask = ((x_lattice < image_shape[0]/2.0)
             & (x_lattice > -image_shape[0]/2.0))
    mask = mask & ((y_lattice < image_shape[1]/2.0)
                    & (y_lattice > -image_shape[1]/2.0))
    x_lattice = x_lattice[mask]
    y_lattice = y_lattice[mask]
    # Translate to the centre pix.
    x_lattice += center_pix[0]
    y_lattice += center_pix[1]
    # Make output compatible with original version.
    out = np.empty((len(x_lattice), 2), dtype=float)
    out[:, 0] = y_lattice
    out[:, 1] = x_lattice
    return out

# create the environment for the cells. 
# Just a regular lattice 
env_type = 3
in_tissue = 1
class Environ(object):
    def __init__(self, lx, ly, packing, envtype=2, envbtype=3, centre=np.array([0,0]) ): 
        self.L = (lx, ly)
        self.xbounds = (-lx/2, lx/2)
        self.ybounds = (-ly/2, ly/2)
        self.packing = packing
        self.centre = centre


        self.envtype = envtype
        self.envbtype = envbtype
        self.particles = []

        # the attributes that an environment particle needs 
        self.ekeys = ['id', 'type', 'x', 'y', 'z', 'nx', 'nvx', 'nvy', 'nvz', 'in_tissue']
        # set everything else to zero

    def generate_hex_lattice(self, zpos):
 
        # define the hex lattice vectors
        packing = self.packing
        hv1 = [packing, 0.]
        hv2 = [packing/2, np.sqrt(packing**2 - (packing/2)**2)]
        hx = generate_lattice(self.L, [hv1, hv2])
        self.N = len(hx)
        # at the moment CellInput.prepend_env expects a bpos attribute with positions
        self.bpos = np.column_stack([hx, np.full(self.N, zpos)])

        # next we want to know lists of adjacent particles so that we can create bonds
        # and also the boundary particles so that we can set special types

        # going to find adjacencies using the celllist
        from CellList2D import CellList2D
        rcut = hv1[0]*1.1
        cl = CellList2D(self.L, rcut)

        for i, pt in enumerate(hx):
            cl.add_particle(pt, i)
        neighbours = {}
        hxids = np.arange(len(hx))
        for i, pt in enumerate(hx):
            prox = lambda ptt: np.linalg.norm(pt - ptt) <= rcut
            neigh = cl.get_neighbours(pt)
            neigh.remove(i)
            rneigh = [k for k, ptt in enumerate(hx[neigh]) if prox(ptt)]
            neigh = [neigh[k] for k in rneigh]
            neighbours[i] = neigh

        # create a set of pairs
        hpairs = set()
        for i, neigh in neighbours.items():
            for j in neigh:
                hpairs.add(frozenset([i, j]))

        btype = 1
        bonds = [[idx, btype, i, j] for idx, (i, j) in enumerate(list(hpairs))]
        self.bonds = bonds

        # Ok now find the extrema and set their types. Also construct particles.
        tparr = np.full(self.N, self.envtype, dtype=int)
        # overwrite elements of this array that are on the boundary with self.envbtype
        lxlims = np.min(hx[:,0]) + hv1[0]/2. , np.max(hx[:,0]) - hv1[0]/2.
        lylims = np.min(hx[:,1]) + hv2[1]*0.1, np.max(hx[:,1]) - hv2[1]*0.1

        boundaryn = 0
        for i, pt in enumerate(hx):
            x,y = pt
            if x <= lxlims[0] or x >= lxlims[1]  \
             or y <= lylims[0] or y >= lylims[1]:
                tparr[i] = self.envbtype
                boundaryn += 1

        print 'no. boundary substrate:', boundaryn
        print 'total number in substrate:', len(hx)
        self.tparr = tparr
        return hx, tparr, bonds

    def dump_bonds(self, out):
        # naming convention <inputname>.bonds
        out = path.splitext(out)[0] +'.bonds'
        with open(out, 'w') as f:
            for bond in self.bonds:
                f.write( ('{:5d} '*4 + '\n').format(*bond) )


    def generate_square_box(self):
        # square lattice
        lx, ly =self.L
        self.nx = int(round(self.packing * lx))
        self.ny = int(round(self.packing * ly))
        self.N = self.nx * self.ny

        xline = np.linspace(self.xbounds[0], self.xbounds[1], self.nx, True)
        yline = np.linspace(self.ybounds[0], self.ybounds[1], self.ny, True)

        xaxis = np.zeros(len(xline))
        yaxis = np.zeros(len(yline))
        bottom = np.column_stack([xline, xaxis+yline[0]]) 
        top = np.column_stack([xline, xaxis+yline[-1]]) 
        left = np.column_stack([yaxis+xline[0], yline])[1:-1] # exclude corner points
        right = np.column_stack([yaxis+xline[-1], yline])[1:-1] # exclude corner points

        right = right[::-1]; bottom = bottom[::-1] # orient the boundary clockwise

        self.sides = OrderedDict(zip(['top', 'right', 'bottom', 'left'], [top, right, bottom, left]))
        self.bpos = np.concatenate([ top , right, bottom, left ] )
        self.bpos = np.column_stack([self.bpos[:,0], self.bpos[:,1], np.full(self.bpos.shape[0], 0.)])

    # delete 2n + 1 environment particle from the right edge of the box
    # run after generate_square_box()
    def breakbox(self, n=4):
        top, right, bottom, left = self.sides.values()
        dex = len(right)/2

        tdel = np.arange(dex-n, dex+n+1, 1, dtype=int)
        tcorner, bcorner = right[tdel[0]], right[tdel[-1]]

        liplen = 10.
        lipn = self.packing*liplen
        tlip = np.column_stack([np.linspace(tcorner[0], tcorner[0]+liplen, lipn), np.full(lipn, tcorner[1])])
        blip = np.column_stack([np.linspace(bcorner[0], bcorner[0]+liplen, lipn), np.full(lipn, bcorner[1])])


        right = np.delete(right, tdel, axis=0)
        #print map(len, [top, right, bottom, left])

        self.bpos = np.concatenate([ top , right, bottom, left, tlip, blip ] )
        self.bpos = np.column_stack([self.bpos[:,0], self.bpos[:,1], np.full(self.bpos.shape[0], 0.)])

           
### global
# preferred order
pref = ['id', 'type', 'boundary', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'nx', 'ny', 'nz', 'nvx', 'nvy', 'nvz', 'area']

# This class operates on input files for Samos. Including adding environments to the system.
class CellInput(object):
    def __init__(self, ifile):
        # Read any input file to start with
        self.ifile = ifile
        # The ReadData object is a bit clunky so convert to an ordered dictionary
        self.outd = rwoperate(ifile)
        # The number of particles
        self.N = len(self.outd.values()[0])
        # keep track of the tissue ids
        self.all = None
        if 'id' in self.outd:
            self.all = self.outd['id']
        else:
            print 'Warning: This file has no ids, assume range starting form zero'
            self.all = np.arange(len(self.outd['x']))



    # The reason we prepend the environment is because when cells are removed
    #  the particle ids change.
    # We need the particle ids to stay the same for the substrate particles.
    # We will be connecting them with harmonic springs
    def prepend_env(self, env):
        outd = self.outd
        otherattrs = set(outd.keys()) - set(env.ekeys)
        print 'attrs to set to zero', otherattrs
        # create id array
        nbp = len(env.bpos)
        startid = len( outd.values()[0] )
        idx = np.arange(startid, startid + nbp, 1, dtype=int)
        x = env.bpos[:,0]
        y = env.bpos[:,1]
        z = env.bpos[:,2]
        zero = np.zeros(nbp)
        zeroi = np.zeros(nbp, dtype=int)
        one = np.full(nbp, 1.)
        # not all environment particles have the same type
        if hasattr(env, 'tparr'):
            envt = env.tparr
        else:
            envt = np.full(nbp, env_type, dtype=int)


        # match these up
        #self.ekeys = ['id', 'type', 'x', 'y', 'z', 'nz', 'nvx', 'nvy', 'nvz', 'in_tissue']
        ekarr = [idx, envt, x, y, z, one, zero, zero, one, zeroi]

        for i, k in enumerate(env.ekeys):
            outd[k] = np.concatenate([ekarr[i], outd[k]])

        for att in otherattrs:
            outd[att] = np.concatenate([zeroi, outd[att]])

        # fix ids
        outd['id'] = np.arange(len(outd['id']))

        # keep track the environment particle ids
        self.envids = idx


    # generic set method
    def set(self, param, val):
        self.outd[param] = np.full(self.N, val)

    # fix the direction vector (nx, ny, nz) to be non-zero and have length 1.
    # for environment
    def fixn(self):
        assert hasattr(self.outd, 'in_tissue')
        in_tissue = self.outd['in_tissue']
        nx = self.outd['nx']
        ny = self.outd['ny']
        nz = self.outd['nz']
        for i in range(len(in_tissue)):
            if np.linalg.norm(np.array([nx[i], ny[i], nz[i]])) == 0.:
                nz = 1. # choose e_z for the direction

    # keep this corner for methods which return selections of the data
    def circle_fetch(self, rad=5.15, condition='inside'):
        # fetch the ids of all the particles 
        outd = self.outd
        pts = np.column_stack([outd['x'], outd['y']])
        condition = lambda x,y:x<y if 'inside' else lambda x,y:x>y
        ids = []
        for idx, pt in zip(outd['id'], pts):
            if condition(lg.norm(pt), rad):
                ids.append(idx)
        return np.array(ids)

    def sample(self, init, perc):
        types = self.outd['type']
        ids = self.outd['id']
        targets = []
        for i, tp in zip(ids, types):
            if tp == init:
                targets.append(i)
        ntc = int(perc*len(targets))
        print ntc
        print 'changing the type of {} particles'.format(ntc)
        return rn.sample(targets, ntc)

    def get_sides(self):
        outd = self.outd
        left = []
        right = []
        for i, idx in enumerate(self.outd['id']):
            x = outd['x'][i]
            if outd['boundary'][i]: continue # don't touch boundary guys
            if x < 0.:
                left.append(idx)
            else:
                right.append(idx)
        return left, right
    
    def select_type(self, tp):
        types = self.outd['type']
        ids = self.outd['id']
        # I got fed up with np.where and np.argwhere 
        return np.array( [idx for idx in ids if types[idx] == tp] )

    def set_types(self, ids, tptarget):
        # set the types of particles [ids] to tptarget

        # Don't assume that particle ids are a range starting from 0
        base = np.array(self.outd['id'])[ids] # numpy indexing
        types = self.outd['type']
        for i in base:
            types[i] = tptarget
        self.outd['type'] = types

    def add_boundary_type(self):
        # create a type column 
        cell_type= 2
        types = np.full(self.N, cell_type, dtype=int)
        for i, boundary  in enumerate(self.outd['boundary']):
            if boundary != 0:
                types[i] = 1
        self.outd['type'] = types

    def add_radius(self, r=1., renv=2.):
        rads = np.zeros(self.N)
        for i, idx in enumerate(self.outd['id']):
            rads[i] = r if self.outd['in_tissue'][i] == 1 else renv

        self.outd['radius']= rads
        
    def make_input(self):
        # A bit different. Take a .dat file and add nvz column. 
        # This should make it a valid input file.
        nvx, nvy = np.zeros(self.N), np.zeros(self.N)
        nvz = np.full(self.N, 1.)
        self.outd['nvx'] = nvx
        self.outd['nvy'] = nvy
        self.outd['nvz'] = nvz
        # also get rid of some data that doesn't belong in the input column
        todel = ['fx', 'fy', 'fz', 'cell_area', 'cell_perim', 'cont_num']
        for key in todel:
            if key in self.outd:
                del self.outd[key]

    def make_spherical_input(self):
        outd = self.outd
        nvx = np.zeros(self.N)
        nvy = np.zeros(self.N)
        nvz = np.zeros(self.N)
        x, y, z = outd['x'], outd['y'], outd['z']
        for i, xyz in enumerate(zip(x, y, z)):
            xt, yt, zt = xyz
            ve = np.array([xt, yt, zt])
            nvx[i], nvy[i], nvz[i] =  ve/norm(ve)
        self.outd['nvx'] = nvx
        self.outd['nvy'] = nvy
        self.outd['nvz'] = nvz

    def sphere_prepend(self):
        #needs to be added to the start of the dict so great a new object
        noutd = OrderedDict()
        ids = np.arange(self.N)
        radius = np.full(self.N, 1.)
        tp = np.full(self.N, 1, dtype='int32')
        noutd['id'] = ids
        noutd['type'] = tp
        noutd['radius'] = radius
        for k, v in self.outd.items():
            noutd[k] = v
        self.outd = noutd

    def boundary_flags(self, boundaries):
        bd = self.outd['boundary']
        for b_id, boundary in enumerate(boundaries):
            # label boundaries starting from 1, cells are 0
            b_id += 1   
            for i in boundary:
                bd[i] = b_id
        self.outd['boundary'] = bd

    def addintissue(self):
        self.outd['in_tissue'] = np.full(self.N, 1, dtype=int)

    def addenv(self, environ):
        self.make_input()
        for p in environ.particles:
            pass # todo

    def point_outwards(self):
        centre = [0., 0.]
        cx, xy = centre
        bd = self.outd['boundary']
        allbd = np.where(bd==1)[0]
        xy = zip(np.array(self.outd['x'])[allbd], np.array(self.outd['y'])[allbd])
        for i, pt in zip(allbd, xy):
            x, y = pt[0] - cx, pt[1] - cy
            norm = np.sqrt(x**2 + y**2)
            nx, ny = x/norm, y/norm
            self.outd['nx'][i] = nx
            self.outd['ny'][i] = ny

    def add_bflag(self):
        boundary = np.zeros(self.N, dtype=int)  
        for idx, tp in zip( self.outd['id'], self.outd['type'] ):
            boundary[idx] = 1 if tp == 1 else 0
        self.outd['boundary'] = boundary



    # Finally. A consistent output routine
    def dump(self, fout=None):

        ifile = self.ifile
        os.system(' '.join(['mv', ifile, ifile+'.save']))
        if fout: ifile = fout
        with open(ifile, 'w') as fo:
            io.dump(self.outd, fo, htag='keys:')


    #Fix obstacle
    def orient_xy(self):

        direction = setup_for_obstacle(self)
        r, n, rot_angle = rotate_xyz(self, direction)

        # now pick an id which is on the x-y plane
        print np.absolute(r[:,2])
        seedid = np.argmin(np.absolute(r[:,2]))
        print 'using seed particle in the band', seedid, 'at position', self.rvals[seedid]
        
        #hardcode constraint radius
        R = 28.2094791774 
        #geom = Geometry.GeometrySphere(param(R))
        rcut = 2.4
        # create celllist and Populate it with all the particles:
        # This is just for greater efficiency
        #clist = CellList(geom, rcut)
        #x, y, z = self.outd['x'], self.outd['y'], self.outd['z']
        #for i, rval in enumerate(rvals):
            #self.clist.add_particle(rval,i)

        seed = self.rvals[seedid]
        return seed, r, n, rot_angle
   
    
    def whatobstacle(self, seed, obsr):

        ids = self.outd['id']
        # define boolean function
        def isnear(idx):
            r = self.rvals[idx]
            if norm(seed - r) < obsr:
                return 1
            return 0

        obst = []
        for idx in ids:
            if isnear(idx):
                obst.append(idx)
        print 'obstacle particles'
        print obst

        return obst 
 

    def fixobstacle(self, obsr):
        seed, r, n, rot_angle = self.orient_xy(self)
        obst = self.whatobstacle(seed, obsr)
        #now just change the type of those particles to fix them
        for idx in obst:
            self.outd['type'][idx] = 2

    def insertobstacle(self, obsr):
        # add 0.5 to obsr and create space for obstacle
        seed, r, n, rot_angle = self.orient_xy()
        # or ... choose (R, 0, 0) for the obstacle
        R = 28.2094791774 
    
        # throw away the seed we have
        seed = np.array([R, 0, 0])
        rseed = rotate_vectorial(np.array([seed]), np.array([n[0]]), rot_angle)[0]
        print rseed
        obst = self.whatobstacle(rseed, obsr+0.5)

        # remove these particles
        mask = np.zeros(self.N, dtype='bool')
        mask[np.array(obst)] = 1
        for k, v in self.outd.items():
            # so all the data is in numpy array form now
            # we should make it like this at the read-in step
            self.outd[k] = np.array(v)[~mask]
        
        # we didn't rotate the velocities and directors only the positions
        # Best to construct the obstacle and then rotate it 
        ll = len(self.outd.values()[0])
        print 'number of particles remaining', ll

        # add obstacle particles
        # needed fields for adding a particle
        # [id, type, radius, x, y, z, vx, vy, vz, nx, ny, nz, nvx, nvy, nvz]
        obs = Obstacle(seed, obsr)
        obsd = obs.build_particles(R, n, rot_angle)


        # append these ordered dictionaries
        for k in self.outd:
            self.outd[k] = np.append(self.outd[k], obsd[k])

        # reset ids
        self.N = len(self.outd['id'])
        self.outd['id'] = np.arange(self.N, dtype='int32')
        print self.outd['id']

        

### Obstacle operations

def resettypes(ifile):
    ci= CellInput(ifile)
    ci.outd['type'] = np.full(ci.N, 1, dtype=int)
    ci.dump(ifile)

def insertobstacle(ifile, obsr=5.):
    ci = CellInput(ifile)
    ci.insertobstacle(obsr)
    ci.dump(ifile)

# This is the one that adds ids and types, etc. to Rastko data file
def sphereinput(ifile):
    Ci = CellInput(ifile)
    Ci.sphere_prepend()
    Ci.make_spherical_input()
    Ci.dump()

def fixobstacle(ifile, obsr=5.):
    Ci = CellInput(ifile)
    Ci.fixobstacle(obsr)
    Ci.dump()


### Environment operations
def addsquareenv(ifile, dens=1.5):
    ci = CellInput(ifile)
    # add the in_tissue flag
    ci.addintissue()
    # Create the enivronment
    env = Environ(60., 60, dens)
    env.generate_square_box()
    #outd = env.append(ci.outd)
    ci.prepend_env(env)

    ci.dump(ifile)

# A broken box
def addbrokensquare(ifile, dens=2.5, n=10):
    ci = CellInput(ifile)
    # add the in_tissue flag
    ci.addintissue()
    # Create the enivronment
    env = Environ(40., 40, dens)
    env.generate_square_box()
    env.breakbox(n)
    #outd = env.append(ci.outd)
    ci.prepend_env(env)

    ci.dump(ifile)

def addsubstrate(ifile):
    packing = 1.4
    ci = CellInput(ifile)
    ci.addintissue()
    env = Environ(65., 65., packing)
    zpos = -0.5
    env.generate_hex_lattice(zpos)
    ci.prepend_env(env)

    env.dump_bonds(ifile)
    ci.dump(ifile)

def fixn(ifile):
    ci = CellInput(ifile)
    ci.fixn()
    ci.dump()

### Basic operations

# very basic add back radius to the input with all particles radius of 1.
def addradius(ifile):
    Ci = CellInput(ifile)
    Ci.add_radius()
    Ci.dump()

# sample the particles of a given type and change their types 
def sampletype(ifile, init=2, tptarget=3, perc=0.5):
    Ci = CellInput(ifile)
    cids = Ci.sample(init, perc)
    Ci.set_types(cids, tptarget)
    Ci.dump()

# Use the CellInput object to quickly construct the operation I need
def coretype(ifile, tp=3):
    Ci = CellInput(ifile)
    cids = Ci.circle_fetch()
    Ci.set_types(cids, tptarget=tp)
    Ci.dump()

# set the types of all the particles on one side
def sidetype(ifile):
    Ci = CellInput(ifile)
    left, right = Ci.get_sides()
    Ci.set_types(left, 2)
    Ci.set_types(right, 3)
    Ci.dump()

def addtype(ifile):
    Ci = CellInput(ifile)
    Ci.add_boundary_type()
    Ci.dump()

def onetype(ifile):
    Ci = CellInput(ifile)
    Ci.set_types(Ci.all, 1)
    Ci.dump()
    
def makeinput(ifile, name=None):
    if not name: name = ifile
    Ci = CellInput(ifile)
    Ci.make_input()
    #Ci.add_radius()
    #Ci.addintissue()
    Ci.dump(name)

def pointout(ifile):
    Ci = CellInput(ifile)
    Ci.point_outwards()
    Ci.dump()

def settype(ifile, tp=2, target=4):
    ci = CellInput(ifile)
    tpi = ci.select_type(tp)
    ci.set_types(tpi, target)
    ci.dump()

def swapactivetypes(ifile, tp_1=2, tp_2=3):
    Ci = CellInput(ifile)
    tpid1 = Ci.select_type(tp_1)
    tpid2 = Ci.select_type(tp_2)
    Ci.set_types(tpid1, tp_2)
    Ci.set_types(tpid2, tp_1)
    Ci.dump()

def addbflag(ifile):
    Ci = CellInput(ifile)
    Ci.add_bflag()
    Ci.dump()


def piarea(ifile):
    Ci = CellInput(ifile)
    if 'area' in Ci.outd:
        Ci.set('area', round(np.pi, 6))
    else:
        Ci.outd['area'] = np.full(Ci.N, round(np.pi, 6))
    Ci.dump()

# old style method (doesn't use the CellInput class)
# extract_boundary
# Method to extract a .boundary file form a FRESHLY GENERATED .input file
#  where we assume that the boundary particles are all ordered and placed at the end of the file.
# If you want to extract a boundary from an simulation output for restarting simulations, you need to run
# analyse_cells -i <input_file_name>.input mesh --boundary
# this will produce a file named <input_file_name>.boundary
def eb(ifile):
    outd = rwoperate(ifile)
    idx = outd['id']
    boundary = np.array(outd['boundary'], dtype=int)
    bdids = np.where(boundary==1)[0]
    #idids = idx[bdids]
    idids = bdids
    lines = []
    for newid, (i, j) in enumerate(zip(idids, np.roll(idids,-1))):
        lines.append([newid, i, j])

    # write out
    base, _ = path.splitext(ifile)

    def write_boundary(fo, lines):
        # write the header
        fo.write('#\n')
        for line in lines:
            fo.write(' '.join(map(str, line)) + '\n')

    with open(base+'.boundary', 'w') as fo:
        write_boundary(fo, lines)


#####################################################################################
# Non object oriented code

# Triangulation needs positions and associated areas 
def makeout(rvals, parea):
    outd = OrderedDict()
    keys = ['id', 'x', 'y', 'z', 'area']
    for k in keys:
        outd[k] = []

    for i, r in enumerate(rvals):
        x, y, z = r
        outd['x'].append( x )
        outd['y'].append( y )
        outd['z'].append( z ) 
        outd['id'].append( i )
        outd['area'].append( parea[i] )
    return outd

npnormal= np.array([0.,0.,1.])
def fillout(outd, normal=npnormal):
    ld= len(outd.values()[0])
    outd['nvz'] = np.full(ld, 1.)
    for hh in ['vx', 'vy', 'vz', 'nvx', 'nvy']:
        outd[hh] = np.full(ld, 0.)
    return outd

def _extend_ptset(pts, npts):
    rpts = []
    for npt in npts:
        isnew = True
        for pt in pts:
            if np.allclose(np.array(pt), np.array(npt), rtol=0.1):
                isnew = False
                break
        if isnew:
            rpts.append(npt)
    return pts.extend(rpts)

# minimal hexagonal configuration
def areahex(area=3.0):
    # first find the radius length of a hexagon with the correct area
    l = np.sqrt((2*area)/(3*np.sqrt(3)))
    # distance between hex centers
    lc = np.cos(np.pi/6) * 2 * l
    rads = np.linspace(0, 2*np.pi, 6, False)
    def hexes(rad):
        return [np.cos(rad), np.sin(rad), 0.]
    hpts = lc*np.array(map(hexes, rads))
    pts = list(hpts)

    R = wr.rotation_2d(np.pi/4)
    for pt in pts:
        pt[1] *= 2
        pt[:2] = np.einsum('mn,n->m', R, pt[:2])

    crads = lc*np.array(map(hexes, rads+(np.pi/6)))
    for hcenter in hpts:
        _extend_ptset(pts, [hcenter+hpt for hpt in hpts])
    outd = makeout(pts, [area]*len(pts))


    outd =fillout(outd)
    return outd

# Single cell, n sides
def single(nsides, prefarea=5., origin=[0.,0.,0.]):
    parea = prefarea
    rvals = np.array(origin).reshape((1,3))
    rad = np.linspace(-np.pi, np.pi, nsides, endpoint=False)
    nb = len(rad)
    x = np.cos(rad)
    y = np.sin(rad)
    z = np.zeros(nb)
    area = np.full(nb+1, parea)
    rbound = np.column_stack([x, y, z])
    rvals = np.vstack([rvals, rbound])
    #print np.column_stack([x, y, z])

    outd = makeout(rvals, area)
    outd = fillout(outd)
    return outd
single.defaults = ['6']
single.call = ['single( nsides )']


# Manipulating existing input files

def random_direction():
    rtheta = random()* 2*np.pi
    return np.cos(rtheta), np.sin(rtheta)

# read the dat files and then adjust the preferred areas to a uniform random distrubution
# decpreciated
def randomise_area(ifile, minmax=[3.0, 3.5]):
    mn,mx = minmax
    def rarea(outd, i):
        # set the area
        uarea = uniform(mn, mx) 
        #print 'setting area ', uarea
        outd['area'][i] = uarea
    outd = rwoperate(ifile, rarea)

    with open('areaset.input', 'w') as fo:
        io.dump(outd, fo, htag='keys:')

def shift_area(ifile, shift):
    def fshift(outd, i):
        outd['area'][i] += shift
    outd= rwoperate(ifile, fshift)
    with open('shifted.input', 'w') as fo:
        io.dump(outd, fo, htag='keys:')


# of the particles of type x randomly choose a percentage of them to be of type y
def addactive(ifile, init=2, perc=0.1, fin=3):
    outd = rwoperate(ifile, None)
    types = map(int, outd['type'])
    linit = [i for i, ty in enumerate(types) if ty==init]
    llinit = len(linit)
    tochange = int(perc * llinit)
    print 'changing the type of {} particles'.format(tochange)

    change = rn.sample(range(llinit), tochange)
    for i in change:
        types[i] = fin
    outd['type'] = types

    with open('chtyped.input', 'w') as fo:
        io.dump(outd, fo, htag='keys:')
        
def addarea(ifile, area=3.0):
    outd = rwoperate(ifile)
    k1 = outd.keys()[0]
    areas = np.zeros(len(outd[k1]))
    areas.fill(area)
    outd['area']  = areas

    os.system(' '.join(['mv', ifile, ifile+'.save']))
    with open(ifile, 'w') as fo:
        io.dump(outd, fo, htag='keys:')

# cut a piece of an existing input dat file
#constructmesh.py circleslice epithelial_equilini.dat 
def circleslice(ifile, radius=5.15,  condition='square', boundary=True,  center=[0.,0.,0.]):
    rdat = ReadData(ifile)
    keys = rdat.keys
    x = np.array(rdat.data[keys['x']])
    y = np.array(rdat.data[keys['y']])
    z = np.array(rdat.data[keys['z']])
    rlists = [rdat.data[k] for k in keys.values()]
    outd = OrderedDict()
    for k in keys:
        outd[k] = []
    for i, xval in enumerate(x):
        cd= norm(center - np.array([x[i],y[i],z[i]]))
        if condition=='circle':
            if cd > radius:
                continue
        elif condition=='annulus':
            if cd < radius:
                continue
        elif condition=='square':
            xpt, ypt = x[i], y[i]
            if abs(xpt) > radius or abs(ypt) > radius:
                continue

        # keep this line
        for kn, kv in keys.items():
            outd[kn].append( rdat.data[kv][i] )
    # fix the ids
    ids = range(len(outd.values()[0]))
    lz = ids[-1] +1
    outd['id'] = ids
    # set all the velocities to zero
    for kn in ['vx', 'vy']:
        outd[kn] = list(np.zeros(lz))
    if 'type' in outd:
        outd['type'] = map(int, outd['type'])
    outd['boundary'] = map(int, outd['boundary'])

    # add circular boundary of cells
    if boundary:
        rr = radius + 0.7
        bcelldensity = 0.9
        nbcells = round( (2* np.pi * rr) * bcelldensity)
        thetas = np.linspace(0, 2*np.pi, nbcells, False)
        print 'number of boundary particles', len(thetas)

        xb = rr *np.cos(thetas)
        yb = rr *np.sin(thetas)
        area = rdat.data[keys['area']][0]
        for istepped, xval in enumerate(xb):
            i = lz + istepped
            outd['id'].append(i)
            outd['radius'].append(1.)
            outd['x'].append(xb[istepped])
            outd['y'].append(yb[istepped])
            outd['nvz'].append(1.)
            outd['area'].append(area)
            # make a random direction
            nx, ny = random_direction()
            outd['nx'].append(nx)
            outd['ny'].append(ny)
            outd['type'].append(1)
            outd['boundary'].append(1)
            for kn in ['z', 'vx', 'vy', 'vz', 'nz', 'nvx', 'nvy']:
                outd[kn].append(0.)

    # dump
    with open('slice.input', 'w') as fo:
        io.dump(outd, fo)

if __name__=='__main__':

    args = sys.argv[1:]
    nargs = len(args)
    if nargs is 0:
        print 'python constructmesh <function> [*arguments]'
        sys.exit()

    fname = args[0]
    fargs = args[1:]
    try:
        f_using = locals()[fname]
    except KeyError:
        print 'I don\'t have a function \'%s\'', fname
        raise

    if len(fargs) is 0:
        print 'No arguments given to constructmesh.%s' % fname

        if hasattr(f_using, 'defaults'):
            fargs = f_using.defaults
            print 'Using defaults ', f_using.defaults
        else:
            fargs = []

    ff = fargs
    for i, f in enumerate(ff):
        try:
            ff[i] = eval(f)
        except:
            # This string has a '.' and look like an object
            #ff[i] = f
            pass
    print ff

    outd = f_using(*ff)
    if outd:
        with open('test.dat', 'w') as fo:
            io.datdump(outd, fo)



