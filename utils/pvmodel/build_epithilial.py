# * *************************************************************
# *
# *   Soft Active Mater on Surfaces (SAMoS)
# *
# *   Author: Rastko Sknepnek
# *
# *   Division of Physics
# *   School of Engineering, Physics and Mathematics
# *   University of Dundee
# *
# *   (c) 2013, 2014
# *
# *   School of Science and Engineering
# *   School of Life Sciences
# *   University of Dundee
# *
# *   (c) 2015
# *
# *   Author: Silke Henkes
# *
# *   Department of Physics
# *   Institute for Complex Systems and Mathematical Biology
# *   University of Aberdeen
# *
# *   (c) 2014, 2015
# *
# *   This program cannot be used, copied, or modified without
# *   explicit written permission of the authors.
# *
# * ***************************************************************

# Utility code for building random initial configuration on xy plane
# This code was spawned from a configuation script called plane_circle_epithelial.py
#  which was written by Silke Henkes.

# I expanded it to cover producing one or two patches of cells with square boundaries.
# in addition to it's original function of producing a circular patch.

# In order to produce all of my cell configuration files I then operated on the output from 
# this program using Samos to relax the random configuration and parts of /utils/pvmodel/constructmesh.py 
# to modify cell types, areas, etc. as well as to add environment particles to the input file.

# When it comes to obtaining the boundary connectivity information (.boundary file) which cell simulations
# require, the command: 
#   constructmesh.py eb <inputfile>
# Will produce <inputfile>.boundary ONLY for single square and circular patches outputted by this program.
# In order to obtain the boundary information for more complicated cases we can try and triangulate
# using numpy, construct the halfedge data structure and then recover the boundary.
# The command is:
#   analyse_cells -i <inputfile> --triangulate mesh --boundary
# (mesh is a subcommand which changes the behaviour of the program, --triangulate is just a flag.)
# This won't work for patches with multiple boundaries or multiple patches.

# --Daniel Barton 11/10/16


from datetime import *
from random import uniform
from math import *
import argparse
import numpy as np
import sys
import os

pi = np.pi


from particle import *

cell_type= 2
boundary_type = 1

class Plane:

    def __init__(self, Lx, Ly, N, v, bpacking):
        self.L = (Lx,Ly)
        self.N = N
        self.buff = 0.5
        self.bpacking = bpacking
        self.Nbound = int(round( (Lx + self.buff) * (Ly + self.buff)/(1/bpacking) ))
        self.particles = [Particle(i) for i in range(N)]
        self.__generate_pos()
        self.__generate_vel(v)
        self.__generate_director()
        self.__generate_posbound()
        self.Nt = self.N + len(self.boundparticles)

    def __generate_pos(self):
        for i in range(self.N):
            x = uniform(-0.5*self.L[0],0.5*self.L[0])
            y = uniform(-0.5*self.L[1],0.5*self.L[1])
            z = 0
            self.particles[i].r = [x,y,z]

    def __generate_posbound(self):
        lx, ly = self.L
        lx, ly = lx/2., ly/2.

        Nlinex = int(round( (lx + 2*self.buff) * (1/self.bpacking) ) )
        Nliney = int(round( (ly + 2*self.buff) * (1/self.bpacking) ) )
        
        rvals = []
        xline = np.linspace(-lx-self.buff, lx+self.buff, Nlinex, True)
        yline = np.linspace(-ly-self.buff, ly+self.buff, Nliney, True)
        xaxis = np.zeros(Nlinex)
        yaxis = np.zeros(Nliney)
        bottom = np.column_stack([xline, xaxis+yline[0]]) 
        top = np.column_stack([xline, xaxis+yline[-1]]) 
        left = np.column_stack([yaxis+xline[0], yline])[1:-1] # exclude corner points
        right = np.column_stack([yaxis+xline[-1], yline])[1:-1] # exclude corner points

        right = right[::-1]; bottom = bottom[::-1] # orient the boundary clockwise

        bparticles = np.concatenate([ top , right, bottom, left ] )

        self.boundparticles = [Particle(i) for i in range(self.N,self.N+len(bparticles))]
        for i in range(len(bparticles)):
            x, y = bparticles[i]
            z = 0.
            self.boundparticles[i].r = [x, y, z]

            sval = uniform(0,2*pi)
            self.boundparticles[i].n=[cos(sval),sin(sval),0]


    def __generate_vel(self,vav=1.0):
        for p in self.particles:
            phi = uniform(0,2*pi)
            p.v = [vav*cos(phi),vav*sin(phi),0.0]

    def __generate_director(self):
        for p in self.particles:
            phi = uniform(0,2*pi)
            p.n = [cos(phi),sin(phi),0.0]


    def writefo(self,out, header=True):
        gentime = datetime.now()

        headers = ['id', 'type', 'radius', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'nx', 'ny', 'nz', 'nvx', 'nvy', 'nvz', 'boundary']
        htypes = [int, int, float, float, float, float, float, float, float, float, float, float, float, float, float, int]
        dataline = ''
        for ht in htypes:
            if ht == int: dataline += '{:>10}'
            elif ht == float: dataline += '{:11.5f}'
        dataline += '\n'
        if header:
            #out.write('# Total of %s particles\n' % str(self.N+self.Nbound))
            #out.write('# Generated on : %s\n' % str(gentime))
            out.write('# ' + ('{:<8}  '*len(headers)).format(*headers) + '\n')

        for p in self.particles:
            p.tp = cell_type
            boundary = 0
            x, y, z = p.r
            vx, vy, vz = p.v
            nx, ny, nz = p.n
            nvx, nvy, nvz = [0.,0.,1.]
            out.write(dataline.format(*[p.idx,p.tp,p.R,x,y,z,vx,vy,vz,nx,ny,nz,nvx,nvy,nvz,boundary]))
        for p in self.boundparticles:
            p.tp = boundary_type
            p.R = 2.0
            boundary = 1
            x, y, z = p.r
            vx, vy, vz = p.v
            nx, ny, nz = p.n
            nvx, nvy, nvz = [0.,0.,1.]
            out.write(dataline.format(*[p.idx,p.tp,p.R,x,y,z,vx,vy,vz,nx,ny,nz,nvx,nvy,nvz,boundary]))
    
    def translate(self, pt):
        x, y = pt
        # move the x,y positions of all the particles by pt
        for p in self.particles:
            p.r = [p.r[0]+x, p.r[1]+y, p.r[2]]
        for p in self.boundparticles:
            p.r = [p.r[0]+x, p.r[1]+y, p.r[2]]

    def shiftids(self, shift):
        for p in self.particles:
            p.idx += shift
        for p in self.boundparticles:
            p.idx += shift

    def write(self, outfile):
        with open(outfile, 'w') as fo:
            self.writefo(fo)

    # boundary representation
    # [edge_id, p_id_1, p_id_2]
    # assume that self.boundparticles is ordered
    def bdrep(self, ishift=0, idxshift=0):
        pids = [p.idx for p in self.boundparticles]
        zpids = zip(range(len(pids)), pids, np.roll(pids,-1))
        return [[i+ishift, pi+idxshift, pj+idxshift] for i, pi, pj in zpids]

# need a way to combine two planes together
# How about a general way of combining multiple configurations
# time to generalise

def combinewrite(outfile, psets, ptt=[(-10,0),(10,0)] ):
    # translate all sets
    for i, ps in enumerate(psets):
        ps.translate(ptt[i])
    #
    lastid = 0
    # boundary files
    bdls = []
    with open(outfile, 'w') as out:
        last = 0
        for i, seti in enumerate(psets):
            # update ids
            seti.shiftids(last)
            if i == 0:
                seti.writefo(out, header=True)
            else:
                seti.writefo(out, header=False)
            bdls.append( seti.bdrep(last) )
            last += seti.Nt

    bdfile = os.path.splitext(outfile)[0] + '.boundary'
    def writebd(bd):
        form = '{:<5d} {:>5d} {:>5d}\n'
        with open(bdfile, 'w') as bdo:
            for line in bd:
                bdo.write(form.format(*line))
    
    # finally need to concatenate the lists bda and bdb
    writebd(sum(bdls, []))


class Circle(Plane):
    def __init__(self, R, N, v, poly,boundary,bpacking):
        self.R = R
        self.N = N
        self.poly = poly
        self.boundary=boundary
        self.Nbound = int(bpacking*2.0*pi*R/2.0) # Dense boundary wall, linear packing fraction 1.5
        self.particles = [Particle(i,1) for i in range(N)]
        self.boundparticles = [Particle(i,2) for i in range(N,N+self.Nbound)]
        self.__generate_posinside()
        self.__generate_posbound()
        self.__generate_vel(v)
        self.__generate_director()

    def __generate_posinside(self):
        for i in range(self.N):
            rpart=4*self.R**2
            while rpart>(self.R-1.0)**2:
                x = uniform(-self.R+1.0,self.R-1.0)
                y = uniform(-self.R+1.0,self.R-1.0)
                z = 0
                rpart=x**2+y**2
            self.particles[i].r = [x,y,z]
            self.particles[i].R = uniform(1-0.5*self.poly,1+0.5*self.poly)

    def __generate_posbound(self):
        for i in range(self.Nbound):
            # fractional position along circumference (angle)
            sval = 1.0*i/self.Nbound*2*pi
            self.boundparticles[i].r=[self.R*cos(sval),self.R*sin(sval),0]
            # introduce a little noise there; no ratchets here
            self.boundparticles[i].R = uniform(1-0.5*self.poly,1+0.5*self.poly)

    def __generate_vel(self,vav=1.0):
        for p in self.particles:
            phi = uniform(0,2*pi)
            p.v = [vav*cos(phi),vav*sin(phi),0.0]

    def __generate_director(self):
        for p in self.particles:
            phi = uniform(0,2*pi)
            p.n = [cos(phi),sin(phi),0.0]
        # This is actually intriguing: make them point inwards, outwards or random
        for i in range(self.Nbound):
            if self.boundary==2:
                sval = 1.0*i/self.Nbound*2*pi+pi
            elif self.boundary==1:
                sval = 1.0*i/self.Nbound*2*pi
            else:
                sval = uniform(0,2*pi)
            self.boundparticles[i].n=[cos(sval),sin(sval),0]

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    #parser.add_argument("-r", "--R", type=float, default=10.0, help="glued circle radius")
    parser.add_argument("-l", "--L", type=float, default=10.0, help="Length of the planar box")
    parser.add_argument("-lx", type=float, default=10.0, help="Length of the planar box in x")
    parser.add_argument("-ly", type=float, default=10.0, help="Length of the planar box in y")
    parser.add_argument("-n", "--N", type=int, default=1000, help="number of particles")
    parser.add_argument("-p", "--poly", type=float, default=0.0, help="polydispersity fraction")
    parser.add_argument("-f", "--phi",  type=float, default=1.0, help="packing fraction")
    parser.add_argument("-o", "--output", type=str, default='out.input', help="output file")
    parser.add_argument("-v", "--vavr", type=float, default=0.1, help="average velocity")
    parser.add_argument("-b", "--boundary", type=int, default=0, help="boundary n: 0 random 1 outward 2 inward")
    parser.add_argument("-bp", "--bpacking", type=float, default=1.0, help="Boundary packing density")
    parser.add_argument("-a", "--area", type=float, default=np.pi, 
            help="Decide the number of particles N using approx area per particle")
    parser.add_argument("-rounded", action='store_true', help='Round the corners of the squares')
    valid_types= ['Circle', 'Plane']
    parser.add_argument("-type", type=str, default='Circle', help='Shape of patch')
    args = parser.parse_args()

    # overwrite -n argument
    # N for plane
    #args.N = int(round((args.lx * args.ly)/args.area))
    # N for circle
    args.N = int(round(args.L**2 )) # An individual particle has area pi R^2 = 4 pi

    print 'adding {} particles'.format(args.N)

    print
    print "\tActive Particles on Curved Spaces (APCS)"
    print "\tBuilding a glued circle random configuration on the plane"
    print
    print "\tRastko Sknepnek"
    print "\tUniversity of Dundee"
    print "\t(c) 2013"
    print
    print "\tSilke Henkes"
    print "\tUniversity of Aberdeen"
    print "\t(c) 2014"
    print "\t----------------------------------------------"
    print
    #print "\tRadius : ", args.R
    print "\tpolydispersity: ", args.poly


    print "\tPacking fraction : ", args.phi
    print "\tNumber of particles : ", args.N

    #R = sqrt(args.N/args.phi)
    #print "\tRadius : ", R

    print "\tRadius : ", args.L
    print "\tAverage velocity : ", args.vavr
    print "\tBoundary type (0 random 1 outward 2 inward) : ", args.boundary
    print "\tOutput file : ", args.output
    print

    start = datetime.now()
    #p = Circle(R, args.N, args.vavr,args.poly,args.boundary,args.bpacking)
    if args.type == 'Plane':
        p = Plane(args.lx, args.ly, args.N, args.vavr, args.bpacking)
    elif args.type == 'Circle':
        p = Circle(args.L, args.N, args.vavr, args.poly, args.boundary, args.bpacking)

    p.write(args.output)

    end = datetime.now()

    total = end - start

    print
    print "  *** Completed in ", total.total_seconds(), " seconds *** "
    print
