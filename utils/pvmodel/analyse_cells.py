#!/usr/bin/env python

# 

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

from senario import *
import cmplotting as cm


def process_commands():
    topparser = argparse.ArgumentParser()
    # General arguments
    topparser.add_argument("-i", "--input", type=str, default='cell*.dat', help="Input dat file")
    topparser.add_argument("-last", action='store_true', help='just the last input file ')
    topparser.add_argument("-first", action='store_true', help='just the first input file ')
    topparser.add_argument("-c", "--conf", type=str, default='*.conf', help="Input dat file")
    topparser.add_argument("-d", "--dir", type=str, default='tmp/', help="Output directory")
    topparser.add_argument("-p", "--plot", type=str, default='plots/', help="plots directory")
    topparser.add_argument('-step', type=int, default=1, help='Step through time')
    topparser.add_argument('-start', type=int, default=None, help='Start number')
    topparser.add_argument('-stop', type=int, default=None, help='Stop number')
    topparser.add_argument("-a0",  type=float, default=3., 
            help='Give the target area explicitely')
    topparser.add_argument("-bt", '--boundary_type', type=int, default=1, help='boundary type')
    topparser.add_argument("--all", type=str, default=None,
            help='Run some analysis on every sub-directory')
    topparser.add_argument("--triangulate", action='store_true', 
            help='use numpy for triangulation')
    topparser.add_argument("-f", '--faces', type=str, help='specify the faces file')
    topparser.add_argument("-ref", type=str, default=None, 
            help='specify the number of the reference configuration')

    subparsers = topparser.add_subparsers(dest='subcommand')
 
    # Analysis on the mesh, Particularly for calculating stresses
    parser = subparsers.add_parser('stress', help='Analysis involving stress and forces')
    parser.add_argument("-wl", type=float, default=1., help='smoothing length')
    parser.add_argument("-s", action='store_true')
    parser.add_argument("-w", type=str, default='', 
            help='A list containing start finish and step for wl values')
    parser.add_argument('-kf', "--kernel", type=str, default='quartic',
            help = 'the type of smoothing function to use')
    parser.add_argument('--scale', type=float, default=1., 
            help='scale stress ellipses')
    parser.add_argument('-adj', type=int, default=0, 
            help='The number of adjacent rings of cells to include while averaging')

    parser.add_argument('-lt', '--ltransition', type=float, default=0.02, 
            help = 'Threshold length for a t1 transition')
    parser.add_argument('--test', action='store_true')

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


    # Want a sub parser for simple operations on mesh which don't involve stress calculation
    # Examples. Calculating the structure tensor, calculating cell neighbour topology
    # Associated with the 'Property' senario
    meshparser = subparsers.add_parser('mesh', help='Analysis on the mesh')
    meshparser.add_argument('--areas', action='store_true', 
            help = 'plot the distribution of cell areas')
    meshparser.add_argument('--texture', action='store_true', help='Calculate texture tensor')
    meshparser.add_argument('--top', action='store_true', help='Cell neighbour topology')
    meshparser.add_argument('--boundary', action='store_true', help='Extract the boundaries')

    # Want a new sub parser for calculations which don't involve the mesh
    datparser = subparsers.add_parser('dat', help='Simple analysis on the particle data')
    datparser.add_argument('-type', type=int, default=None, help='specify a particular type')

    # just print some information about the shape index
    shapeparser = subparsers.add_parser('shape') 

    # Look at the connectivity of the mesh
    transparser = subparsers.add_parser('trans', help='Tracking or counting transitions')
    
    # Look at the .dat data 
    transparser = subparsers.add_parser('basic', help='basic operations on .dat')


    # just testing
    testparser = subparsers.add_parser('test')

    args = topparser.parse_args()
    # save plots to 'plots/' and plot data to 'plots/data/'
    args.data = path.join(args.plot, 'data/')
    # Just make the directories for plots ahead of time
    cm.result_dirs()

    
    # cleanup
    if hasattr(args, 'exclude'):
        args.exclude = not args.include 

    return args 

# Process the configuration file.
# We want to know:
# (i) What groups are defined.
# (ii) What what 'vp' potentials are defined, are there more than 1? Their attributes.
# (iii) Check how many integrators are defined. 

def configuration(args):
    # find that configuration file
    if '*' in args.conf:
        # probably gave a regular expression to find the configuration file
        confs = glob.glob(args.conf)
        if not confs:
            return args
        print 'found configuration file(s)', confs
        args.conf = confs[0]


    from read_conf import ReadConf
    rc = ReadConf(args.conf)

    # groups
    
    # let this flag determine whether there are multiple types of active cells
    is_active_types = False
    if 'pair_type_param' in rc.key_words: 
        is_active_types= True

    # Conversion between the names I gave the parameters and the SAMoS names
    samparams = ['lambda', 'K', 'gamma']
    myparams = ['L', 'k', 'gamma']
    pmap = dict(zip(samparams, myparams))

    # This will get the default attributes but not necessarily the correct values 
    # which are set by pair_type_param
    attrs = rc.key_words['pair_potential'][0].attributes
    for attr in attrs:
        if attr.name in samparams:
            setattr(args, pmap[attr.name], float(attr.val))

    # Try and determine what the active types are
    # {int type: {vp attributes dictionary}}
    if is_active_types:
        active_types = {}
        for param in rc.key_words['pair_type_param']:
            tp, k, gamma, L = None, None, None, None
            for attr in param.attributes:
                if attr.name == 'type':
                    tp = int(attr.val)
                elif attr.name == 'K':
                    k = float(attr.val)
                elif attr.name == 'gamma':
                    gamma = float(attr.val)
                # shouldn't ever use this one but leave it in for now
                elif attr.name == 'lambda':
                    L = float(attr.val)
            active_types[tp] = {'k':k,'gamma':gamma,'L':L}
        args.active_types = active_types

    # For the lambda parameter we set the values for pairs of types
    lpairs = {'default':args.L}
    if 'pair_param' in rc.key_words:
        for param in rc.key_words['pair_param']:
            if param.name == 'vp':
                # fill up the lpairs dictionary
                # attrs are always type_1, type_2, lambda
                atype_1, atype_2, alambda = param.attributes
                type_1 = int(atype_1.val)
                type_2 = int(atype_2.val)
                lam = float(alambda.val)
                lpairs[frozenset([type_1, type_2])] = lam
                print 'setting lambda to {} for types {} and {}'.format(type_1, 
                        type_2, lam)


            elif param.name == 'soft':
                pass # for now

    args.lpairs = lpairs
    print 'default lambda is {}'.format(lpairs['default'])

    # We never consider the boundary particles as cells anymore.
    args.exclude = True
    
    return args

def shapeindex(args):
    #target perimeter
    try:
        tperim = -args.L/args.gamma
    except ZeroDivisionError:
        tperim = 0.
        print 'warning, target perimeter undefined'
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
    print

if __name__=='__main__':

    args = process_commands()

    # update the args object with information from the configuration file
    configuration(args)
    
    def meshmain():
        shapeindex(args)
        # Original method of running analysis
            
        if not os.path.exists(args.dir):
            os.mkdir(args.dir)
        print_parameters()
        stressrun = Stress_Senario(args)
        for _ in stressrun: pass

    def shapemain():
        shapeindex(args)

    # Use the MSD senario
    def datmain():
        datread = MSD(args)
        # read the whole dataset
        for _ in datread: pass
        timeline, msd = datread.allmsd()
        cm.plotmsd(timeline, msd)

    def meshproperty():
        prop = Property(args)
        for _ in prop: pass

    def transitions():
        if not args.all:
            trans = Transitions(args)
            for _ in trans: pass

            helist = trans.helist
            tt = TransAnalyse(helist, trans.tlist)
            # compare stress and texture tensor
            #tt.structure_alignment(trans.outnum)

        elif args.all:
        # plot the number of transitions
        #Want to be able to iterate over all folders in a directory, run analyse and collect data
            print 'ntt'
            ntt = OrderedDict()
            import opdir as op
            def macro():
                configuration(args)
                shapeindex(args)
                print_parameters()
                trans = Transitions(args)
                for _ in trans: pass
                ntrans = trans.tlist.nt
                ntt[args.sindex] = ntrans
            op.diriterate(macro, args.all)
            print ntt

            cm.nttplot(ntt, log=True)

    def basicmain():
        sen = Basic(args)
        for _ in sen: pass

        timeline = sen.timeline
        nghosts = sen.numghosts
        #cm.plotnghosts(timeline, nghosts)

        # count n cells
        timeline, n_cells = map(int, sen.timeline), sen.n_cells

        cm.growthexp(timeline, n_cells)

        #plt.clf()
        #plt.plot(timeline, sen.avgshape)
        #plt.show()

        #plt.clf()
        #plt.plot(timeline, sen.avgarea)
        #target = 3.
        #plt.plot(timeline, np.full(len(timeline), 3.), linestyle='--')

        #plt.xlabel('timestep')
        #plt.ylabel('average cell area')

        #tnum = [10, 100, 200, 500, 700, 1000]
        #for tn in tnum:
            #print 'timestep {}, n_cells {}'.format(tn*1000, n_cells[tn])



    if args.subcommand == 'stress':
        meshmain()
    elif args.subcommand == 'dat':
        datmain()
    elif args.subcommand == 'shape':
        shapemain()
    elif args.subcommand == 'mesh':
        meshproperty()
    elif args.subcommand == 'trans':
        transitions()
    elif args.subcommand == 'basic':
        basicmain()

    ## additional commands for testing

    def testconfiguration():
        print args.k
        print args.active_types

    if args.subcommand == 'test':
        testconfiguration()
