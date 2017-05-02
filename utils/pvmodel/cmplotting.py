#!/usr/bin/env python

import numpy as np
import ioutils as io
#from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib import rc
from matplotlib import pyplot as plt
import time
from glob import glob
import os.path as path
import os
import pickle
from collections import OrderedDict, Sequence
import json

# default font size
mpl.rcParams.update({'font.size': 18, 'legend.fontsize':10})

# default axis
rc('text', usetex=True)
ax = plt.gca()

# helper for defaultsave
pdir = 'plots/'
ddir = os.path.join(pdir, 'data/')
def result_dirs():
    if not os.path.exists(pdir):
        os.mkdir(pdir)

    if not os.path.exists(ddir):
        os.mkdir(ddir)

# descriptor for saving the plots to a default name under plots/ and the data under plots/data
def defaultsave(f):
    def saved(*args, **kw):
        ret = f(*args, **kw)
        # compatibiity
        if isinstance(ret , Sequence):
            outd, ax = ret 
        else:
            outd = ret

        result_dirs()
        #out = os.path.join(pdir, f.__name__+'.png')
        out = os.path.join(pdir, f.__name__+'.svg')
        if outd:
            datout = os.path.join(ddir, f.__name__+'.dat')
            with open(datout, 'w') as fo:
                io.dump(outd, fo)

        plt.savefig(out, format='svg')
        return ret
        #plt.show()
    return saved

def metasave(fmeta, metad):
    # Saves a meta data file of key=value pairs to fmeta
    kvline = '{}={}\n'
    with open(fmeta, 'w') as fo:
        for k,v in metad.items():
            fo.write(kvline.format(k, v))
def metaread(fmeta):
    # let (#) denote comments
    outd = OrderedDict()
    with open(fmeta, 'r') as fi:
        for line in fmeta:
            if line[0] == '#':
                continue
            k, v = line.split('=')
            outd[k] = v
    return outd

# convenience for dealing the .pkl files
def _justifynum(num):
    return '{:0>10}'.format(num)
def _makename(num, fn, ext='.pkl'):
    fn = '_'.join([fn, _justifynum(num)])  + ext
    return fn
def _loadpkl(fname):
    return pickle.load(open(fname, 'rb'))
def _outnum(name):
    name = path.basename(name)
    bname, _ = path.splitext(name)
    return int(bname.split('_')[-1])

# Compile data from all the directories into plots/data
# We needn't always recreate a new, compiled data file. The original data is sufficient.
# Need to couple plot data with metadata... 
# Metadata can come from conf file

# Ok fine lets make a class for this purpose.
# We want to gather simple column data with headers from several files 
# presenting them in an easy format to plot while preserving information about
# where they came from and giving ways to gather the associated metadata
# from .conf files or .meta files or .json configuration file
# pdir and ddir global variables in this file at the moment

# If necessary to use configuration object then try and use Silke read_param.Param class

# If we are in a simulation directory then I want to be able to just return one dataset
incself = False # ugly global variable for running on current directory
class Plotcompile(object):
    def __init__(self, dataname, dirglob='*/'):
        # Write now just stick to compiling data from files with the same name 
        self.dataname = dataname
        self.dirglob = dirglob
        self.pdata = [] 
        self.json = 'configuration.json'
        self._gather(dataname)
        # data structure list of dicts 
        # [{'data':<OrderedDict>{header:nparray}, 'dir':<dirname>, 'conf':<path_to_conf>,
        #       'meta':{}}]

    def _gather(self, dataname):
        if incself: 
            # only interested in current directory
            dirs=['./']
        else:
            # start by using io.finddirs to search through all the current directory sub-dirs
            dirs = io.finddirs(self.dirglob)

        for di in dirs:
            print di
            datapath = path.join(di, ddir)
            dataf = path.join(datapath, dataname+'.dat')
            cglob = glob( path.join(di, '*.conf') )
            assert len(cglob) == 1
            confpath = cglob[0]
            ind = io.freaddump(dataf)
            js = self.load_json(di)
            entry = {'data':ind, 'dir':di, 'conf':confpath, 'json':js}
            self.pdata.append(entry)
        return self.pdata
    
    def load_json(self, di):
        jsonpath = path.join(di, self.json)
        if not os.path.exists(jsonpath):
            return None
        with open(jsonpath, 'r') as fj:
            return json.load(fj)

    def get_type_conf(self, tp):
        # look for the for pair_type_param vp line for this type and read the parameters
        pass

@defaultsave
def tracer(dfile):
    with open(dfile, 'r') as fi:
        ddump = io.readdump(fi)
        x, y = ddump.values()
        x = x[20:800]
        y = y[20:800]
        lx= np.log10(np.asarray(x))
        ly= np.log10(np.asarray(y))
        #plt.xlim([10**3, 10**7])
        #plt.ylim([10**-2, 10**3])
        plt.xscale('log')
        plt.yscale('log')

        # compute slope
        m, c = np.polyfit(lx,ly,1)
        print m, c
        def fitline(lx):
            return 10**(m*lx + c)
        
        #fspace = np.linspace(lx[0], lx[-1], len(x), True)
        fy = map(fitline, lx)
        plt.loglog(x,y,'.b',x,fy,'--g')
        plt.title('gradient fitted line is {}'.format( m ))

        # specific lines
        plt.ylabel('MSD')
        plt.rc('text', usetex=True)
        plt.xlabel(r'$\tau$')


        out = 'plots/'+ path.splitext(path.split(dfile)[1])[0] + '.png'
        plt.savefig(out)

# Just read a dump of a line graph and plot it
def pl(dfile, loglog=False):
    with open(dfile, 'r') as fi:
        ddump = io.readdump(fi)
        x = ddump.values()[0]
        plt.xlabel(ddump.keys()[0])

        for head, data in ddump.items()[1:]:
            if loglog:
                plt.loglog(x, data, label=head)
            else:
                plt.plot(x, data, label=head)



        # compute the slope of log

        out = 'plots/'+ path.splitext(path.split(dfile)[1])[0] + '.png'
        print 'saving plot to ', out
        plt.savefig(out)

# specficially for plotnghosts. 
# dataname = 'plotnghosts'
# The first iteration of this kind of cross directory data collection
def nghosts_compile(dataname):
    # mkdir plots/data in this directory as well
    dataname = 'plots/data/' + dataname + '.dat'
    result_dirs()
    outd= OrderedDict()
    xlines = []
    xkeys = []
    datlines = []
    datakeys = []
    def compile_macro(dataname):
        def macro():
            with open(dataname,'r') as fi:
                ind = io.readdump(fi)
                xlines.append(ind.values()[0])
                xkeys.append(ind.keys()[0])
                datlines.extend(ind.values()[1:])
                datakeys.extend(ind.keys()[1:])
        return macro
        
    dirs = io.diriterate(compile_macro(dataname))

    # Now specific to comparing the ghost numbers
    # assert xkeys and xlines are all the same value
    outd[xkeys[0]] = xlines[0]
    def legend_rule(dd):
        dd = dd.strip('/')
        dd =dd.replace('_', '=')
        return dd[:-1] + '.' + dd[-1]
    heads = map(legend_rule, dirs)
    for head, data in zip(heads, datlines):
        outd[head] = data
    dataout, ext  = path.splitext(dataname)
    dataout = dataout + '_compiled' + ext
    with open(dataout, 'w') as fo:
        io.dump(outd, fo)

# just for plotting, no reading
@defaultsave
def nttplot(ntt, log=False):
    x = ntt.keys()
    y = ntt.values()
    #plt.plot(x, y)
    if log:
        plt.loglog(x, y)
    else:
        plt.plot(x, y)
    plt.xlabel('Shape index')
    plt.ylabel('No. of T1 transitions in 3000 timesteps with v0=0.2')
    plt.show()
    outd =  OrderedDict()
    outd['sindex'] = x
    outd['No._T1_transitions'] = y
    return outd

# want to collect data for every subdirectory
def nttcompare():
    ptname= 'nttplot'
    pldirs = sorted(glob('./*/plots/data/'))
    print 'reading data from', pldirs
    def reader(di):
        datf = os.path.join(di, ptname+'.dat')
        pld = {}
        with open(datf, 'r') as fo:
            ind = io.readdump(fo)
            pld['x'], pld['y']  = ind.values()
        return pld
    pldata = map(reader, pldirs)
    v0l= [0.01, 0.05, 0.10, 0.20]
    for v0, pld in zip(v0l, pldata):
        plt.plot(pld['x'], pld['y'], label='v0 {}'.format(v0))
    plt.show()

# plot for different directories the avg_stresses from stresses_0000001000.pkl
@defaultsave
def adj_compare():
    name = 'stresses_0000001000.pkl'
    adirs = sorted(glob('./adj_*'))
    def reader(di):
        datf = os.path.join(di, name)
        stresses = _loadpkl(datf)
        shear = np.mean(stresses['virial'].sstress.values())
        pressure = stresses['virial'].pressure
        return (pressure, shear)
    adata = map(reader, adirs)
    pressures, shears = zip(*adata)
    adjx = [di[-1] for di in adirs]
    try:
        adjx =map(int ,adjx)
    except ValueError:
        sys.exit('Could not find the values of adj from the directory names')
    plt.plot(adjx, shears)


# helper for plotting the number of boundary particles
@defaultsave
def plotnghosts(x, y):
    plt.xlabel('timestep')
    plt.ylabel('No. Boundary Particles')
    plt.plot(x, y)
    return {'timestep':x, 'N_Boundary_Particles':y}

# helper for quickly plotting the neighbour topology
@defaultsave
def nntopology(nnlist, outnum):
    mn = min(nnlist)
    mx = max(nnlist)
    bins = np.arange(mn-0.5, mx+0.5+0.01, 1.)
    plt.hist(nnlist, bins=bins)
    plt.xlabel('No. of sides')
    plt.title('Neighbour topology of cell mesh. timestep {}'.format(outnum))
    outd = OrderedDict()
    outd['id'], outd['n_number'] = zip(*list(enumerate(nnlist)))
    return outd

# helper for quickly plotting the msd curves for the types of particles
@defaultsave
def plotmsd(timeline, msd):
    plt.xlabel('timestep')
    plt.ylabel('MSD')
    outd= {'timestep':timeline}
    for tp, msdarr in msd.items():
        plt.loglog(timeline, msdarr, label='type {}'.format(tp))
        outd['msd_type_{}'.format(tp)] = msdarr
    #plt.plot(timeline, msd)
    plt.legend()
    return outd

# helper for plotting the number of cells, and hence growth of the tissue

@defaultsave
def growthexp(timeline=None, n_cells=None, ax=ax):
    if not timeline:
        # read timeline and n_cells instead
        gdata = io.freaddump('plots/data/growthexp.dat')
        timeline, n_cells = gdata['timestep'], gdata['N_cells']

    # specific
    fl = [100000/1000, 700000/1000]

    # general
    #ax.ticklabel_format(style='sci', axis='x')
    ax.xaxis.set_ticks([0, 500000, 1000000])
    #ax.xaxis.set_ticklabels([r'$0$',  r'$5\times 10^5$',  r'$10^6$'])
    ax.xaxis.set_ticklabels([r'$0$',  r'2500',  r'5000'])


    ax.set_xlabel(r'time ($\tau$)')
    ax.set_ylabel('Number of cells')
    #ax.set_title('Growth rate of patch')

    ax.semilogy(timeline, n_cells)
    deg = 1

    m, c = np.polyfit(timeline[fl[0]:fl[1]], map(np.log, n_cells[fl[0]:fl[1]]), deg)
    print 'fitting: gradient = {}, c = {}'.format(m, c)
    def fitline(x):
        return np.exp(c) * np.exp(m * x)
    fspace = np.linspace(timeline[0], timeline[-1], 10001)
    fy = map(fitline, fspace)
    ax.semilogy(fspace, fy, linestyle='--')

    outd = OrderedDict()
    outd['timestep']= timeline
    outd['N_cells']  = n_cells
    outd['lin_fit'] = fy

    return outd


# For plotting the radial distribution from a .pkl file
@defaultsave
def radial(num=1000000, stname='virial'):
    plt.xlabel('radial distance')
    plt.ylabel('radial normal stress')

    pklname = _makename(num, 'stresses')
    data = _loadpkl(pklname)

    stress = data[stname]
    plt.clf()

    xs = stress.rspace
    xdiff = xs[1] - xs[0]
    x = xs[:-1] + xdiff/2

    k = 'radial_normal_stress'
    lk = stname + '_' + k
    plt.plot(x, stress.ravg, label=lk, marker='o')

    #plt.plot(x, stress.ravg[k], label=lk, marker='o')

    #k = 'radial_shear_stress'
    #lk = stname + '_' + k
    #plt.plot(x, stress.ravg[k], label=lk, marker='o')
    #plt.legend()


@defaultsave
def radial_all(gdir='tmp/', ax=ax):

    mpl.rcParams.update({'legend.fontsize': 12})
    fig = plt.gcf()
    fig.tight_layout()

    ax.set_xlabel('Radial Distance')
    ax.set_ylabel('Average Pressure')
    ax.set_ylim([0, 3.0])
    #ax.set_title('Radial pressure profiles. Separated by 100,000 timesteps.')


    fstress = sorted(glob(os.path.join(gdir, '*.pkl')))
    print fstress
    for pklname in fstress:

        data = _loadpkl(pklname)

        stress = data['virial']

        xs = stress.rspace
        xdiff = xs[1] - xs[0]
        x = xs[:-1] + xdiff/2

        k = 'radial_normal_stress'

        fd = int(str(_outnum(pklname))[0])
        if fd != 1:
            lk = r'${}\times 10^5$'.format( fd )
        else:
            lk = r'$10^6$'
        ax.plot(x, stress.ravg[k], label=lk)
        ax.plot(x, stress.ravg[k], 'bo', label=lk, )

    #legend =ax.legend(title='timestep')
    #legend.get_title().set_fontsize('18.0')
    #legend.get_frame().set_linewidth(1)

# helper
def shiftrspace(rspace):
    xdiff = rspace[1] - rspace[0]
    x = rspace[:-1] + xdiff/2
    return x

@defaultsave
def arearadial(rspace, radiff, ax=ax):
    ax.set_xlabel('radial distance')
    ax.set_ylabel('$A - A_{max}$')
    x = shiftrspace(rspace)
    ax.plot(x, radiff)
    outd = {'r':x, 'adiff':radiff}
    return outd, ax

#matplotlib.rcParams['text.usetex'] = 'false'
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['xtick.major.size'] = 16
mpl.rcParams['xtick.minor.size'] = 0
mpl.rcParams['ytick.major.size'] = 16
mpl.rcParams['ytick.minor.size'] = 0
mpl.rcParams['font.size']=24.0
mpl.rcParams['legend.fontsize']=18.0

@defaultsave
def grow_radial():
    figsize = (12.0,6)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    #f.subplots_adjust(bottom=0.2)
    #fig.subplots_adjust(hspace=.5)

    growthexp(ax=ax1)
    #ax1.set_title('Sharing Y axis')
    radial_all('tmp/', ax=ax2)



## plotting texture tensor

@defaultsave
def dev_texture(xx_shear):
    plt.clf()
    x = xx_shear.keys()
    ymean = map(np.mean, xx_shear.values())
    yupper = map(max, xx_shear.values())
    ylower = map(min, xx_shear.values())
    plt.plot(x, ymean, marker='o')
    plt.plot(x, yupper, linestyle='--', color='g', marker='o')
    plt.plot(x, ylower, linestyle='--', color='g', marker='o')

    outd = OrderedDict()
    outd['adjn'] = x
    outd['ymean'] = ymean
    outd['yupper'] = yupper
    outd['ylower'] = ylower
    return outd


########################################################################################
# pretty plotting
# These methods operate on data produced by the analysis routines and are intended to be 
# adapted to produce production quality plots.

# collect any available msd data (name is plotmsd) and add to it details about 
# what the types of particles actually mean
# 

# global convenience
#last_title =

# store plot configurations in dictionaries, one for each of the simulations
#shape = {'msd_type_2':'lambda 0.58', 'msd_type_3':'lambda 0.7'}
#k = {msd_type_2:'K 1.0', msd_type_3:'K 20.0'}
#gamma = {msd_type_2:'gamma 0.1', msd_type_3:'gamma 0.2'}

#@defaultsave
#def pretty_msd(dirglob='*/', loglog=True, typeconf=shape):
    #dataname = 'plotmsd'

    ##returning one dataset
    #pdata = Plotcompile(dataname).pdata

    #typeconf = {'msd_type_2':'lambda 0.58', 'msd_type_3':'lambda 0.7'}
    #thisdata = pdata[0]['data']

    #x = thisdata['timestep']
    #y1 = thisdata['msd_type_2']
    #y2 = thisdata['msd_type_3']

    #fig, ax = plt.subplots()

    #if loglog:
        #pltf = ax.loglog
        #print 'using loglog plot'
    #else:
        #pltf = ax.plot
        #print 'using ax.plot function'
    #pltf(x, y1, label=typeconf['msd_type_2'])
    #pltf(x, y2, label=typeconf['msd_type_3'])

    #ax.set_ylabel('MSD')
    #ax.set_xlabel('timestep')
    #title = 'High p0 cells travel through low p0 tissue'
    #ax.set_title(title)

    #ax.legend(loc='lower right')
#pretty_msd.evfs = [str, eval]




# Determine shape of the neighbour topology histogram and plot as bar chart with spaces
@defaultsave
def pretty_nnt(compare='active'):
    #hack
    mpl.rcParams['xtick.major.size'] = 4
    mpl.rcParams['ytick.major.size'] = 4
    #plt.gcf().subplots_adjust(bottom=0.5)

    # search up the data
    dataname = 'nntopology'
    glo = 're_growA0'
    pdata = Plotcompile(dataname, glo ).pdata

    if compare == 'growth':
        def retrieve_growth(ddir):
            return float(ddir.strip('/').split('_')[1])/1000
        def get_legend(grs):
            keys = ['{}'.format(gr) for gr in grs]
            return keys
    elif compare == 'active':
        def retrieve_growth(ddir): return ddir.strip('/')
        def get_legend(grs): 
            legmap = {'divonly':'divison only', 'activediv':'division + activity'}
            return [legmap[gr] for gr in grs]

    nttall = []
    grs = []
    for dd in pdata:
        grs.append( retrieve_growth(dd['dir']) )
        ntts = np.array(dd['data'].values()[1], dtype=int)
        nttall.append(ntts)
    mn = min([ntt.min() for ntt in nttall])
    mx = max([ntt.max() for ntt in nttall])
    sides = np.array(range(mn, mx+1))
    nbars = len(sides)

    #want to normalise the results and give fractions for each side number
    nttnormed = []
    barmap = dict(zip(sides, range(nbars)))
    for ntt in nttall:
        norm = np.zeros(nbars)
        # bin the data
        for n in ntt:
            norm[barmap[n]] += 1
        N = ntt.size
        nttnormed.append(norm/N)

    # 
    fig, ax = plt.subplots(figsize=(10,10))
    fig.subplots_adjust(bottom=0.2)

    barspace = 0.7 # 0.3 free space
    step = barspace/len(nttnormed)
    colors = ['b', 'r', 'g', 'y']

    rcts = []
    #keys = get_legend(grs)
    for i, normed in enumerate(nttnormed):
        shift = i * step
        sides = sides[1:-3]
        normed = normed[1:-3]
        print normed
        rct = ax.bar(sides + shift, normed, step, color=colors[i%4])
        rcts.append(rct)

    ax.set_xticks(sides + barspace/2.)
    ax.set_xticklabels(sides)
    ax.set_ylabel('Fraction of cells')
    ax.set_xlabel('Number of neighbours')
    #ax.set_title(raw_input('title: '))
    
    #legend = ax.legend(rcts, keys, title=r'growth rate ($\eta$)')
    #legend.get_title().set_fontsize('18.0')
    #legend.get_frame().set_linewidth(1)


if __name__=='__main__':

    import sys
    import argparse

    # This is a cute command line function that takes the first argument to be the name
    #  of a function defined locally and any other arguments as arguments to the function
    args = sys.argv[1:]
    nargs = len(args)
    if nargs is 0:
        print 'cmplotting <function> [arguments]'
        calls = [name for name, loc in locals().items() if callable(loc)
                and name[0] != '_']
        for call in calls:
            print call
        sys.exit()

    # going to also process flags first
    # Not a general solution for this kind of command line tool hacking 
    # looks for flags that start with '-'
    flagcheck = [arg[0] == '-' for arg in args]
    farg= (i for i, b in enumerate(flagcheck) if b == False).next()
    commandargs = args[farg:]
    flags = args[:farg]
    for f in flags:
        if f == '-this':
            incself = True

    args = commandargs
    fname = args[0]
    print 'using function ', fname
    fargs = args[1:]
    try:
        f_using = locals()[fname]
    except KeyError:
        print 'I don\'t have a function \'%s\'', fname
        raise

    if len(fargs) is 0:
        print 'No arguments given to cmplotting.%s' % fname

        if hasattr(f_using, 'defaults'):
            fargs = f_using.defaults
            print 'Using defaults ', f_using.defaults
        else:
            fargs = []

    #print map(eval, fargs)
    if hasattr(f_using, 'evfs'):
        ff = [ev(i) for i, ev in zip(fargs, f_using.evfs)]
        print 'using evfs to convert arguments'
    else:
        ff = fargs
    print ff

    f_using(*ff)



