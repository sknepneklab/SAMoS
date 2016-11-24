
#### OLD plotting functions


# we are using numpys .npz now
def plot_radial(*fils):
    fig = plt.figure()
    axi = fig.add_subplot(111)
    #axi.show()
    axi = plt
    for fi in fils:
        data = np.load(fi)
        xs = data['rspace']
        assert len(xs) > 1
        xdiff = xs[1] - xs[0]
        x = xs[:-1] + xdiff/2
        rad = data['radial_pressure']
        axi.plot(x, rad, marker='o')
        #axi.draw()
        axi.show()
        time.sleep(0.05)

def plot_avg_pressures(fglob='stress_st*.npz'):
    stfiles = sorted(glob(fglob))
    plt.xlabel('step')
    plt.ylabel('average_pressure')
    stdatas = [np.load(st) for st in stfiles]
    xsteps = map(_outnum, stfiles)
    spressure = [np.mean(data['simple_radial_pressure']) for data in stdatas]
    plt.plot(xsteps, spressure, label='simple pressure', marker='o')
    vpressure = [np.mean(data['virial_radial_pressure']) for data in stdatas]
    plt.plot(xsteps, vpressure, label='virial pressure', marker='o')
    #hpressure = [np.mean(data['hardy_radial_pressure']) for data in stdatas]
    #plt.plot(xsteps, hpressure, label='hardy pressure', marker='o')
    plt.legend()
    plt.show()


from ioutils import _nanmean

@defaultsave
def avg_pressures(fglob='stresses*.pkl'):
    sfiles = sorted(glob(fglob))
    datas = map(_loadpkl, sfiles)
    plt.xlabel('step')
    plt.ylabel('average_pressure')
    xsteps = map(_outnum, sfiles)
    stn = datas[0]
    if 'simple' in stn:
        spressure= [_nanmean(stress['simple'].pressure) for stress in datas]
        plt.plot(xsteps, spressure, label='simple pressure', marker='o')
    if 'virial' in stn:
        vpressure= [_nanmean(stress['virial'].pressure) for stress in datas]
        plt.plot(xsteps, vpressure, label='virial pressure', marker='o')
    if 'hardy' in stn:
        hpressure= [_nanmean(stress['hardy'].pressure) for stress in datas]
        plt.plot(xsteps, np.array(hpressure), label='hardy pressure', marker='o')
    if 'hardy_vertices' in stn:
        hpressure= [_nanmean(stress['hardy_vertices'].pressure) for stress in datas]
        plt.plot(xsteps, np.array(hpressure), label='hardy vertex pressure', marker='o')
    if 'vertices' in stn:
        hpressure= [_nanmean(stress['vertices'].pressure) for stress in datas]
        plt.plot(xsteps, np.array(hpressure), label='vertex pressure', marker='o')
    # go ahead and plot the average force magnitude as well
    avgfs = [st['avg_force'] for st in datas]
    plt.plot(xsteps, avgfs, label='average force magnitude', marker='o')
    if 'avg_vertex_force' in stn:
        avgfs = [st['avg_vertex_force'] for st in datas]
        plt.plot(xsteps, avgfs, label='average vertex force', marker='o')
    plt.legend()


@defaultsave
def avg_stresses(stname='virial', fglob='stresses*.pkl'):
    sfiles = sorted(glob(fglob))
    datas = map(_loadpkl, sfiles)
    plt.xlabel('step')
    xsteps = map(_outnum, sfiles)
    if not stname in datas[0] and stname != 'all':
        print 'didn\'t find that stress name'
        return 
    if stname == 'all':
        names = datas[0].keys()
    else:
        plt.title(stname)
        names= [stname]
    for name in names:
        vpressure= [_nanmean(stress[name].sstress) for stress in datas]
        plt.plot(xsteps, vpressure, label=name+' shear stress', marker='o')
        hpressure= [_nanmean(stress[name].nstress) for stress in datas]
        plt.plot(xsteps, np.array(hpressure), label=name+' normal stress', marker='o', linestyle='-')
    plt.legend()


#cmplotting.py compare_radial 4000 stress_st . ../../rA_2.5/pressure_test/ ../../vpotential_only/pressure_test/
@defaultsave
def compare_radial(num, fname, *fdirs):
    plt.clf()
    fname = '_'.join([fname, _justifynum(num)])  + '.npz'
    paths = [path.join(fd, fname) for fd in fdirs]
    plt.xlabel('radial distance')
    plt.ylabel('pressure')
    for fi in paths:
        data = np.load(fi)
        xs = data['rspace']
        assert len(xs) > 1
        xdiff = xs[1] - xs[0]
        x = xs[:-1] + xdiff/2
        rad = data['radial_pressure']
        plt.plot(x, rad, marker='o', label=fi)

#cmplotting.py compare_radial 4000 stress_st . ../../rA_2.5/pressure_test/ ../../vpotential_only/pressure_test/
def compare_radial(num, fname, *fdirs):
    plt.clf()
    fname = '_'.join([fname, _justifynum(num)])  + '.npz'
    paths = [path.join(fd, fname) for fd in fdirs]
    plt.xlabel('radial distance')
    plt.ylabel('pressure')
    for fi in paths:
        data = np.load(fi)
        xs = data['rspace']
        assert len(xs) > 1
        xdiff = xs[1] - xs[0]
        x = xs[:-1] + xdiff/2
        rad = data['radial_pressure']
        plt.plot(x, rad, marker='o', label=fi)
    #http://matplotlib.org/users/legend_guide.html
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)
    plt.show()


def plot_wls(yname='radial_pressure', fglob='stress_wl*.npz'):
    #fname = makename(num, fname)
    npzl = sorted(glob(fglob))
    plt.xlabel('radial distance')
    plt.ylabel('radial pressure')
    print npzl
    for npz in npzl:
        data = np.load(npz)
        xs = data['rspace']
        xdiff = xs[1] - xs[0]
        x = xs[:-1] + xdiff/2
        y= data[yname]
        plt.plot(x, y, label=npz, marker='o')

    plt.legend()
    plt.show()


