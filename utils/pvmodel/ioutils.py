import numpy as np

### These are general methods copied from my command.py module

#print square data to file, first column is int and rest are floats.
def dump(dd, fo, htag='#'):
    nc = len(dd.keys())
    fo.write(''.join([htag+' ', '%s\t'*nc, '\n']) % tuple(dd.keys()))
    ddv = dd.values()
    nr = len(ddv[0]) # assumption
    outstr = '%d\t' + '%f\t'*(nc-1) + '\n' # assumption
    for i in range(nr):
        tup = tuple([ddvi[i] for ddvi in ddv])
        fo.write(outstr % tup)

def datdump(dd, fo):
    htag = 'keys:'
    dump(dd, fo, htag=htag)

def readdump(fo):
    headers = fo.next()[1:].split()
    dd = {}
    for h in headers:
        dd[h] = []
    for line in fo:
        for i, dat in enumerate(line.split()):
            ev = float
            if i is 0:
                ev = int
            dd[headers[i]].append( ev(dat) )
    return dd
    

# debugging

def dirk(A):
    print A
    print dir(A)
    sys.exit()
def shiv(al):
    for a in al:
        print a
        print eval(a)

#np.set_printoptions(threshold=np.nan)

import contextlib
import cStringIO
@contextlib.contextmanager
def nostdout():
    save_stdout = sys.stdout
    sys.stdout = cStringIO.StringIO()
    yield
    sys.stdout = save_stdout

# want to print a vector object
def dumpvec(vec):
    print omvec(vec)

def scatter(mesh, mesh2):
    arrl = []
    for vh in mesh.vertices():
        arrl.append(omvec(mesh.point(vh)))
    npts = np.column_stack(arrl)
    x = npts[:][0]
    y = npts[:][1]
    plt.scatter(x, y, color='red')

    mesh = mesh2
    arrl = []
    for vh in mesh.vertices():
        arrl.append(omvec(mesh.point(vh)))
    npts = np.column_stack(arrl)
    x = npts[:][0]
    y = npts[:][1]
    plt.scatter(x, y, color='blue')

    plt.show()


# how to print a dictionary containing serious data
def stddict(dd):
    for k, v in dd.items():
        print k
        print v
        print 
    



# Important!
# this is shared code for cellmesh and writemesh for interacting with openmesh

# openmesh has a vector object
# too lazy to use this to convert to numpy arrays
# fixed to three dimensions...

def omvec(vec):
    return np.array([vec[0], vec[1], vec[2]])

def idtopt(mesh, rmuid):
    return omvec(mesh.point(mesh.vertex_handle(rmuid)))


