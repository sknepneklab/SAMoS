from build_epithilial import *

def stripmain():
    lx = 30
    ly = 100
    area = 3.0
    bp = 0.5
    N = int(round((lx * ly)/area))
    print 'number of cells in strip', N
    planea = Plane(lx,ly,N,0.0,bp)
    planeb = Plane(lx,ly,N,0.0,bp)
    planec = Plane(lx,ly,N,0.0,bp)
    outfile = 'cout.dat'
    # the strip separation
    sep = 45
    x = lx/2 + sep
    ptt = [(-x, 0), (0,0),  (x, 0)]

    pset = [planea, planeb, planec]
    combinewrite(outfile, pset, ptt=ptt)


if __name__=='__main__':
    stripmain()


