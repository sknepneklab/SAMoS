from constructmesh import *

def test_line():
    l1 = Line((0,0), (0,1))
    l2 = Line((1,1), (1,0))
    pti = l1.intersect(l2)
    print pti

def test_hexlattice():
    hx = hexlattice()
    with open('hex.tmp', 'w') as f:
        for pt in hx:
            print pt
            f.write(' '.join(map(str, pt))+ '\n')


if __name__=='__main__':
    #test_line()
    test_hexlattice()

