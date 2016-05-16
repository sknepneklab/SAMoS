import numpy as np
from cellmesh import *

from collections import OrderedDict

def test_bond_intersection():
    case = OrderedDict()
    x_c = np.array([20.,20.])
    a = np.array([25., 20.])
    b = np.array([16., 20.])
    case['flat'] = np.array([[20.,20.], [25.,20.], [16., 20.]])
    case['inner'] = np.array([[0.,0.], [-1.,-1.], [1.,1.]])
    case['half'] = np.array([[0.,0.], [-1/2., 0.], [5.,1.]])
    wl = 2.
    for ck, cv in case.items():
        print 
        print 'case ', ck
        x_c, a, b = cv

        inter = bond_intersection(x_c, wl, a, b)
        if inter == None:
            return None
        else:
            ma, mb, line = inter
            print ma, mb
            print line(ma), line(mb)

if __name__=='__main__':
    test_bond_intersection()

