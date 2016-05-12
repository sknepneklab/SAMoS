import numpy as np
from cellmesh import *

def test_bond_intersection():
    x_c = np.array([20.,20.])
    a = np.array([25., 20.])
    b = np.array([16.,20.])
    wl = 2.
    inter = bond_intersection(x_c, wl, b-a, a)
    print inter

if __name__=='__main__':
    test_bond_intersection()


