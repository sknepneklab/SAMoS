from particle import *

class Obstacle:

    def __init__(self, Rc, dens):
        self.pos = np.array([])
        if Rc.size != 3:
            raise Exception('Obstacle center has to be NumPy array with three elements.')
        self.Rc = Rc
        if dens < 0:
            raise Exception('Particle density has to be positive.')
        self.dens = dens
        self.obstacle_type = 3
        self.in_tissue = 0
        self.boundary = 0
        self.pos = []

    def make(self):
        raise NotImplementedError()