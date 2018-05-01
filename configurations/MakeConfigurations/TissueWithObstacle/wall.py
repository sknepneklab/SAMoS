from particle import *
import numpy as np

class Wall:

    def __init__(self,L,W,density,wall_type = 3):
        if L <= 0:
            raise Exception('Wall length has to be positive.')
        self.L = L
        if W <= 0:
            raise Exception('Wall width has to be positive.')
        self.W = W 
        if density <= 0:
            raise Exception('Wall particle density has to be positive.')
        self.den = density 
        self.wall_type = wall_type
        self.pos = []

    def make(self, offset = 0, mask = [1,1,1,1]):
        self.pos = []
        length = self.L*(mask[0]+mask[2]) + self.W*(mask[1]+mask[3])
        int_offset = 0
        # make top
        if mask[0] == 1:
            y = 0.5*self.W 
            N = int(np.round(self.den*self.L)) 
            xmin = -0.5*self.L
            for i in range(N):
                p = Particle(i + offset + int_offset)
                x = xmin + i*self.L/float(N)
                p.radius = 1.0 
                p.type = self.wall_type
                p.r = np.array([x,y,0])
                p.v = np.array([0,0,0])
                p.n = np.array([0,0,0])
                p.boundary = 0
                p.in_tissue = 0
                self.pos.append(p)
            int_offset += N
        # make right
        if mask[1] == 1:
            x = 0.5*self.L 
            N = int(np.round(self.den*self.W)) 
            ymax = 0.5*self.W 
            for i in range(N):
                p = Particle(i + offset + int_offset)
                y = ymax - i*self.W/float(N)
                p.radius = 1.0 
                p.type = self.wall_type
                p.r = np.array([x,y,0])
                p.v = np.array([0,0,0])
                p.n = np.array([0,0,0])
                p.boundary = 0
                p.in_tissue = 0
                self.pos.append(p)
            int_offset += N
        # make bottom
        if mask[2] == 1:
            y = -0.5*self.W 
            N = int(np.round(self.den*self.L)) 
            xmax = 0.5*self.L
            for i in range(N):
                p = Particle(i + offset + int_offset)
                x = xmax - i*self.L/float(N)
                p.radius = 1.0 
                p.type = self.wall_type
                p.r = np.array([x,y,0])
                p.v = np.array([0,0,0])
                p.n = np.array([0,0,0])
                p.boundary = 0
                p.in_tissue = 0
                self.pos.append(p)
            int_offset += N
        # make left 
        if mask[3] == 1:
            x = -0.5*self.L 
            N = int(np.round(self.den*self.W)) 
            ymin = -0.5*self.W 
            for i in range(N):
                p = Particle(i + offset + int_offset)
                y = ymin + i*self.W/float(N)
                p.radius = 1.0 
                p.type = self.wall_type
                p.r = np.array([x,y,0])
                p.v = np.array([0,0,0])
                p.n = np.array([0,0,0])
                p.boundary = 0
                p.in_tissue = 0
                self.pos.append(p)
            