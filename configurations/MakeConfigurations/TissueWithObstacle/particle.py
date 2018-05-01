import numpy as np

class Particle:

    def __init__(self, id):
        self.id = id 
        self.r = np.array([0.0,0.0,0.0])
        self.radius = 1.0
        self.v = np.array([0.0,0.0, 0.0])
        self.type = 1
        self.n = np.array([1.0,0.0, 0.0])
        self.boundary = 0
        self.in_tissue = 0
        self.area = 0.0
    
    def write(self,out):
        x, y, z = self.r 
        vx, vy, vz = self.v 
        nx, ny, nz = self.n
        out.write('{:4d} {:2d} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.6f} {:1d}  {:1d}\n'.format(self.id,self.type,self.radius,x,y,z,nx,ny,nz,vx,vy,vz,0.0,0.0,1.0,self.area,self.boundary,self.in_tissue))

    def __add__(self,p):
        return self.r + p.r 
    
    def __sub__(self,p):
        return self.r - p.r

    def dist(self,r):
        dr = self - r 
        return np.sqrt(np.sum(dr**2))