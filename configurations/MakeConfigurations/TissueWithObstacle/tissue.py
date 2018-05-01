from particle import *
from polygon_test import *

class Tissue:

    def __init__(self, outer_boundary, inner_boundary, density):
        self.outer_boundary = outer_boundary
        self.inner_boundary = inner_boundary
        self.dens = density   # bulk density
        self.perim = perim(outer_boundary) + perim(inner_boundary)
        self.area = area(outer_boundary) - area(inner_boundary)
        self.stem_type = 2
        self.tissue_type = 1
        self.boundary_type = 4
        self.outer_boundary_normals = normals(outer_boundary)
        self.inner_boundary_normals = -normals(inner_boundary)
        self.pos = []
        self.boundary_tuples = []
        self.cell_area = np.pi

    def make_boundary(self, boundary_density, offset = 0):
        self.Nboundary = 0
        # make outer_boundary 
        Nseg = self.outer_boundary.shape[0] 
        outer_pos = []
        for i in range(Nseg):
            j = (i+1) % Nseg 
            p1 = self.outer_boundary[i,:]
            p = Particle(self.Nboundary + offset)
            p.radius = 1.0
            p.type = self.boundary_type
            p.r = np.array([p1[0],p1[1],0])
            p.v = np.array([0,0,0])
            p.n = self.outer_boundary_normals[i]
            p.boundary = 1
            p.in_tissue = 1
            self.pos.append(p)
            outer_pos.append(p.id)
            p2 = self.outer_boundary[j,:]
            l = np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)
            N = int(np.round(l*boundary_density)) 
            if N > 2:
                t = np.linspace(1/float(N),1,N-2,endpoint=False)
                x1, y1 = p1 
                x2, y2 = p2
                self.Nboundary += 1
                k = i + 1
                for tt in t:
                    x = (x2 - x1)*tt + x1 
                    y = (y2 - y1)*tt + y1 
                    p = Particle(self.Nboundary + offset)
                    p.radius = 1.0
                    p.type = self.boundary_type
                    p.r = np.array([x,y,0])
                    p.v = np.array([0,0,0])
                    p.n = np.array([(y2-y1)/l,-(x2-x1)/l,0])
                    p.boundary = 1
                    p.in_tissue = 1
                    self.pos.append(p)
                    outer_pos.append(p.id)
                    k += 1
                    self.Nboundary += 1
        # make inner boundary 
        Nseg = self.inner_boundary.shape[0] 
        inner_pos = []
        for i in range(Nseg):
            p1 = self.inner_boundary[i,:]
            p = Particle(self.Nboundary + offset)
            p.radius = 1.0
            p.type = self.boundary_type
            p.r = np.array([p1[0],p1[1],0])
            p.v = np.array([0,0,0])
            p.n = self.inner_boundary_normals[i]
            p.boundary = 1
            p.in_tissue = 1
            self.pos.append(p)
            inner_pos.append(p.id)
            self.Nboundary += 1
        for i in range(len(outer_pos)-1):
            self.boundary_tuples.append([i,outer_pos[i],outer_pos[i+1]])
        self.boundary_tuples.append([i,outer_pos[-1],outer_pos[0]])
        off = len(self.boundary_tuples)
        for i in range(len(inner_pos)-1):
            self.boundary_tuples.append([i+off,inner_pos[i],inner_pos[i+1]])
        self.boundary_tuples.append([len(self.boundary_tuples),inner_pos[-1],inner_pos[0]])
        

    def make_stem(self, p1, p2, stem_density = 1.0, offset = 0):
        x1, y1 = p1
        x2, y2 = p2 
        self.Nstem = 0
        l = np.sqrt((x2-x1)**2 + (y2-y1)**2)
        N = stem_density*l 
        t = np.linspace(0,1,N) 
        i = 0
        for tt in t:
            x = (x2 - x1)*tt + x1 
            y = (y2 - y1)*tt + y1 
            if not (inside_polygon(x, y, self.outer_boundary) and not inside_polygon(x, y, self.inner_boundary)):
                raise Exception('Both end points for the stem cell line have to be inside the tissue')
            p = Particle(i + offset)
            p.radius = 1.0
            p.type = self.stem_type
            p.r = np.array([x,y,0])
            p.v = np.array([0,0,0])
            p.n = np.array([(y2-y1)/l,-(x2-x1)/l,0])
            p.outer_boundary = 0
            p.in_tissue = 1
            p.area = np.random.normal(self.cell_area, 0.5)
            self.pos.append(p)
            i += 1
            self.Nstem += 1

    def make_bulk(self, offset = 0, max_attempt_factor = 10):
        N = int(np.round(self.dens*self.area) - self.Nstem)
        if N <= 0:
            raise Exception('Density is too low to add any bulk cells.')
        min_dist = np.sqrt(1.0/float(self.dens))
        xmin, xmax = np.min(self.outer_boundary[:,0]), np.max(self.outer_boundary[:,0])
        ymin, ymax = np.min(self.outer_boundary[:,1]), np.max(self.outer_boundary[:,1])
        i = 0
        attempt = 0
        while i < N:
            if attempt > max_attempt_factor*N:
                print('Warning! Could not add all particles.')
                break
            can_add = True 
            x = np.random.uniform(xmin,xmax)
            y = np.random.uniform(ymin,ymax)
            if inside_polygon(x,y,self.outer_boundary) and not inside_polygon(x,y,self.inner_boundary):
                p = Particle(i + offset)
                p.r = np.array([x,y,0])
                for pp in self.pos:
                    dr = p.dist(pp)
                    if dr < min_dist: 
                        can_add = False
                        break
                if can_add:
                    p.radius = 1.0
                    p.type = self.tissue_type
                    p.r = np.array([x,y,0])
                    p.v = np.array([0,0,0])
                    phi = np.random.uniform(0,2*np.pi)
                    p.n = np.array([np.cos(phi), np.sin(phi), 0])
                    p.outer_boundary = 0
                    p.in_tissue = 1
                    p.area = np.random.normal(self.cell_area, 0.5)
                    self.pos.append(p)
                    i += 1
            attempt += 1            
