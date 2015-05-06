#
#  Cell class
#
#  Cell list class
#
#

from copy import deepcopy

class Cell:
  
  def __init__(self, idx, r, L):
    self.idx = idx
    self.r = r
    self.L = L
    self.indices = []
    self.neighbors = []
    
  def add_vertex(self, idx):
    self.indices.append(idx)
    
  
class CellList:
  
  def __init__(self, r_cut, box):
    self.r_cut = r_cut
    self.box = box
    self.cell_indices = {}
    Lx, Ly, Lz = self.box
    self.Lx, self.Ly, self.Lz = Lx, Ly, Lz 
    nx = int(Lx/r_cut)
    ny = int(Ly/r_cut)
    nz = int(Lz/r_cut)
    dx = Lx/float(nx)
    dy = Ly/float(ny)
    dz = Lz/float(nz)
    self.nx, self.ny, self.nz = nx, ny, nz
    self.dx, self.dy, self.dz = dx, dy, dz
    n_cell = nx*ny*nz;
    self.cell_list = [None for i in range(n_cell)]
    for i in range(nx):
      x = -0.5*Lx + float(i)*dx
      for j in range(ny):
        y = -0.5*Ly + float(j)*dy
        for k in range(nz):
          z = -0.5*Lz + float(k)*dz
          idx = ny*nz*i + nz*j + k
          self.cell_list[idx] = Cell(idx,(x,y,z),(dx,dy,dz))
          for ix in range(-1,2):
            for iy in range(-1,2):
              for iz in range(-1,2):
                iix, iiy, iiz = i + ix, j + iy, k + iz
                if (iix < 0):  iix = nx - 1
                if (iix == nx): iix = 0
                if (iiy < 0):  iiy = ny - 1
                if (iiy == ny): iiy = 0
                if (iiz < 0):  iiz = nz - 1
                if (iiz == nz): iiz = 0
                self.cell_list[idx].neighbors.append(ny*nz*iix + nz*iiy + iiz)
                
  def get_cell_idx(self, v):
    xmin, ymin, zmin = -0.5*self.Lx, -0.5*self.Ly, -0.5*self.Lz
    xmax, ymax, zmax = 0.5*self.Lx, 0.5*self.Ly, 0.5*self.Lz
    x, y, z = v
    if (x >= xmax): x -= self.Lx
    elif (x < xmin): x += self.Lx
    if (y >= ymax): y -= self.Ly
    elif (y < ymin): y += self.Ly
    if (z >= zmax): z -= self.Lz
    elif (z < zmin): z += self.Lz
    i, j, k = int((x-xmin)/self.dx), int((y-ymin)/self.dy), int((z-zmin)/self.dz) 
    cell_idx = self.ny*self.nz*i + self.nz*j + k
    return cell_idx
    
  def add_vertex(self, v, idx, cell_index = None):
    if cell_index == None:
      cell_idx = self.get_cell_idx(v)
    else:
      cell_idx = cell_index
    self.cell_list[cell_idx].add_vertex(idx)
    self.cell_indices[idx] = cell_idx
    
 
  def wipe(self):
    for cell in self.cell_list:
      cell.vertices = []
      
  def get_neighbours(self,v):
      cell_index = self.get_cell_idx(v)
      neighbors = []
      for idx in self.cell_list[cell_index].neighbors:
          neighbors.extend(deepcopy(self.cell_list[idx].indices))
      return neighbors
      