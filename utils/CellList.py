# * *************************************************************
# *  
# *   Soft Active Mater on Surfaces (SAMoS)
# *   
# *   Author: Rastko Sknepnek
# *  
# *   Division of Physics
# *   School of Engineering, Physics and Mathematics
# *   University of Dundee
# *   
# *   (c) 2013, 2014
# * 
# *   School of Science and Engineering
# *   School of Life Sciences 
# *   University of Dundee
# * 
# *   (c) 2015
# * 
# *   Author: Silke Henkes
# * 
# *   Department of Physics 
# *   Institute for Complex Systems and Mathematical Biology
# *   University of Aberdeen  
# * 
# *   (c) 2014, 2015
# *  
# *   This program cannot be used, copied, or modified without
# *   explicit written permission of the authors.
# * 
# * ***************************************************************


#
#  Cell class
#
#  Cell list class
#
#

from copy import deepcopy

# Single cell of the CellList object
# it has its own index (label), position r (==(rx,ry,rz)) and extension L (==(Lx,Ly,Lz))
class Cell:
  
  def __init__(self, idx, r, L):
    self.idx = idx
    self.r = r
    self.L = L
    self.indices = []
    # Labels (indices) of the neighbouring cells
    self.neighbors = []
    
  def add_vertex(self, idx):
    self.indices.append(idx)
    
  
class CellList:
  # Create my boxes
  def __init__(self,geom, r_cut, box):
    self.geom=geom
    self.r_cut = r_cut
    self.box = box
    self.cell_indices = {}
    # The size of the box collection is always set by the system box size
    # Which is what is used in the C++ code as well
    Lx, Ly, Lz = self.box
    self.Lx, self.Ly, self.Lz = Lx, Ly, Lz 
    # Number and linear extensions of the boxes in all three dimensions
    nx = int(Lx/r_cut)
    ny = int(Ly/r_cut)
    nz = int(Lz/r_cut)
    dx = Lx/float(nx)
    dy = Ly/float(ny)
    dz = Lz/float(nz)
    self.nx, self.ny, self.nz = nx, ny, nz
    self.dx, self.dy, self.dz = dx, dy, dz
    # total number of cells
    n_cell = nx*ny*nz
    print n_cell
    # Cell list is a python list
    self.cell_list = [None for i in range(n_cell)]
    for i in range(nx):
      x = -0.5*Lx + float(i)*dx
      for j in range(ny):
        y = -0.5*Ly + float(j)*dy
        for k in range(nz):
          z = -0.5*Lz + float(k)*dz
          # Cell labeling scheme: for each x, do all y, and for all y, do all z
          idx = ny*nz*i + nz*j + k
          # Create new cell with index, position and size
          self.cell_list[idx] = Cell(idx,(x,y,z),(dx,dy,dz))
          # Find neighbours: up to one above and below in all directions, moduly PBC
          # Note that a cell is its own neigbhbour ...
          for ix in range(-1,2):
            for iy in range(-1,2):
              for iz in range(-1,2):
                iix, iiy, iiz = i + ix, j + iy, k + iz
                if (iix == nx): iix = 0
                if (iiy == ny): iiy = 0
                if (iiz == nz): iiz = 0
                if self.geom.periodic:
                  if (iix < 0):  iix = nx - 1
                  if (iiy < 0):  iiy = ny - 1
                  if (iiz < 0):  iiz = nz - 1
                # A neighbour is just a label of the neighbouring cell
                self.cell_list[idx].neighbors.append(ny*nz*iix + nz*iiy + iiz)
  
  # Get the cell label for a given position vector v
  def get_cell_idx(self, v):
    x, y, z = v
    xmin, ymin, zmin = -0.5*self.Lx, -0.5*self.Ly, -0.5*self.Lz
    xmax, ymax, zmax = 0.5*self.Lx, 0.5*self.Ly, 0.5*self.Lz
    if self.geom.periodic:
      if (x >= xmax): x -= self.Lx
      elif (x < xmin): x += self.Lx
      if (y >= ymax): y -= self.Ly
      elif (y < ymin): y += self.Ly
      if (z >= zmax): z -= self.Lz
      elif (z < zmin): z += self.Lz
    i, j, k = int((x-xmin)/self.dx), int((y-ymin)/self.dy), int((z-zmin)/self.dz) 
    cell_idx = self.ny*self.nz*i + self.nz*j + k
    return cell_idx
  
  # Add a particle to a cell: This means compute its cell index (if not given already)
  # Then add it to ??? recursive call (I think we may have meant add_vertex of cell, not cell list?)
  # 
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
      