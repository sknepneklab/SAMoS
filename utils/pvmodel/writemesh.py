
import sys
from openmesh import *


import numpy as np
import vtk

from collections import OrderedDict

from ioutils import omvec, idtopt

# Generic parts

def vtk_write(polydata, outfile):
    writer = vtk.vtkXMLPolyDataWriter()
    if vtk.VTK_MAJOR_VERSION <= 5:
        writer.SetInput(polydata)
    else:
        writer.SetInputData(polydata)

    print 'saving ', outfile 

    writer.SetFileName(outfile)
    writer.SetDataModeToAscii()
    writer.Write()


# take an openmesh object and write it as .vtk with faces
def writemesh(mesh, outfile):


    Points = vtk.vtkPoints()

    Faces = vtk.vtkCellArray()

    for vh in mesh.vertices():
        pt =omvec(mesh.point(vh))
        Points.InsertNextPoint(pt)

    for fh in mesh.faces():
        vhids = []
        for vh in mesh.fv(fh):
            # need to store all the relevant vertex ids
            vhids.append(vh.idx())
        n = len(vhids)
        Polygon = vtk.vtkPolygon()
        Polygon.GetPointIds().SetNumberOfIds(n)
        for i, vhid in enumerate(vhids):
            Polygon.GetPointIds().SetId(i, vhid)
        Faces.InsertNextCell(Polygon)


    polydata = vtk.vtkPolyData()
    polydata.SetPoints(Points)
    polydata.SetPolys(Faces)

    polydata.Modified()
    writer = vtk.vtkXMLPolyDataWriter()
    if vtk.VTK_MAJOR_VERSION <= 5:
        writer.SetInput(polydata)
    else:
        writer.SetInputData(polydata)

    writer.SetFileName(outfile)
    writer.SetDataModeToAscii()
    writer.Write()



# take an openmesh object and write it as .vtk with faces
def writemeshenergy(pv, outfile):
    
    tri = pv.tri
    mesh = pv.mesh
    Points = vtk.vtkPoints()
    Faces = vtk.vtkCellArray()

    energy = vtk.vtkDoubleArray()
    energy.SetNumberOfComponents(1)
    energy.SetName("energy")
    enprop = pv.enprop

    boundary = vtk.vtkDoubleArray()
    boundary.SetNumberOfComponents(1)
    boundary.SetName("boundary")

    for vh in mesh.vertices():
        pt =omvec(mesh.point(vh))
        Points.InsertNextPoint(pt)

    # throw the edge vertices on the end
    nmv = mesh.n_vertices()
    boundary_vmap  = {}
    i = 0 # the boundary vertice count
    for bd in pv.boundaries:
        for vhid in bd:
            boundary_vmap[vhid] = nmv + i
            pt = omvec(tri.point(tri.vertex_handle(vhid)))
            Points.InsertNextPoint(pt)
            i += 1

    #print Points.GetNumberOfPoints()

    # Can actually add the boundary polygons here and show them 
    nfaces =0
    for vh in tri.vertices():
        vhids = []
        is_boundary = tri.is_boundary(vh)

        fh = mesh.face_handle(vh.idx())
        for meshvh in mesh.fv(fh):
            # need to store all the relevant vertex ids
            vhids.append(meshvh.idx())

        n = len(vhids)
        Polygon = vtk.vtkPolygon()
        Polygon.GetPointIds().SetNumberOfIds(n)
        for i, vhid in enumerate(vhids):
            Polygon.GetPointIds().SetId(i, vhid)
        Faces.InsertNextCell(Polygon)
        nfaces += 1

        fen = tri.property(enprop, vh)
        energy.InsertNextValue(fen)

        b_id = tri.property(pv.boundary_prop, vh)
        boundary.InsertNextValue(b_id)

    print 'added faces', nfaces



    polydata = vtk.vtkPolyData()
    polydata.SetPoints(Points)
    polydata.SetPolys(Faces)

    polydata.GetCellData().AddArray(energy)
    polydata.GetCellData().AddArray(boundary)

    polydata.Modified()
    writer = vtk.vtkXMLPolyDataWriter()
    if vtk.VTK_MAJOR_VERSION <= 5:
        writer.SetInput(polydata)
    else:
        writer.SetInputData(polydata)

    print 'saving ', outfile 

    writer.SetFileName(outfile)
    writer.SetDataModeToAscii()
    writer.Write()



# take an openmesh object and write it as .vtk with faces
def writetriforce(pv, outfile):

    tri = pv.tri

    Points = vtk.vtkPoints()
    Faces = vtk.vtkCellArray()

    force = vtk.vtkDoubleArray()
    force.SetNumberOfComponents(3)
    force.SetName("force")
    fprop = pv.fprop

    idtriv = vtk.vtkDoubleArray()
    idtriv.SetNumberOfComponents(1)
    idtriv.SetName("id")


    for vh in pv.tri.vertices():
        pt =omvec(pv.tri.point(vh))
        Points.InsertNextPoint(pt)
        fov = pv.tri.property(fprop, vh)

        force.InsertNextTuple3(*fov)
        idtriv.InsertNextValue(vh.idx())

    for fh in tri.faces():
        vhids = []
        for vh in tri.fv(fh):
            # need to store all the relevant vertex ids
            vhids.append(vh.idx())
        n = len(vhids)
        Polygon = vtk.vtkPolygon()
        Polygon.GetPointIds().SetNumberOfIds(n)
        for i, vhid in enumerate(vhids):
            Polygon.GetPointIds().SetId(i, vhid)
        Faces.InsertNextCell(Polygon)


    polydata = vtk.vtkPolyData()
    polydata.SetPoints(Points)
    polydata.SetPolys(Faces)

    polydata.GetPointData().AddArray(force)
    polydata.GetPointData().AddArray(idtriv)

    polydata.Modified()

    vtk_write(polydata, outfile)


def rotation_2d(theta):
    c, s = np.cos(theta), np.sin(theta)
    R = np.matrix('{} {}; {} {}'.format(c, -s, s, c))
    return R

def add_ellipse(Points, evals, shift, R, res):
    x, y = shift
    a, b = evals
    thetas = np.linspace(0, 2*np.pi, res, True)
    for i, th in enumerate(thetas):
        rot = np.einsum('mn,n->m',  R, np.array([a * cos(th), b * sin(th) ]))
        xe, xy = rot
        Points.InsertPoint(i, xe + x, xy + y, 0.)
        print 'inserting point', xe + x, xy + y, 0.
        
 
from numpy.linalg import eig
import ioutils as io
def write_stress_ellipses(pv, outfile, res=10):
    alpha = 1. # scaling factor

    Points = vtk.vtkPoints()

    e_x = np.array([1., 0., 0.])
    io.stddict(pv.stress)

    # calculate principle stresses
    evalues, evectors = eig(pv.stress.values())
    for vh in pv.tri.vertices():
        vhid = vh.idx()
        pt = io.idtopt(pv.tri, vhid)
        # cut down to two dimensions 
        shift =  pt[:2]
        print evalues
        evals =  evalues[vhid][:2]
        ea, eb, _ = evectors[vhid]
        theta = np.arccos(np.dot(e_x, ea))
        R = rotation_2d(theta)
        add_ellipse(Points, evals, shift, R, res)

        polyline = vtk.vtkPolyLine()
        polyline.GetPointIds().SetNumberOfIds(res)
        ptstart = vhid*res
        for j in range(ptstart, ptstart + res):
            polyline.GetPointIds().SetId(vhid, j)

        ells = vtk.vtkCellArray()
        ells.InsertNextCell(polyline)
        break

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(Points)
    
    polydata.SetLines(ells)
    polydata.Modified()

    vtk_write(polydata, 'test_ellipse.vtp')



from math import cos, sin
def write_test_ellipse(res=50):

    Points = vtk.vtkPoints()
    a = 1
    b = 2

    thetas = np.linspace(0, 2*np.pi, res, True)
    for i, th in enumerate(thetas):
        Points.InsertPoint(i, a*cos(th)+ 0,b*sin(th)+ 0, 0)
        
    aPolyLine = vtk.vtkPolyLine()
    aPolyLine.GetPointIds().SetNumberOfIds(res)
    for j in range(res):
        aPolyLine.GetPointIds().SetId(j, j)

    ells = vtk.vtkCellArray()
    ells.InsertNextCell(aPolyLine)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(Points)
    
    polydata.SetLines(ells)
    polydata.Modified()

    vtk_write(polydata, 'test_ellipse.vtk')

if __name__=='__main__':
    #write_test_ellipse()
    pass

