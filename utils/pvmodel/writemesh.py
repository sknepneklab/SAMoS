
import sys
from openmesh import *


import numpy as np
import vtk

from collections import OrderedDict

from numpy.linalg import eig, norm
import ioutils as io
from ioutils import omvec, idtopt


from math import cos, sin
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

    cradius = vtk.vtkDoubleArray()
    cradius.SetNumberOfComponents(1)
    cradius.SetName("cradius")

    for i, vh in enumerate(mesh.vertices()):
        pt =omvec(mesh.point(vh))
        Points.InsertNextPoint(pt)
        is_tri = pv.mesh.property(pv.tri_prop, vh)
        crnv = 0. if is_tri else pv.cradius[i] 
        cradius.InsertNextValue(crnv)

    ## throw the edge vertices on the end
    #nmv = mesh.n_vertices()
    #boundary_vmap  = {}
    #i = 0 # the boundary vertice count
    #for bd in pv.boundaries:
        #for vhid in bd:
            #boundary_vmap[vhid] = nmv + i
            #pt = omvec(tri.point(tri.vertex_handle(vhid)))
            #Points.InsertNextPoint(pt)
            #cradius.InsertNextValue(0.)
            #i += 1

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
    polydata.GetPointData().AddArray(cradius)

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

def get_internal_stress(pv, pvstress):
    stress = OrderedDict()
    for vh in pv.tri.vertices():
        boundary = pv.tri.is_boundary(vh)
        if not boundary:
            stress[vh.idx()] = pvstress[vh.idx()]
    return stress

# take an openmesh object and write it as .vtk with faces
def writetriforce(pv, outfile, stress=True):

    tri = pv.tri

    Points = vtk.vtkPoints()
    Faces = vtk.vtkCellArray()

    idtriv = vtk.vtkDoubleArray()
    idtriv.SetNumberOfComponents(1)
    idtriv.SetName("id")

    force = vtk.vtkDoubleArray()
    force.SetNumberOfComponents(3)
    force.SetName("force")
    fprop = pv.fprop

    stress_1 = vtk.vtkDoubleArray()
    stress_1.SetNumberOfComponents(3)
    stress_1.SetName("stress_1")

    stress_2 = vtk.vtkDoubleArray()
    stress_2.SetNumberOfComponents(3)
    stress_2.SetName("stress_2")

    # Stress calculation fails for boundaries
    # Stress at boundaries is not symmetric
    if stress:
        stress = get_internal_stress(pv, pv.n_stress)
        evalues, evectors = eig(stress.values())
    #io.stddict(pv.stress)
    #print evalues
    #print 
    #print evectors
    
    for vh in pv.tri.vertices():
        vhid = vh.idx()
        pt =omvec(pv.tri.point(vh))
        Points.InsertNextPoint(pt)
        fov = pv.tri.property(fprop, vh)

        idtriv.InsertNextValue(vh.idx())
        force.InsertNextTuple3(*fov)
        if stress:
            if not pv.tri.is_boundary(vh):
                s_1 = evalues[vhid][0] * evectors[vhid][0]
                s_2 = evalues[vhid][1] * evectors[vhid][1]
                stress_1.InsertNextTuple3(*s_1)
                stress_2.InsertNextTuple3(*s_2)
            else:
                stress_1.InsertNextTuple3(*np.zeros(3))
                stress_2.InsertNextTuple3(*np.zeros(3))

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

    polydata.GetPointData().AddArray(idtriv)
    polydata.GetPointData().AddArray(force)
    if stress:
        polydata.GetPointData().AddArray(stress_1)
        polydata.GetPointData().AddArray(stress_2)

    polydata.Modified()
    vtk_write(polydata, outfile)


def rotation_2d(theta):
    c, s = np.cos(theta), np.sin(theta)
    R = np.matrix('{} {}; {} {}'.format(c, -s, s, c))
    return R

def add_ellipse(Points, evals, shift, R, res):
    x, y = shift
    a, b = evals
    a, b = abs(a), abs(b)
    thetas = np.linspace(0, 2*np.pi, res, True)
    for i, th in enumerate(thetas):
        xy_i = np.array([a * cos(th), b * sin(th) ])
        xy_i = np.einsum('mn,n->m',  R, xy_i)
        xe, xy = xy_i
        Points.InsertNextPoint(xe + x, xy + y, 0.)
        #print 'inserting point', xe + x, xy + y, 0.
 
def write_stress_ellipses(pv, outfile, pvstress, res=20):
    alpha = 1. # scaling factor

    Points = vtk.vtkPoints()

    e_x = np.array([1., 0., 0.])

    # calculate principle stresses
    stress = get_internal_stress(pv, pvstress)
    evalues, evectors = eig(stress.values())

    #maxe1 =max(map(abs, list(evalues[:][0])))
    #maxe2 =max(map(abs, list(evalues[:][1])))
    #normf = maxe1
    normf= 1.

    ells = vtk.vtkCellArray()
    for vhid, _  in enumerate(evalues):
        vh = pv.tri.vertex_handle(vhid)

        vhid = vh.idx()
        pt = io.idtopt(pv.tri, vhid)
        # cut down to two dimensions 
        shift =  pt[:2]
        evals =  evalues[vhid][:2]
        av, bv, = evals
        a, b = av/normf, bv/normf
        ea, eb, _ = evectors[vhid]
        theta = np.arccos(np.dot(e_x, ea)) 
        sg = 1. if np.dot(np.cross(e_x, ea),pv.normal) > 0 else -1.
        R = rotation_2d(sg * theta)  
        add_ellipse(Points, (a,b), shift, R, res)

        polyline = vtk.vtkPolyLine()
        polyline.GetPointIds().SetNumberOfIds(res)
        ptstart = vhid*res
        for j, ptj in enumerate(range(ptstart, ptstart + res)):
            polyline.GetPointIds().SetId(j, ptj)

        ells.InsertNextCell(polyline)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(Points)
    
    polydata.SetLines(ells)
    polydata.Modified()

    vtk_write(polydata, outfile)

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

