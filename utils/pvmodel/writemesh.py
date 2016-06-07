
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
    
    mesh = pv.mesh
    Points = vtk.vtkPoints()
    Faces = vtk.vtkCellArray()

    #energy = vtk.vtkDoubleArray()
    #energy.SetNumberOfComponents(1)
    #energy.SetName("energy")
    #enprop = pv.enprop

    #cradius = vtk.vtkDoubleArray()
    #cradius.SetNumberOfComponents(1)
    #cradius.SetName("cradius")

    #pressure = vtk.vtkDoubleArray()
    #pressure.SetNumberOfComponents(1)
    #pressure.SetName("pressure")

    meshpt = mesh.pym
    for i, vh in enumerate(mesh.vertices()):
        pt = meshpt[vh.idx()]
        Points.InsertNextPoint(pt)
    
    # Can actually add the boundary polygons here and show them 
    for mf in mesh.faces():
        mvhs = list(mesh.fv(mf))
        n = len(mvhs)
        Polygon = vtk.vtkPolygon()
        Polygon.GetPointIds().SetNumberOfIds(n)
        for i, mvh in enumerate(mvhs):
            Polygon.GetPointIds().SetId(i, mvh.idx())
        Faces.InsertNextCell(Polygon)
        
        #fen = tri.property(enprop, vh)
        #energy.InsertNextValue(fen)

        #pr = pv.stresses['virial'].pressure[vh.idx()]
        #prfacevalue = 0. if np.isnan(pr) else pr
        #pressure.InsertNextValue(prfacevalue)


    polydata = vtk.vtkPolyData()
    polydata.SetPoints(Points)
    polydata.SetPolys(Faces)

    #polydata.GetCellData().AddArray(energy)
    #polydata.GetCellData().AddArray(pressure)
    #polydata.GetPointData().AddArray(cradius)

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

# depreciated
#def get_internal_stress(pv, pvstress):
    #stress = OrderedDict()
    #for vh in pv.tri.vertices():
        #boundary = pv.tri.is_boundary(vh)
        #if not boundary:
            #stress[vh.idx()] = pvstress[vh.idx()]
    #return stress

def get_stress(pv, pvstress):
    stress = OrderedDict()
    for i, pvs in pvstress.items():
        if pvs is not None:
            stress[i]= pvs
    return stress
def get_stress(a,b): return b

# take an openmesh object and write it as .vtk with faces
def writetriforce(pv, outfile):

    tri = pv.tri
    tript = tri.pym

    Points = vtk.vtkPoints()
    Faces = vtk.vtkCellArray()

    idtriv = vtk.vtkDoubleArray()
    idtriv.SetNumberOfComponents(1)
    idtriv.SetName("id")

    force = vtk.vtkDoubleArray()
    force.SetNumberOfComponents(3)
    force.SetName("force")

    stress_1 = vtk.vtkDoubleArray()
    stress_1.SetNumberOfComponents(3)
    stress_1.SetName("stress_1")

    stress_2 = vtk.vtkDoubleArray()
    stress_2.SetNumberOfComponents(3)
    stress_2.SetName("stress_2")

    # Stress calculation fails for boundaries
    # Stress at boundaries is not symmetric
    #if pv.forces:
        #stress = get_stress(pv, pv.n_stress)
        #evalues, evectors = eig(stress.values())
    #io.stddict(pv.stress)
    #print evalues
    #print 
    #print evectors
    
    for vh in pv.tri.vertices():
        vhid = vh.idx()
        pt = tript[vhid]
        Points.InsertNextPoint(pt)
        idtriv.InsertNextValue(vh.idx())
        
        fov = pv.tri.forces[vhid]
        force.InsertNextTuple3(*fov)
        if False:
            if not pv.tri.is_boundary(vh):
                s_1 = evalues[vhid][0] * evectors[vhid][0]
                s_2 = evalues[vhid][1] * evectors[vhid][1]
                stress_1.InsertNextTuple3(*s_1)
                stress_2.InsertNextTuple3(*s_2)
            else:
                stress_1.InsertNextTuple3(*np.zeros(3))
                stress_2.InsertNextTuple3(*np.zeros(3))

    for fh in tri.faces():
        vhids = [vh.idx() for vh in tri.fv(fh)]
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
    if False:
        polydata.GetPointData().AddArray(stress_1)
        polydata.GetPointData().AddArray(stress_2)

    polydata.Modified()
    vtk_write(polydata, outfile)


def rotation_2d(theta):
    c, s = np.cos(theta), np.sin(theta)
    R = np.matrix('{} {}; {} {}'.format(c, -s, s, c))
    return R

def add_ellipse(Points, evals, shift, R, res, pressure, pr=0.):
    #print 'add ellipse ', evals
    x, y = shift
    a, b = evals
    a, b = abs(a), abs(b)
    thetas = np.linspace(0, 2*np.pi, res, True)
    for i, th in enumerate(thetas):
        xy_i = np.array([a * cos(th), b * sin(th) ])
        xy_i = np.einsum('mn,n->m',  R, xy_i)
        xe, xy = xy_i
        Points.InsertNextPoint(xe + x, xy + y, 0.)
        pressure.InsertNextValue(pr)
 
def write_stress_ellipses(pv, outfile, pvstress, res=20):
    # any visualisation has to do something with None stresses
    alpha = 1. # scaling factor

    Points = vtk.vtkPoints()

    tri = pv.tri
    mesh = pv.mesh
    tript = tri.pym
    meshpt = mesh.pym
    e_x = np.array([1., 0., 0.])

    # calculate principle stresses
    stress = pvstress.stress
    clist = pvstress.clist
    evalues, evectors = {}, {}
    for i in clist:
        ss = stress[i]
        if ss[0] is not np.nan:
            evalues[i], evectors[i] = eig(ss)

    pressure = vtk.vtkDoubleArray()
    pressure.SetNumberOfComponents(1)
    pressure.SetName("pressure")

    #maxe1 =max(map(abs, list(evalues.values())))
    #maxe2 =max(map(abs, list(evalues[:][1])))

    maxe = np.max(np.absolute(np.array(evalues.values())))
    if maxe == 0.:
        print evalues
    #normf = maxe/0.5
    normf = 1.
    print 'adjusting stress ellipses by a factor of ', normf
    
    ells = vtk.vtkCellArray()
    for ellid, vhid in enumerate(evalues.keys()):
        pt = tript[vhid]
        shift =  pt[:2]

        # cut down to two dimensions 
        evals =  evalues[vhid][:2]
        av, bv, = evals
        a, b = av/normf, bv/normf
        ea, eb, _ = evectors[vhid]
        theta = np.arccos(np.dot(e_x, ea)) 
        sg = 1. if np.dot(np.cross(e_x, ea),mesh.normal) > 0 else -1.
        R = rotation_2d(sg * theta)  


        pr = pvstress.pressure[vhid]
        prfacevalue = 0. if np.isnan(pr) else pr
        pressure.InsertNextValue(prfacevalue)
        add_ellipse(Points, (a,b), shift, R, res,  pressure, pr=prfacevalue )

        polyline = vtk.vtkPolyLine()
        polyline.GetPointIds().SetNumberOfIds(res)
        ptstart = ellid*res
        for j, ptj in enumerate(range(ptstart, ptstart + res)):
            polyline.GetPointIds().SetId(j, ptj)

        ells.InsertNextCell(polyline)


    polydata = vtk.vtkPolyData()
    polydata.SetPoints(Points)
    
    polydata.GetPointData().AddArray(pressure)
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

