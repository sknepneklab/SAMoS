
import sys
from openmesh import *

import numpy as np
import vtk

from collections import OrderedDict

# copied from cellmesh.py
def omvec(vec):
    return np.array([vec[0], vec[1], vec[2]])

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
def writemeshenergy(mesh, outfile):


    Points = vtk.vtkPoints()
    Faces = vtk.vtkCellArray()

    energy = vtk.vtkDoubleArray()
    energy.SetNumberOfComponents(1)
    energy.SetName("energy")
    enprop = FPropHandle()
    mesh.get_property_handle(enprop, 'energy')

    force = vtk.vtkDoubleArray()
    force.SetNumberOfComponents(3)
    force.SetName("force")
    fprop = FPropHandle()
    assert mesh.get_property_handle(fprop, 'force')

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

        fen = mesh.property(enprop, fh)
        energy.InsertNextValue(fen)

        fov = mesh.property(fprop, fh)
        force.InsertNextTuple3(*fov)


    polydata = vtk.vtkPolyData()
    polydata.SetPoints(Points)
    polydata.SetPolys(Faces)

    polydata.GetCellData().AddArray(energy)
    polydata.GetCellData().AddArray(force)

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
def writetriforce(pv, outfile):

    mesh = pv.tri

    Points = vtk.vtkPoints()
    Faces = vtk.vtkCellArray()


    force = vtk.vtkDoubleArray()
    force.SetNumberOfComponents(3)
    force.SetName("force")
    fprop = FPropHandle()
    assert pv.mesh.get_property_handle(fprop, 'force')


    imforce = vtk.vtkDoubleArray()
    imforce.SetNumberOfComponents(3)
    imforce.SetName("imforce")

    nnforce = vtk.vtkDoubleArray()
    nnforce.SetNumberOfComponents(3)
    nnforce.SetName("nnforce")

    idtriv = vtk.vtkDoubleArray()
    idtriv.SetNumberOfComponents(1)
    idtriv.SetName("id")

    for vh in pv.tri.vertices():
        pt =omvec(pv.tri.point(vh))
        Points.InsertNextPoint(pt)
        # hack to deal with boundaries
        fov = [0., 0., 0.]
        fim = [0., 0., 0.]
        fnn = [0., 0., 0.]
        if not pv.tri.is_boundary(vh):
            # get appropriate face
            #fh = pv.mesh.face_handle(vh.idx())
            fh = pv.mesh.face_handle(pv.vfmap[vh.idx()])
            fov = pv.mesh.property(fprop, fh)
            fim = pv.mesh.property(pv.imfprop, fh)
            fnn = pv.mesh.property(pv.nnfprop, fh)
            #print fh.idx()
            #print fov
        force.InsertNextTuple3(*fov)
        imforce.InsertNextTuple3(*fim)
        nnforce.InsertNextTuple3(*fnn)
        idtriv.InsertNextValue(vh.idx())

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

    polydata.GetPointData().AddArray(force)
    polydata.GetPointData().AddArray(imforce)
    polydata.GetPointData().AddArray(nnforce)
    polydata.GetPointData().AddArray(idtriv)

    polydata.Modified()
    writer = vtk.vtkXMLPolyDataWriter()
    if vtk.VTK_MAJOR_VERSION <= 5:
        writer.SetInput(polydata)
    else:
        writer.SetInputData(polydata)

    writer.SetFileName(outfile)
    writer.SetDataModeToAscii()
    writer.Write()


### These are general methods copied from my command.py module

#print square data to file, first column is int and rest are floats.
def dump(dd, fo, htag='#'):
    nc = len(dd.keys())
    fo.write(''.join([htag+' ', '%s\t'*nc, '\n']) % tuple(dd.keys()))
    ddv = dd.values()
    nr = len(ddv[0]) # assumption
    outstr = '%d\t' + '%f\t'*(nc-1) + '\n' # assumption
    for i in range(nr):
        tup = tuple([ddvi[i] for ddvi in ddv])
        fo.write(outstr % tup)

def datdump(dd, fo):
    htag = 'keys:'
    dump(dd, fo, htag=htag)

def readdump(fo):
    headers = fo.next()[1:].split()
    dd = {}
    for h in headers:
        dd[h] = []
    for line in fo:
        for i, dat in enumerate(line.split()):
            ev = float
            if i is 0:
                ev = int
            dd[headers[i]].append( ev(dat) )
    return dd
    


