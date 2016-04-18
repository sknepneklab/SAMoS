
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
        if not is_boundary:
            fh = mesh.face_handle(vh.idx())
            for meshvh in mesh.fv(fh):
                # need to store all the relevant vertex ids
                vhids.append(meshvh.idx())
        else:
            # We append the vertex id that we are assigning to the edge vertex above
            vhids = list(pv.halfcells[vh.idx()])
            vhids.append( boundary_vmap[vh.idx()] )
            #print vhids
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

    mesh = pv.tri

    Points = vtk.vtkPoints()
    Faces = vtk.vtkCellArray()

    force = vtk.vtkDoubleArray()
    force.SetNumberOfComponents(3)
    force.SetName("force")
    fprop = pv.fprop

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
        #fim = [0., 0., 0.]
        #fnn = [0., 0., 0.]
        if not pv.tri.is_boundary(vh):
            # get appropriate face
            #fh = pv.mesh.face_handle(vh.idx())
            fh = pv.mesh.face_handle(vh.idx())
            fov = pv.tri.property(fprop, vh)
            #fim = pv.mesh.property(pv.imfprop, fh)
            #fnn = pv.mesh.property(pv.nnfprop, fh)
        # tmp hack
        if fov is None:
            fov = np.zeros(3)

        force.InsertNextTuple3(*fov)

        #imforce.InsertNextTuple3(*fim)
        #nnforce.InsertNextTuple3(*fnn)
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
    #polydata.GetPointData().AddArray(imforce)
    #polydata.GetPointData().AddArray(nnforce)
    polydata.GetPointData().AddArray(idtriv)

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
    


