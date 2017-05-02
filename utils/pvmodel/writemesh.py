
import sys
from openmesh import *


import numpy as np
import vtk

from collections import OrderedDict

from numpy.linalg import eig, norm
import ioutils as io
from ioutils import omvec, idtopt
import transitions as tr


from math import cos, sin
# Generic parts

verb = False

def vtk_write(polydata, outfile):
    writer = vtk.vtkXMLPolyDataWriter()
    if vtk.VTK_MAJOR_VERSION <= 5:
        writer.SetInput(polydata)
    else:
        writer.SetInputData(polydata)

    if verb: print 'saving ', outfile 

    writer.SetFileName(outfile)
    writer.SetDataModeToAscii()
    writer.Write()

def writepoints(pts, outfile):
    Points = vtk.vtkPoints()
    for pt in pts:
        Points.InsertNextPoint(pt)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(Points)
    polydata.Modified()
    vtk_write(polydata, outfile)

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

    #print 'writing polygons'
    vtk_write(polydata, outfile)

# take an openmesh object and write it as .vtk with faces
# the dual
def writemeshenergy(pv, outfile):
    
    mesh = pv.mesh
    Points = vtk.vtkPoints()
    Faces = vtk.vtkCellArray()

    #energy = vtk.vtkDoubleArray()
    #energy.SetNumberOfComponents(1)
    #energy.SetName("energy")
    #enprop = pv.enprop

    tp = vtk.vtkDoubleArray()
    tp.SetNumberOfComponents(1)
    tp.SetName("type")

    force = vtk.vtkDoubleArray()
    force.SetNumberOfComponents(3)
    force.SetName("force")

    area = vtk.vtkDoubleArray()
    area.SetNumberOfComponents(1)
    area.SetName("area")

    pressure = vtk.vtkDoubleArray()
    pressure.SetNumberOfComponents(1)
    pressure.SetName("pressure")

    sshear = vtk.vtkDoubleArray()
    sshear.SetNumberOfComponents(1)
    sshear.SetName("shear stress")

    # need to add an abitrary number of texture tensors
    # naming convention t_trace_<n> where n is the averaging number
    texs = []
    trace_arrs = []
    stexs = []
    shear_arrs = [] # for deviatoric part
    is_xx_trace = True if hasattr(pv, 'xx_trace') else False
    if is_xx_trace:
        for adjn, trace in pv.xx_trace.items():
            t_name = 't_trace_{}'.format(adjn)
            texs.append(t_name)
            tmarr = vtk.vtkDoubleArray()
            tmarr.SetNumberOfComponents(1)
            tmarr.SetName(t_name)
            trace_arrs.append(tmarr)

        for adjn, trace in pv.xx_shear.items():
            t_name = 't_shear_{}'.format(adjn)
            stexs.append(t_name)
            tsarr = vtk.vtkDoubleArray()
            tsarr.SetNumberOfComponents(1)
            tsarr.SetName(t_name)
            shear_arrs.append(tsarr)


    meshpt = mesh.pym
    vforces =mesh.vertex_force
    for i, vh in enumerate(mesh.vertices()):
        pt = meshpt[vh.idx()]
        Points.InsertNextPoint(pt)
        force.InsertNextTuple3(*vforces[vh.idx()])
    
    stress = pv.stresses['virial']
    press = stress.pressure; sstress = stress.sstress
    ar= pv.mesh.areas
    types = pv.tri.types
    for vhid in pv.tri.bulk:
        mfid = pv.tri.to_mesh_face[vhid]
        mf = pv.mesh.face_handle(mfid)

        mvhs = list(mesh.fv(mf))
        n = len(mvhs)
        Polygon = vtk.vtkPolygon()
        Polygon.GetPointIds().SetNumberOfIds(n)
        for i, mvh in enumerate(mvhs):
            Polygon.GetPointIds().SetId(i, mvh.idx())
        Faces.InsertNextCell(Polygon)

        tp.InsertNextValue(types[vhid])
        area.InsertNextValue(ar[mfid])
        pressure.InsertNextValue(press[vhid])
        sshear.InsertNextValue(sstress[vhid])

    if is_xx_trace:
        for i, tarr in enumerate(trace_arrs):
            adjn = pv.xx_trace.keys()[i]
            for traceval in pv.xx_trace[adjn]:
                tarr.InsertNextValue(traceval)

        for i, tarr in enumerate(shear_arrs):
            adjn = pv.xx_shear.keys()[i]
            for shearval in pv.xx_shear[adjn]:
                tarr.InsertNextValue(shearval)
        

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(Points)
    polydata.SetPolys(Faces)

    #polydata.GetCellData().AddArray(energy)
    polydata.GetCellData().AddArray(tp)
    polydata.GetCellData().AddArray(area)
    polydata.GetCellData().AddArray(pressure)
    polydata.GetCellData().AddArray(sshear)
    polydata.GetPointData().AddArray(force)
    if is_xx_trace:
        for tarr in trace_arrs:
            polydata.GetCellData().AddArray(tarr)
        for tarr in shear_arrs:
            polydata.GetCellData().AddArray(tarr)

    polydata.Modified()
    writer = vtk.vtkXMLPolyDataWriter()
    if vtk.VTK_MAJOR_VERSION <= 5:
        writer.SetInput(polydata)
    else:
        writer.SetInputData(polydata)

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
    force.SetName("Vertex Force")

    aforce = vtk.vtkDoubleArray()
    aforce.SetNumberOfComponents(3)
    aforce.SetName("Active Force")

    stress_1 = vtk.vtkDoubleArray()
    stress_1.SetNumberOfComponents(3)
    stress_1.SetName("stress_1")

    stress_2 = vtk.vtkDoubleArray()
    stress_2.SetNumberOfComponents(3)
    stress_2.SetName("stress_2")

    # add arrows for the the boundary pulling direction
    extpull = vtk.vtkDoubleArray()
    extpull.SetNumberOfComponents(3)
    extpull.SetName("bd force direction")
    bdnormals = tri.get_boundary_normal(tri.rcm)


    for vh in pv.tri.vertices():
        vhid = vh.idx()
        pt = tript[vhid]
        Points.InsertNextPoint(pt)
        idtriv.InsertNextValue(vh.idx())
        
        fov = pv.tri.forces[vhid]
        force.InsertNextTuple3(*fov)
        factive = pv.tri.factive[vhid]
        aforce.InsertNextTuple3(*factive)
        if False:
            if not pv.tri.is_boundary(vh):
                s_1 = evalues[vhid][0] * evectors[vhid][0]
                s_2 = evalues[vhid][1] * evectors[vhid][1]
                stress_1.InsertNextTuple3(*s_1)
                stress_2.InsertNextTuple3(*s_2)
            else:
                stress_1.InsertNextTuple3(*np.zeros(3))
                stress_2.InsertNextTuple3(*np.zeros(3))
        if vhid in bdnormals:
            extpull.InsertNextTuple3(*bdnormals[vhid])
        else:
            extpull.InsertNextTuple3(*np.zeros(3))
            

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
    polydata.GetPointData().AddArray(extpull)
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
 
# Need a 'write ellipses' function which takes matrix dictionary
def write_stress_ellipses(pv, outfile, pvstress, res=20, usecentres=True, scale=1., 
        normalise=False):
    
    Points = vtk.vtkPoints()

    tri = pv.tri
    mesh = pv.mesh
    tript = tri.pym
    meshpt = mesh.pym
    e_x = np.array([1., 0., 0.])
    e_y = np.array([0., 1., 0.])

    # calculate principle stresses
    stress = pvstress
    evalues, evectors = {}, {}
    for i in stress.keys():
        ss = stress[i]
        evalues[i], evectors[i] = eig(ss)

    pressure = vtk.vtkDoubleArray()
    pressure.SetNumberOfComponents(1)
    pressure.SetName("pressure")

    if normalise:
        evs = [v for row in list(evalues.values()) for v in row]
        maxe1 =max(map(abs, evs))
        scale = maxe1

    maxe = np.max(np.absolute(np.array(evalues.values())))
    if maxe == 0.:
        print evalues
    normf = pv.tri.kproperty[0]
    if normf == 0.:
        normf = 1.
    normf *= scale

    print 'adjusting stress ellipses by a factor of ', normf
    
    ells = vtk.vtkCellArray()
    for ellid, vhid in enumerate(stress.keys()):
        if usecentres:
            pt = tript[vhid]
        else:
            pt = mesh.pym[vhid]
        shift =  pt[:2]

        # cut down to two dimensions 
        evalues, evectors = eig(stress[vhid])

        av, bv, _ = evalues
        a, b = av/normf, bv/normf
        ea, eb = evectors[:,0], evectors[:,1]

        # arrange so the larger eigenvalue corresponds to the longer ellipse axis
        swap = abs(b) >= abs(a)
        if swap:
            a,b=b,a
            ea,eb=eb,ea
        pvec = a * ea 

        theta = np.arccos(np.dot(e_x, pvec/norm(pvec))) 
        sg = 1. if np.cross(e_x, pvec)[2] > 0 else -1.
        R = rotation_2d(sg * theta)  

        pr = 1/2. * np.trace(stress[vhid])
        prfacevalue = pr if vhid in stress else 0.
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

