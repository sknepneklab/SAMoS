
from openmesh import *
import writemesh as wr
import ioutils as io
from ioutils import omvec, OrderedSet

import numpy as np
from numpy.linalg import norm, eig
import scipy.integrate as integrate


import matplotlib.pyplot as plt
from collections import OrderedDict
from scipy.spatial import Delaunay
from tmp_read_data import ReadData

from math import pi 
import math as m
import sys


# Samos
from CellList2D import CellList2D


# Yes, these are defined here! I use them a lot
npI = np.identity(3)
def rmmultiply(v, M):
    return np.einsum('n,nm->m', v, M)

# debugging 
def diagnose(mesh):
    print 'nedges', mesh.n_edges()
    print 'nvertices', mesh.n_vertices()
    print 'nfaces', mesh.n_faces()

from matplotlib import pyplot as plt
def scat(xyz):
    plt.scatter(xyz[:,0],xyz[:,1])
    plt.show()

### IGNORE. Utility methods for overly complicated hardy stress calculation  

# A function for integration which takes and arraylike and the space that it is defined on
#  and initial and final values 
def hashintegrate(arromega, arrspace, mm, mp):
    arra = arrspace[0]; arrz = arrspace[-1]
    arrange = arrz-arra
    ne = len(arrspace)
    # mm and mp are values in arrange
    dx = (arrz - arra)/(ne-1)
    imm = int(round(mm/dx))
    imp = int(round(mp/dx)) 
    # cut the array and then integrate with simpsons method
    #print imm, imp
    print mm, mp
    return integrate.simps(arromega[imm:imp], x=arrspace[imm:imp])

def in0_1(x):
    return x >= 0 and x <= 1

def bond_intersection(x_c, wl, x_a, x_b):
    # x_c circle centre
    # wl circle radius
    # bond is the vector from particle x_a to particle x_b
    # x_a is the position of particle a, 
    # x_a is going to be the second vertex in the bond pairs

    bond = x_b-x_a
    nca = norm(x_c - x_a)
    ncb = norm(x_c - x_b)

    def line(m):
        return m*bond + x_a

    if nca < wl and ncb < wl:
        return 0, 1, line
    
    # using formula for quadratic ax^2 + bx + c = 0
    a = np.dot(bond, bond)
    b = 2 * np.dot(bond, x_a-x_c)
    c = np.dot(x_c,x_c) + np.dot(x_a,x_a) -2*np.dot(x_c,x_a) - wl**2
    disc = b**2 - 4 * a * c
    if disc <= 0:
        return None 
    else:
        outa, outb = None, None
        sdisc = np.sqrt(disc)
        m_plus = (-b + sdisc)/(2*a)
        m_minus = (-b - sdisc)/(2*a)
        # worry about the exactly which of these should be >,< and >=,<=
        if nca < wl and ncb >= wl:
            # m_plus and m_minus should be one positive one negative in these two cases
            assert m_plus * m_minus  <= 0
            m = m_plus if m_plus > 0 else m_minus
            y= line(m)
            outa, outb = 0, m
            #print 'found a line partly in the smoothing region'
        elif ncb < wl and nca >= wl:
            # the smaller of the two absolute values 
            assert m_plus >= 0.
            assert m_minus >= 0.
            m_plus = abs(m_plus); m_minus = abs(m_minus)
            m = m_plus if m_plus < m_minus else m_minus
            y = line(m)
            outa, outb = m, 1
            #print 'found a line partly in the smoothing region'
        elif nca >= wl and ncb >= wl:
            outa, outb= m_minus, m_plus
            # The case where a bond line cuts the averaging zone but 
            # the bond itself does not
            if not in0_1(m_minus) or not in0_1(m_plus):
                return None

        if outa is None:
            sys.exit('failed to determine bond intersection')
        return outa, outb, line

# NOT USED
def henonarea(areatri):
    s = sum(areatri)/2
    a,b,c = areatri
    area = np.sqrt(s*(s-a)*(s-b)*(s-c))
    return area

# Want an object which wraps the Openmesh objects
# We want to easily get list of boundary and bulk elements
# We want to easily map from vertices of one mesh to the faces of another
# It's hard to say whether Openmesh is offering advantages over 
# a pure python mesh implementation --- 

# e_z vector GLOBALLY defined
normal = np.array([0.,0., 1.])
class Pymesh(object):
    # This object just holds methods which are commmon the triangle mesh and polymesh sub classes 

    # return the pair of vertices associated with an edge in any order
    def getvpair(self, eid):
        eh =self.edge_handle(eid)
        assert eh.idx() >= 0
        heh =self.halfedge_handle(eh, 0)
        aux = frozenset([self.to_vertex_handle(heh).idx(), self.from_vertex_handle(heh).idx()])
        return aux

    # get the vertex ids associated with a halfedge
    def vids_from_to(self, heh):
        to = self.to_vertex_handle(heh).idx()
        fr = self.from_vertex_handle(heh).idx()
        return (to, fr)

    def is_boundary_edge(self, heh):
        tovh = self.to_vertex_handle(heh)
        frvh = self.from_vertex_handle(heh)
        return (self.is_boundary(tovh) and self.is_boundary(frvh))

    # Create a mapping between pairs of vertices and edges
    def vedgemap(self):
        tri = self
        vedge = {}
        for eh in tri.edges():
            heh = tri.halfedge_handle(eh, 0)
            vt = tri.to_vertex_handle(heh).idx()
            vf = tri.from_vertex_handle(heh).idx()
            vedge[(vt, vf)] = heh.idx()
            heho = tri.halfedge_handle(eh, 1)
            vedge[(vf, vt)] = heho.idx()
        self.vedge = vedge

    # This is called during the construction of mesh objects to create some useful numpy arrays and data structures
    def finalize(self):
        self._pythonise()
        self._helengths()
        self.normal = normal

    # Openmesh stores vertex positions in such an awkward way, here we duplicate them in to numpy array.
    def _pythonise(mesh):
        meshpt = OrderedDict()
        for vh in mesh.vertices():
            meshpt[vh.idx()] = omvec(mesh.point(vh))
        mesh.pym = meshpt
        
    # Trading memory for speed by storing edge lengths and vectors that we can re-use
    # lvec - maps halfedge ids onto a vector
    # ll - maps halfedge ids onto halfedge lengths
    # lle - maps edge ids onto edge lengths
    def _helengths(self):
        tri = self
        pym = tri.pym
        lvec = {}
        ll = {}
        lle = {}
        for eh in tri.edges():
            heh = tri.halfedge_handle(eh, 0)
            heho = tri.opposite_halfedge_handle(heh)
            vt = tri.to_vertex_handle(heh)
            vf = tri.from_vertex_handle(heh)
            npvt = pym[vt.idx()]
            npvf = pym[vf.idx()]
            vedge = npvt - npvf
            lvec[heh.idx()] = vedge
            lvec[heho.idx()] = -vedge
            nve = norm(vedge)
            ll[heh.idx()] = nve
            ll[heho.idx()] = nve
            lle[eh.idx()] = nve
        self.lvec = lvec 
        self.ll = ll
        self.lle = lle

    # Storing squared lengths of halfedges
    def store_llsquared(self):
        ll = self.ll
        llsquared = dict(zip(ll.keys(), map(lambda x: x*x, ll.values())))
        self.llsquared = llsquared
        
    # create a cell list in order to group the vertices in space
    def clorganise(self, wl):
        Lx, Ly, mzero = np.amax(np.abs(self.pym.values()),axis=0)
        # get the maximum cell spacing
        cut = max(tri.ll.values()) + 0.1
        if wl < cut:
            print 'Warning. smoothing length is less than the maximum interparticle spacing'
        cl = CellList2D([2*Lx+1, 2*Ly+1], cut)
        for i, ipt in self.pym.items():
            cl.add_particle(ipt[:2], i)
        self.cl = cl
        return cl

    # -- some tools for iterating through parts of the openmesh --
    # going to assume that we constuct pymesh objects from Samos output and don't mutate the mesh
    # Return value is a python iterable which yeilds halfedge ids.
    def iterable_boundary(self, vh):
        # vh is the starting vertex which should be on the boundary
        mesh = self
        vhidstart = vh.idx()
        # A quick check for a boundary halfedge using the fact that it doesn't have a face 
        def is_boundary_he(mesh, he):
            fh = mesh.face_handle(he)
            return fh.idx() is -1
        start_he = None
        # We need to find the boundary halfedge specifically
        for heh in mesh.voh(vh):
            if is_boundary_he(mesh, heh):
                #print 'found boundary edge for vhid', vh.idx()
                start_he = heh
                break
        assert start_he is not None
        def iterable(start_he, mesh):
            he_i = start_he
            while True:
                yield he_i
                he_i = mesh.next_halfedge_handle(he_i)
                if he_i.idx() == start_he.idx():
                    raise StopIteration
        self._iterable_boundary = iterable(start_he, mesh)
        return self._iterable_boundary

    # -- Finding the boundary and bulk 
    @property
    def boundary(self):
        if hasattr(self, '_boundary'):
            return self._boundary
        mesh = self
        for vh in mesh.vertices():
            if mesh.is_boundary(vh):
                break
        itb = self.iterable_boundary(vh)
        self._boundary = [mesh.to_vertex_handle(heh).idx() for heh in list(itb)]
        return self._boundary
    
    # support for multiple boundaries is patchy at best!
    # Returns a list of boundaries, each is a list of halfedge ids
    def get_boundaries(self):
        # iterate through all the vertices and use sets to find all the boundaries
        self.boundaries = []
        allboundary = []
        for vh in self.vertices():
            if self.is_boundary(vh):
                allboundary.append(vh.idx())
        allboundary = set(allboundary)
        while allboundary: # is not empty
            vhid = allboundary.pop()

            itb = self.iterable_boundary(self.vertex_handle(vhid))
            boundary = [self.from_vertex_handle(heh).idx() for heh in list(itb)]
            self.boundaries.append(boundary)

            allboundary = allboundary - set(boundary) 

        return self.boundaries

    # get the unique pairs of connected vertices on the boundary of the triangulation
    def get_boundary_connectivity(self):
        self.get_boundaries()
        bconnect = [] # [(i, j)]
        for boundary in self.boundaries:
            for i, j in zip(boundary, np.roll(boundary, -1)):
                bconnect.append((i, j))
        return bconnect


    # map between boundary vertices and their outgoing boundary halfedges?
    @property
    def boundary_hegdes(self):
        if hasattr(self, '_boundary_hedges'):
            return self._boundary_hedges
        tboundary = self.boundary
        vhi = tboundary[0]
        vh = self.vertex_handle(vhi)
        itb = self.iterable_boundary(vh)
        heids = [heh.idx() for heh in list(itb)]
        self._boundary_hedges = OrderedDict(zip(tboundary, heids))
        return self._boundary_hedges

    # just returns the face ids of all the faces adjacent to a boundary vertex
    def iterate_boundary_vertex(self, hehp):
        tri = self
        assert tri.is_boundary(hehp)
        facels = []
        while True:
            heho = tri.opposite_halfedge_handle(hehp)
            if tri.is_boundary(heho):
                break
            facels.append(tri.face_handle(heho).idx())
            hehp = tri.prev_halfedge_handle(heho)
        return facels
 
    # Simply iterate the vertex faces but make sure it is done anticlockwise
    # Can also do tri.vf(vh) and then reverse it 
    def iterate_vertex(tri, vh):
        start = tri.halfedge_handle(vh)
        facels = []
        he = start
        while True:
            fh = tri.face_handle(he)
            facels.append(fh.idx())
            hehp = tri.next_halfedge_handle(he)
            heho = tri.opposite_halfedge_handle(hehp)
            he = heho
            if he.idx() == start.idx():
                break
        return facels

    # outgoing halfedges from a vertex
    def get_star_hedges(self, vhi):
        polygons= self
        hehids = []
        for fh in polygons.vf(vhi):
            for heh in polygons.fh(fh):
                hehids.append(heh.idx())
        return hehids

    # outgoing edges from a vertex
    def get_star_edges(self, vhi):
        polygons = self
        ehids = []
        for fh in polygons.vf(vhi):
            for heh in polygons.fh(fh):
                eh = polygons.edge_handle(heh)
                ehids.append(eh.idx())
        return set(ehids)

# I just subclass Openmesh's Trimesh / Polymesh

# Class for the dual mesh
class NPolyMesh(PolyMesh, Pymesh):
    
    # overriding pymesh to add a vertex to edge mapping (looks like I didn't do this in the end)
    def finalize(self):
        self._pythonise()
        self._helengths()
        self.normal = normal
        
    def _set_face_properties(self):
        # organise the properties we will use
        # mesh area, mesh perimeter, angular defecit
        mesh = self
        areas = OrderedDict()
        prims = OrderedDict()
        centroid = OrderedDict()

        meshpt = mesh.pym
        for fh in mesh.faces():
            area = 0.
            prim = 0.
            ag = 1. # Angular defecit. Set it to one for internal vertices
            # the mesh points
            fpts = [meshpt[mvh.idx()] for mvh in mesh.fv(fh)]
            n = len(fpts)
            for i, fpt in enumerate(fpts):
                ip = (i+1) % n
                area += np.dot(np.cross(fpt, fpts[ip]), mesh.normal)
                l = norm(fpt-fpts[ip])
                prim += l
            area = area/2.
            fhid = fh.idx()
            areas[fhid] = area
            prims[fhid] = prim
            #print 'perimeter', prim
            # centroid
            cent = self._centroid(fpts, area)
            centroid[fhid] = cent
        mesh.areas = areas
        mesh.prims = prims
        mesh.centroids = centroid

    # We don't use the cell centroids for the stress calculation. 
    def _centroid(self, fpts, area):
        csum = np.zeros(3)
        n = len(fpts)
        for i, fpt in enumerate(fpts):
            ip = (i+1) % n
            xi, yi, _ = fpt; xu, yu, _ = fpts[ip]
            cr = (xi*yu - xu*yi)
            csx= (xi + xu) * cr
            csy= (yi + yu) * cr
            csum += np.array([csx, csy, 0.])
        cent = 1/(6. * area) * csum
        #print cent
        return cent
    
# This objects construction is handled by the constructor 
# of PVmesh which is the last method of the last object in the file
class NTriMesh(TriMesh, Pymesh):
    
    #def __init__(self, *args, **kw):
        #super(NTriMesh, self).__init__(*args, **kw)

    # overwriting this method which is used in the construction of the object
    def finalize(self):
        self._pythonise()
        self._helengths()
        self.normal = normal
        self.calc_rcm()

    # the triangulation vertex ids but without the boundary 
    @property
    def bulk(self):
        # instead of calculating this from self.boundary, use bflag which is defined in the constructor at the end of the class
        if hasattr(self, '_bulk'):
            return self._bulk
        assert hasattr(self, 'bflag')
        self._bulk = np.where(self.bflag==0)[0]
        return self._bulk

    # NOT USED?
    def get_boundary_normal(self, rcm):
        bdnormals = {}
        boundary = self.boundary 

        vh = self.vertex_handle(boundary[0])
        itb = self.iterable_boundary(vh)
        heids = [heh.idx() for heh in list(itb)]

        for heid, heidprev in zip(heids, np.roll(heids, 1)):
            
            lh =self.lvec[heid]
            lhp = self.lvec[heidprev]
            bn = (lh-lhp)/norm(lh - lhp)

            vh = self.from_vertex_handle(self.halfedge_handle(heid))
            vhid = vh.idx()
            vhpt = self.pym[vhid]
            R = (vhpt- rcm)/norm(vhpt - rcm)
            flip = np.dot(bn, R) < 0 
            if flip: bn = -bn

            # scrap it
            bn = R
            bdnormals[vh.idx()] = bn
        return bdnormals


    # find the centre of mass
    # I want to do this just for the bulk and not the boundary particles
    # The centre of mass is useful for calculating radial stress profiles
    def calc_rcm(self):
        triptarr = np.array(self.pym.values())[self.bulk]
        ltript = len(triptarr)
        #print len(triptarr), len(self.pym)
        self.rcm = np.sum(triptarr, axis=0)/ltript
        self.rcmd = map(norm, triptarr - self.rcm)


    # Construct the dual mesh as an NPolyMesh (OpenMesh Polymesh object)
    # I kept adding mappings between elements of the triangulation and dual meshes as they were needed
    # Other than that the core of the code this function has been the same since approx. May 2016
    def dual(self):
        tri = self
        mesh = NPolyMesh()
        tript = tri.pym

        self.slambda = {}

        # prep for mapping dual and triangulation edges together
        # trimap
        trimap = {}

        # Trimesh vertice ids and mesh faces ids naturally match up
        ccenters = np.zeros((tri.n_faces(),3))
        cradius = np.zeros(tri.n_faces()) # for visualisation
        for j, fh in enumerate(tri.faces()):
            # Calculate circumcentres of mesh
            fhid = fh.idx()
            l_s = np.zeros(3)
            vhs = []
            # prep for mapping dual edges
            trids = [vh.idx() for vh in tri.fv(fh)]
            for vi, vj in zip(trids, np.roll(trids, -1)):
                # j is the mesh vertex id here
                trimap[(vi,vj)] = j 

            for i, heh in enumerate(tri.fh(fh)):
                lvec = tri.lvec[heh.idx()]
                l_s[i] = norm(lvec)**2
                vtmp = tri.to_vertex_handle(heh)
                vhs.append(vtmp)
                
            # match up vertices and edges 
            vhs = np.roll(vhs, -1, axis=0)
            vi, vj, vk = [tript[vh.idx()] for vh in vhs]
            lsi, lsj, lsk = l_s
            lli = lsi*(lsj + lsk - lsi)
            llj = lsj*(lsk + lsi - lsj)
            llk = lsk*(lsi + lsj - lsk)
            # actually want to save lli,llj,llk for later as a face property
            llp = OrderedDict([(vhs[0].idx(),lli), (vhs[1].idx(),llj), (vhs[2].idx(),llk)])
            self.slambda[fhid] = llp

            llnorm = lli + llj + llk
            cc = np.array(lli*vi + llj*vj + llk*vk)/llnorm
            ccenters[j] = cc
            cradius[j] = norm(cc- vi)
        self.cradius = cradius # for visualisation
        # Add cicumcenters to form a new mesh, one for each face
        mverts = np.array(np.zeros(len(ccenters)),dtype='object')
        for j, cc in enumerate(ccenters):
            mverts[j] = mesh.add_vertex(PolyMesh.Point(*cc))
        # Add mesh faces, one for each vertex
        # these mappings are extensively used to go between the triangulation and mesh
        to_mesh_face = {} # mapping between triangulation vertex ids and their cell faces 
        to_tri_vertex = {} # the opposite mapping to 'to_mesh_face'
        for vh in tri.vertices():
            if tri.is_boundary(vh):
                pass
            else: # Internal vertex
                # needs to be anticlockwise
                fhs = [f.idx() for f in list(tri.vf(vh))]
                fhs.reverse()
                vhs = list(mverts[fhs])
                mfh = mesh.add_face(vhs)
                to_mesh_face[vh.idx()] = mfh.idx()
                to_tri_vertex[mfh.idx()] = vh.idx()

        self.trimap = trimap
        # map ordered vertex pairs to halfedges
        mesh.vedgemap()

        # now map dual edges to eachother
        to_mesh_edge = {}
        
        # construct mapping between edges of each mesh
        for eh in tri.edges():
            if tri.is_boundary(eh):
                continue
            heh= tri.halfedge_handle(eh, 0)
            heho= tri.halfedge_handle(eh, 1)
            f = tri.face_handle(heh)
            fo = tri.face_handle(heho)
            mvh = mesh.vertex_handle(f.idx())
            mvho = mesh.vertex_handle(fo.idx())
            eset = [e.idx() for e in mesh.ve(mvh)]
            eseto = [e.idx() for e in mesh.ve(mvho)]
            es = set(eset).intersection(set(eseto))
            if len(es) != 1:
                 print 'no intersection'
                 print eset, eseto
            assert len(es) == 1
            dual = es.pop()
            to_mesh_edge[eh.idx()] = dual

        tri.to_mesh_edge = to_mesh_edge

        #### Gonna get rid of most of this at some point

        #mesh.bulk
        mesh.boundary

        tboundary = tri.boundary
        tb_hedges = tri.boundary_hegdes
        halfcells = OrderedDict()
        tript = tri.pym

        # construct the halfcells object which represents the boundary faces
        # Add boundary triangulation points to the mesh

        for vhid in tboundary:

            hehnid = tb_hedges[vhid] # So this is the outgoing edge
            hehn = tri.halfedge_handle(hehnid)
            hehp = tri.prev_halfedge_handle(hehn)
            # write custom function for iterating around a vertex starting from hehp
            halfcell= self.iterate_boundary_vertex(hehp)
            vh = tri.vertex_handle(vhid)
            halfcells[vhid] = halfcell

        # Assign the dictionaries to their appropriate object
        self.to_mesh_face = to_mesh_face
        mesh.to_tri_vertex = to_tri_vertex
        mesh.halfcells = halfcells
        mesh.finalize()
        return mesh

# Passing an ordered dictionary of stress values -- { triangulation vertex id : stress tensor }
#  -- to this object triggers calculation of shear and normal stress and pressure
class Vstress(object):
    def __init__(self, stress, clist, name):
        
        self.stress = stress
        assert set(clist) == set(self.stress.keys())
        self.clist = clist
        self.name = name
        self._parts()
        self.lcl = len(self.clist)

    def _parts(self):
        # normal stress and pressure are the same but with opposite sign
        self.nstress = OrderedDict()
        self.sstress = OrderedDict()
        self.pressure = OrderedDict()
        # still ignoring kinetic part for now
        stress = dict([(i, self.stress[i]) for i in self.clist])
        evalues, evectors = {}, {}
        for i in self.clist:
            st = stress[i]

            est = eig(st[:2,:2])

            evalues[i], evectors[i] = est
                
            #pressure is -ve of the average of the diagonal of the stress
            self.pressure[i] = - 1/2. * np.trace(st)
            smax, smin = evalues[i]
            if abs(smin) > abs(smax):
                smin, smax = smax, smin
                # don't lose anything by doing this right?
                evalues[i][0], evalues[i][1] = evalues[i][1], evalues[i][0] 
                evectors[i][0], evectors[i][1] = evectors[i][1], evectors[i][0] 
            smax, smin = evalues[i]
            self.nstress[i] = (smax + smin)/2
            
            self.sstress[i] = (smax - smin)/2
        # these now ordered max and then min
        self.evalues = evalues
        self.evectors = evectors
        self.avg_pressure = np.mean(self.pressure.values())
        self.avg_shear= np.mean(self.sstress.values())

    def _numpify(self):
        nstress = np.array(self.nstress.values())
        sstress = np.array(self.sstress.values())
        paxis = np.array([ev[0] for ev in self.evectors.values()])
        return nstress, sstress, paxis

    # Now we have self.stress and self.pressure the next step is to 
    #  reduce to averaged quantities
    # The radial stress profile is calculated by binning the cells based on distance from the centre of mass
    #  Why on earth didn't I just use a sliding window instead?
    def radial(self, rcm, rcmd, nstat=1000.):
        # need to adjust nstat to be a proportion of the number of cells (with stresses)
        nbins = int(m.ceil(self.lcl/nstat))
        if nbins < 4: nbins = 4
        
        print 'number of cells per bin', self.lcl/nbins
        radialq= [self.pressure, self.sstress]
        # going to name the output arrays by prepending self.name
        qn = ['radial_normal_stress', 'radial_shear_stress']
        radials = dict(zip(qn, radialq))
        #now bin the points according to their distance from the centre
        mrc = np.max(rcmd)
        rspace = np.linspace(0, mrc, nbins+1, endpoint=True)
        # the bin index of each point based on distance from rcmb
        npd = np.digitize(rcmd, rspace, right=True)-1
        # It remains to average the contents of each bin

        #assert nbins == len(rspace)-1

        lenr = len(radials)
        ravg = dict([(rname, np.zeros(nbins)) for rname in qn])
        hcount = np.zeros(nbins)
        for i, bn in zip(self.clist, npd):
            for rname, rq in radials.items():
                prs = rq[i]
                ravg[rname][bn] += prs
            hcount[bn] += 1
        for k, ravgit in ravg.items():
            ravg[k] = np.true_divide(ravgit, hcount)
            #print k, ravgit, hcount
            #print ravg[k]
        print ravg
        self.hcount = hcount
        self.rnames = qn
        self.ravg = ravg
        self.rspace =rspace

    ### some old experimental methods for getting information about the stresses
    def compare(self, vst):
        # The stress should have the same keys as self
        assert self.clist == vst.clist
        # normal and shear stress
        nstress, sstress, paxis = self._numpify()
        vnstress, vsstress, vpaxis = vst._numpify()
        ncc = np.abs(nstress - vnstress)
        scc = np.abs(sstress - vsstress)
        pax = np.arccos( [np.dot(a, v) for a, v in zip(paxis, vpaxis)]) 
        nmean, nvar = np.mean(ncc), np.var(ncc)
        smean, svar = np.mean(scc), np.var(scc)
        pmean, pvar = np.mean(pax), np.var(pax)
        print 'comparing %s with %s' % (self.name, vst.name)
        print 'average1 average2 mean variance'
        print 'normal'
        print np.mean(nstress), np.mean(vnstress), nmean, nvar
        print 'shear'
        print np.mean(sstress), np.mean(vsstress),  smean, svar
        print 'principle axis, mean variance'
        print pmean, pvar

    def orientation(self):
        base = np.array([1., 0.])
        ev = [self.evectors[i][0] for i in self.clist]
        def theta(axi):
            return np.arccos(np.dot(base, axi))
        thetas = np.array(map(theta, ev))
        plt.hist(thetas, 20)
        plt.show()


debug = False

# The main object for operating on the triangulation and dual meshes together.
# self.tri and self.mesh are the two meshes
class PVmesh(object):

    def __init__(self, tri):

        # The cell centre triangulation
        self.tri = tri

        # openmesh doesn't store edge lengths?
        if debug:
            print
            print 'trimesh'
            diagnose(self.tri)

        # The cell vertex mesh
        if debug: print 'calculating the dual'
        # Build the dual mesh
        self.mesh = self.tri.dual()

        if debug:
            print 
            print 'mesh'
            diagnose(self.mesh)

        self.mesh._set_face_properties()
        self.stresses = {}

         # want to factor all this out since it is only used for Hardy stress calculation
        self.meshes = {}
        self.meshes[0] = self.tri
        self.meshes[1] = self.mesh
        # perhaps we should construct python dictionaries for retrieving mesh points and lengths
        self.ptmesh = {}
        self.ptmesh[0] = self.tri.pym
        self.ptmesh[1] = self.mesh.pym


    # We create a mapping between halfedge ids and the lambda parameter
    # for general cell types we need to check the type of each cell adjacent to the edge.
    # NEEDS CHECKING
    def _construct_cl_dict(self, lpairs, boundary_type=1):
        # create a dictionary for the L property 
        L = lpairs['default']
        cl_dict = {}
        mesh = self.mesh; tri = self.tri
        types = tri.types
        # we should iterate through cells instead and use mapping to_mesh_face
        # we can use to_tri_vertex mapping (opposite of to_mesh_face)
        for eh in mesh.edges():
            heh = mesh.halfedge_handle(eh, 0)
            heho = mesh.halfedge_handle(eh, 1)
            mf = mesh.face_handle(heh)
            mfo = mesh.face_handle(heho)
            tp_1, tp_2 = 0, 0
            if mf.idx() == -1: # boundary halfedge
                tp_1 = boundary_type
            else:
                trivhid = mesh.to_tri_vertex[mf.idx()]
                tp_1 = types[trivhid]
            # same again
            if mfo.idx() == -1:
                tp_2 = boundary_type
            else:
                trivhido = mesh.to_tri_vertex[mfo.idx()]
                tp_2 = types[trivhido]

            lpkey = frozenset([tp_1, tp_2])
            
            l = lpairs[lpkey] if lpkey in lpairs else L
            cl_dict[heh.idx()] = l
            cl_dict[heho.idx()] = l
        return cl_dict

    # NOT USED
    def _set_angular_defecit(self):
        tri = self.tri; mesh = self.mesh
        meshpt = mesh.pym
        # set all the defecit values to 1 for bulk cells
        bthetas = OrderedDict(zip(tri.pym.keys(), np.full(tri.n_vertices(), 1.)))
        for vhid in tri.boundary:
            mvhid = tri.to_boundary_mesh_vertex[vhid]
            r_p = meshpt[mvhid]
            hc = mesh.halfcells[vhid]
            r_mu_1, r_mu_n = meshpt[hc[0]], meshpt[hc[-1]]
            r_mu_1_p = r_mu_1 - r_p
            r_mu_n_p = r_mu_n - r_p
            dtheta = np.arccos( np.dot(r_mu_1_p, r_mu_n_p) / (norm(r_mu_1_p) * norm(r_mu_n_p)) )
            sg = np.dot( np.cross(r_mu_1_p , r_mu_n_p), mesh.normal) >= 0
            ag = dtheta if sg else 2*pi - dtheta
            ag /= 2*pi
            bthetas[vhid] = ag
        tri.bthetas = bthetas

    # Set the vertex model parameters
    def set_constants(self, K, Gamma, cl_dict):
        # tri vertex id
        kproperty = OrderedDict()
        # tri vertex id
        gammaproperty = OrderedDict()
        # mesh edge id
        clproperty = OrderedDict()
        self.tri.kproperty = K
        self.tri.gammaproperty= Gamma
        self.mesh.clproperty = cl_dict

    # Iterate through the mesh vertices for any cell (even boundary)
    # return a list of vertex ids and a list of outgoing halfedge ids
    def loop(self, trivh):
        boundary = self.tri.is_boundary(trivh)
        mesh = self.mesh
        trivhid = trivh.idx()
        if not boundary:
            fhid = self.tri.to_mesh_face[trivhid]
            fh = mesh.face_handle(fhid)
            vhs = []

            hehs = [heh.idx() for heh in mesh.fh(fh)]
            for outgoing in hehs:
                heh = mesh.halfedge_handle(outgoing)
                vh = mesh.from_vertex_handle(heh)
                vhs.append(vh.idx())
        else:
            vhs = mesh.halfcells[trivhid]
            # need to find the correct halfedge here
            start_he = None
            for heh in mesh.voh(vh):
                if mesh.is_boundary_edge(heh):
                    start_he = heh
                    break
            hehs = []
            for _ in range(len(vhs)-1):
                hehs.append(start_he.idx())
                start_he = mesh.next_halfedge_handle(start_he)

        return vhs, hehs
    
    # NOT USED
    def _energy(self, vhid):
        mesh = self.mesh; tri = self.tri
        prefarea = tri.prefareas[vhid]
        #ag = tri.bthetas[vhid]
        ag = 1.
        k = tri.kproperty[vhid]
        gamma = tri.gammaproperty[vhid]
        mvhid = tri.to_mesh_face[vhid]
        area = mesh.areas[mvhid]
        perim = mesh.prims[mvhid]

        farea = k/2 * (area - ag*prefarea)**2
        fprim = gamma/2 * perim**2

        e_cl = 0.
        mfid = tri.to_mesh_face[vhid]
        mf = mesh.face_handle(mfid) 
        for heh in mesh.fh(mf):
            hehid = heh.idx()
            li = norm(mesh.lvec[hehid])
            cl = mesh.clproperty[hehid]
# IMPORTANT, we made a choice in handling the boundary edges here, don't forget!
            e_cl += 1/2. * cl * li 

        fen = farea + fprim + e_cl
        return fen

    # NOT UPDATED, DO NOT USE
    def calculate_energy(self):
        # Storing the energies on each vertex
        vhids = self.tri.bulk
        energies = map(self._energy, vhids)
        self.tri.energies = OrderedDict(zip(vhids, energies))
        tenergy = sum(energies)
        #if debug: print 'total energy', tenergy
        print 'total energy', tenergy

    # Angle defecit derivative
    # NOT USED
    def calc_dzetadr(self):
        tri = self.tri; mesh= self.mesh
        # dzetadr[boundary vertex][i, j, k vertex]
        dzetadr = {}
        tboundary = tri.boundary
        nbd = len(tboundary)
        for i, vhid in enumerate(tboundary):
            dzetadr[vhid] = {}
            vh = tri.vertex_handle(vhid)
            jm = (i-1) % nbd
            jp = (i+1) % nbd

            vhmid = tboundary[jm]
            vhpid = tboundary[jp]

            # Calculate dzetadr[i][i], setup
            hc = mesh.halfcells[vhid]
            mu1, mun = hc[0], hc[-1] # mesh vertices
            rmu1, rmun = meshpt[hc[0]], meshpt[hc[-1]]
            ri = tript[vhid]
            r1, rn = rmu1 - ri, rmun -ri
            nr1, nrn = norm(r1), norm(rn)
            agarg = np.dot(r1, rn) / (nr1 *nrn)
            sgn = -1 if np.dot( np.cross(r1 , rn), mesh.normal) >= 0. else 1
            pre_fac = sgn *  1/(2*pi) * 1/np.sqrt(1-agarg**2) 

            # derivative 
            d1 =1/( nr1 * nrn)
            d2 = np.dot(r1, rn) * d1**2
            v1a = rmmultiply(rn, drmudrp[mu1][vhid] -npI) 
            v1b = rmmultiply(r1, drmudrp[mun][vhid] -npI)
            v2a = nr1 * rmmultiply( (rn/nrn), drmudrp[mun][vhid] - npI) 
            v2b = nrn*  rmmultiply( (r1/nr1), drmudrp[mu1][vhid] - npI)
            #deriv_X = d1 * v1 - d2 * ( v2a + v2b )
            #deriv_ag = pre_fac * deriv_X
            #dzetadr[vhid][vhid] = deriv_ag

            dzetadr[vhid][vhid] = {}
            dzetadr[vhid][vhid][mu1] = pre_fac * (d1 * v1a - d2 * v2a)
            dzetadr[vhid][vhid][mun] = pre_fac * (d1 * v1b - d2 * v2b)

            # i, j where j is over nearest neighbours
            lp = tri.vv(vh)
            for nnvh in lp:
                nnvhid = nnvh.idx()
                dzetadr[vhid][nnvhid] = {}
                cell_verts = set([mu1, mun]).intersection(interloops[vhid][nnvhid])
                if bool(cell_verts) is False:
                    # This is the case where an adjacent cell shares no boundary vertices
                    continue
                deriv_X = np.zeros(3)
                for mu in cell_verts:
                    assert (mu == mu1) or (mu == mun)
                    rmu, rother = [r1, rn] if mu == mu1 else [rn, r1]
                    nrmu = norm(rmu); nrother = norm(rother)
                    dzX = ( rmmultiply( rother/(nr1 * nrn), drmudrp[mu][nnvhid])
                            - np.dot(r1,rn)/(nr1 * nrn)**2 * nrother
                            * rmmultiply( rmu/nrmu, drmudrp[mu][nnvhid] ) )
                    deriv_X = dzX

                    dzetadr[vhid][nnvhid][mu] = pre_fac * deriv_X
        return dzetadr   

    # The force calculation which reproduces forces calculated by SAMoS
    # Written Summer 2016. Last time checked against SAMoS forces, Winter 2016
    def calculate_forces(self, exclude_boundary=True):
        mesh = self.mesh
        tri = self.tri

        # Start by calculating d[lambda_i]/d[r_p] for {i,j,k} on each face p
        tlvec= tri.lvec
        tript = tri.pym; meshpt = mesh.pym
        slambda = tri.slambda
        drmudrp = {}
        for tri_fh in tri.faces():
            vh = mesh.vertex_handle(tri_fh.idx())
            #dlamqdrp {tri_vh.idx() : np(3,3) }
            dlamqdrp = {}
            drmudrp[vh.idx()] = {}
            # Get the corresponding face
            lq = np.zeros((3,3))
            lqs = np.zeros(3)
            r_q = np.zeros((3,3))
            r_qvh= []
            for i, hehq in enumerate(tri.fh(tri_fh)):
                lq[i] = tlvec[hehq.idx()]
                lqs[i] = norm(lq[i])**2
                vhi = tri.to_vertex_handle(hehq)
                r_qvh.append(vhi.idx())
            r_qvh = np.roll(r_qvh, -1, axis=0)
            r_q = np.array([tript[trivh] for trivh in r_qvh])
            rjk, rki, rij = -lq
            lis, ljs, lks =  lqs
            # Stepping towards calculating the jacobian
            # we do this for each cell around the vertices
            i = r_qvh[0]
            dlamqdrp[i] = np.zeros((3,3))
            dlamqdrp[i][0,:] = 2*lis*(-rki + rij)
            dlamqdrp[i][1,:] = -2*(lis + lks - 2*ljs)*rki + 2*ljs*rij
            dlamqdrp[i][2,:] = 2*(lis + ljs -2*lks)*rij - 2*lks*rki

            j = r_qvh[1]
            dlamqdrp[j] = np.zeros((3,3))
            dlamqdrp[j][0,:] = 2*(ljs + lks - 2*lis)*rjk - 2*lis*rij
            dlamqdrp[j][1,:] = 2*ljs*(rjk-rij)
            dlamqdrp[j][2,:] = -2*(lis + ljs -2*lks)*rij + 2*lks*rjk

            k = r_qvh[2]
            dlamqdrp[k] = np.zeros((3,3))
            dlamqdrp[k][0,:] = -2*(ljs + lks - 2*lis)*rjk + 2*lis*rki
            dlamqdrp[k][1,:] = 2*(lis + lks -2*ljs)*rki - 2*ljs*rjk
            dlamqdrp[k][2,:] = 2*lks*(-rjk + rki)

            dLambdadrp = {}
            for key, arr in dlamqdrp.items():
                dLambdadrp[key] = np.sum(arr, axis=0)
            lambdaq = slambda[tri_fh.idx()]

            gamma = sum(lambdaq.values())
            # The jacobian for each mesh vertex and adjacent face
            for p in dlamqdrp.keys():
                t1 = gamma * np.einsum('qm,qn->mn', r_q, dlamqdrp[p])
                t2 = gamma * lambdaq[p] * np.identity(3)
                lqrq = np.einsum('q,qn->n', lambdaq.values(), r_q)
                t3 = np.outer(lqrq, dLambdadrp[p])
                drmudrp[vh.idx()][p] = (1/gamma**2) * (t1 + t2 - t3)

        #alias
        loop = self.loop
        
        # Calculate loops 
        # loops {fid:set(vhids)}
        loops = {}
        for vh in tri.vertices():
            loops[vh.idx()], _ = loop(vh)
        self.loops = loops
        # Calculate loop interesctions
        # {vhi:{vhj:set(v1_idx,v2_idx..)}}
        # It's natural to use sets and the intersect() method for dealing with loops
        interloops = {}
        for vhi in tri.vertices():
            vhidx = vhi.idx()
            interloops[vhidx] = {}
            for vhj in tri.vv(vhi):
                vhjdx = vhj.idx()
                intset = set(loops[vhidx]).intersection(set(loops[vhjdx]))
                interloops[vhidx][vhjdx] = intset

        # dAdrmu[fhid][vhid]  
        # this should work for boundary cells as well
        dAdrmu = {}
        dPdrmu = {}
        dLdrmu =  {}
        # self.loop doesn't work for boundary vertices
        # wrong. There are contributions of the boundary cells to lambda.
        # even though their perimiter and areas are zero
        #for trivhid in tri.bulk:

        for trivhid in tri.pym.keys():
            trivh = tri.vertex_handle(trivhid)
            dAdrmu[trivhid] = {}
            dPdrmu[trivhid] = {}
            dLdrmu[trivhid] = {}

            # Vertices and halfedges of the loop
            boundary = tri.is_boundary(trivh)
            lp, hp = loop(trivh)
            if boundary:
                hc = mesh.halfcells[trivhid]
            nl = len(lp)
            for i, muid in enumerate(lp):
                vh = mesh.vertex_handle(muid)

                # Caluculate area and perimeter derivatives on the vertices 
                # need next and previous vertices
                ni = (i+1) % nl
                npr = (i-1) % nl
                vhplus = mesh.vertex_handle(lp[ni])
                vhminus = mesh.vertex_handle(lp[npr])

                vplus, vminus = meshpt[vhplus.idx()], meshpt[vhminus.idx()]

                crtotal = np.array([0.,0.,0.])
                crpl = 1/2. *  np.cross(vplus, mesh.normal) 
                if not boundary or (boundary and (i != nl-1)): crtotal += crpl
                crmi = -1/2. * np.cross(vminus, mesh.normal) 
                if not boundary or (boundary and (i != 0)): crtotal += crmi 
                dAdr = crmi + crpl
                dAdrmu[trivhid][muid] = dAdr

                #dPdrmu
                # get lengths
                vhpt = meshpt[muid]
                lvm = vhpt - vminus
                lvtotal = np.array([0.,0.,0.])
                if not boundary or (boundary and (i != 0)): lvtotal += lvm/norm(lvm)
                lvp = vplus - vhpt
                if not boundary or (boundary and (i != nl-1)): lvtotal -= lvp/norm(lvp)
                dPdrmu[trivhid][muid] = lvtotal

                #dLdrmu
                # for this mesh vertex what are the surrounding points in cell (mesh face) trivh 

                if boundary and i == nl-1:
                    mhehid = hp[i-1] # outgoing halfedge
                    mhehm = mesh.halfedge_handle(mhehid)
                    lmv = mesh.lvec[mhehm.idx()]
                    lmv = lmv/norm(lmv)
                    lv = np.zeros(3)
                else:
                    mhehid = hp[i] # outgoing halfedge
                    mheh = mesh.halfedge_handle(mhehid)
                    mhehm = mesh.prev_halfedge_handle(mheh)
                    lv = mesh.lvec[mheh.idx()]
                    lmv = mesh.lvec[mhehm.idx()]
                    lv = lv/norm(lv); lmv = lmv/norm(lmv)

                cl = mesh.clproperty[mheh.idx()]
                clm = mesh.clproperty[mhehm.idx()]

                # 1/2. factor here becuase each edge is iterated over twice in summing the forces
                #the commented line is actually correct, Rastko's code doesn't have this 1/2. factor
                #dLdrmu[trivhid][muid] = 1/2. * ( clm * lmv - cl * lv )

                ldrm = np.array([0.,0.,0.])
                if not boundary or (boundary and i != 0): 
                    ldrm +=  ( clm * lmv )
                if not boundary or (boundary and i != nl-1): 
                    ldrm += -1. * cl * lv 

                dLdrmu[trivhid][muid] = ldrm

        # It remains to do some complicated math over loops of nearest neighbours, etc..
        lbl =mesh.n_vertices()
        vertex_force = OrderedDict(zip(range(lbl), np.zeros((lbl, 3))))

        forces = OrderedDict()

        for trivh in tri.vertices():
            trivhid = trivh.idx()
            boundary = tri.is_boundary(trivh)
            farea, fprim, flc = [np.zeros(3)]*3
            if not (exclude_boundary and boundary): 
                mvhid = tri.to_mesh_face[trivhid]

                prefarea = tri.prefareas[trivhid]
                ag = 1.
                kp = tri.kproperty[trivhid]
                iarea = mesh.areas[mvhid]
                gammap = tri.gammaproperty[trivhid]
                prim = mesh.prims[mvhid]

                trivhpt = tript[trivh.idx()]

                # Immediate contribution
                farea_fac = -(kp) * (iarea - ag *prefarea)  
                fprim_fac = -(gammap) * prim
                asum = np.zeros(3)
                psum = np.zeros(3)
                lsum = np.zeros(3)

                for mu in loops[trivhid]: 

                    ac = rmmultiply(dAdrmu[trivhid][mu], drmudrp[mu][trivhid])
                    pc = rmmultiply(dPdrmu[trivhid][mu], drmudrp[mu][trivhid])
                    lc = rmmultiply(dLdrmu[trivhid][mu], drmudrp[mu][trivhid])

                    vff = farea_fac * dAdrmu[trivhid][mu] 
                    vfp = fprim_fac * dPdrmu[trivhid][mu]
                    vfl = -1 * dLdrmu[trivhid][mu]
                    vertex_force[mu] +=  vff + vfp + vfl

                    asum += ac
                    psum += pc
                    lsum += lc

                farea = farea_fac * asum
                fprim = fprim_fac * psum
                flc = -1 * lsum

            # And the nearest neighbours contribution
            # Some duplicated code
            area_nnsum = np.zeros(3)
            prim_nnsum = np.zeros(3)
            lc_nnsum = np.zeros(3)
            for vhnn, vidset in interloops[trivhid].items():
                nnvh = tri.vertex_handle(vhnn)
                nnboundary = tri.is_boundary(nnvh)
                if not (exclude_boundary and nnboundary): 

                    prefarea = tri.prefareas[vhnn]
                    ag = 1.
                    kp = tri.kproperty[vhnn]
                    gammap = tri.gammaproperty[vhnn]
                    mmvhid = tri.to_mesh_face[vhnn]
                    area = mesh.areas[mmvhid]
                    prim = mesh.prims[mmvhid]

                    farea_fac = -(kp) * (area - ag * prefarea)  
                    fprim_fac = -(gammap) * prim
                    asum = np.zeros(3)
                    psum = np.zeros(3)
                    lsum = np.zeros(3)

                    for mu in vidset:
                        #ac = np.einsum('n,nm->m', dAdrmu[vhnn][mu], drmudrp[mu][trivhid])
                        ac = rmmultiply(dAdrmu[vhnn][mu], drmudrp[mu][trivhid])
                        pc = rmmultiply(dPdrmu[vhnn][mu], drmudrp[mu][trivhid])
                        lc = rmmultiply(dLdrmu[vhnn][mu], drmudrp[mu][trivhid])

                        asum += ac
                        psum += pc
                        lsum += lc

                    area_nnsum += farea_fac * asum
                    prim_nnsum += fprim_fac * psum
                    lc_nnsum += -1 * lsum

            imforce = farea + fprim + flc
            nnforce = area_nnsum + prim_nnsum + lc_nnsum
            totalforce = imforce + nnforce
            forces[trivhid] = totalforce

        tri.forces = forces

        # active plus random force
        tri.factive = np.zeros((len(tri.forces),3))
        #tri.factive = tri.ftotal - np.array(tri.forces.values())

        mesh.vertex_force = vertex_force


    #########
    ### Simple Stress Calculation

    # Calculate the stress parts for each cell, but don't divide by area yet
    #  While we are at it, also calculate the texture tensor and call it 'structure'
    #  We can use this to calculate statistal strain later, as defined by F. Graner
    #  see for example http://arxiv.org/pdf/cond-mat/0301017v1.pdf 
    def makecellparts(self):
        tri, mesh= self.tri, self.mesh
        self.cellparts = {}
        np2i = np.array([[1,0,0],[0,1,0],[0,0,0]])
        structure = {}

        for vhi in tri.bulk:
            stress = np.zeros((3,3))
            prefarea = tri.prefareas[vhi]
            k = tri.kproperty[vhi]
            gamma = tri.gammaproperty[vhi]
            mvhid = tri.to_mesh_face[vhi]
            area = mesh.areas[mvhid]
            prim = mesh.prims[mvhid]

            # need to retrive the mesh halfedge
            #cl = mesh.clproperty.values()[0] # assume constant at the moment


            cortical_pressure = k * area * (area - prefarea)
            #print cortical_pressure
            stress += cortical_pressure * np2i

            vhih = tri.vertex_handle(vhi)
            struct = np.zeros((3,3))
            texture = np.zeros((3,3))
            for jheh in tri.voh(vhih):
                eh = tri.edge_handle(jheh)
                nuehid = tri.to_mesh_edge[eh.idx()]

                #nueh = tri.edge_handle(nuehid)
                #nuheh = tri.halfedge_handle(nueh,0)

                nueh = mesh.edge_handle(nuehid)
                nuheh = mesh.halfedge_handle(nueh,0)

                rnumu = mesh.lvec[nuheh.idx()]
                bondouter = np.outer(rnumu,rnumu)
                texture += bondouter
                bondouter /= norm(rnumu)
                # get the specific value of cl
                cl = mesh.clproperty[nuheh.idx()]
                #print 'cl property, ', cl
                # because samos is off by a factor of 2, lambda_sam = 2 *  lambda ( lambda we set )
                cl = 2 * cl

                # the factor of two here (cl/2.) comes from the fact that we 
                # split the the bond in to a pair, one for each cell.
                bondst = (gamma * prim + cl/2.) * bondouter
                struct += bondst  

            structure[vhi] = texture / 2.  # factor of two avoids double counting
            stress += struct 
            self.cellparts[vhi] = stress

        
            
        self.structure = structure
        return self.cellparts, structure


    # Get rings of adjacent cells, this will include boundaries although we don't want them
    # Not the most efficient use of set operation here
    # Look in this method if optimisation is necessary
    def adjavg(self, vhi, adj=0, cut_boundary=True):

        tri = self.tri
        last_ring = set([vhi])
        rings = [last_ring]
        # recursively find adjacent vertex ids
        def vvext(rings, last_ring, adj):
            if adj == 0:
                return rings 
            else:
                evss = set()
                for vhi in last_ring:
                    evss.update([v.idx() for v in tri.vv(tri.vertex_handle(vhi))])
                svhids = set([])
                for ring in rings:
                    svhids.update(ring)

                last_ring = evss - svhids

                rings.append(last_ring)
                adj -= 1
                vvext(rings, last_ring, adj)

        vvext(rings, last_ring, adj)
        vhids =[]
        for ring in rings:
            vhids.extend(ring)

        if cut_boundary:
            nuvhids = []
            for vhi in vhids:
                vh = tri.vertex_handle(vhi)
                if not tri.is_boundary(vh):
                    nuvhids.append(vhi)
            vhids = nuvhids

        return vhids

    # take a list of cells to use to calculate the virial stress
    # Just average the contributions to stress from those cells and divide by area
    def virial(self, vhids):
        tri = self.tri; mesh = self.mesh
        stress = np.zeros((3,3))
        cellparts = self.cellparts
                       
        atotal = 0.
        for vhi in vhids:
            # cut out boundary vertices, we do this in adjavg now so cut it
            vh = tri.vertex_handle(vhi)
            if tri.is_boundary(vh):
                continue
            #

            mvhid = tri.to_mesh_face[vhi]
            area = mesh.areas[mvhid]
            cp = cellparts[vhi]
            atotal += area
            stress += cp
        #print 'total area averaging over', atotal
        stress *= 1/atotal
        return stress

    # consider using  a circular region to pick appropriate cells to average over.
    # This is an alternative to recursively picking neighbours.
    #  todo,
    # IMPORTANT looks like I left this method half written and never finished it!
    # I simply preferred using the adjacent cell picking method above.
    def circle_virial(self, vhi, wl):
        # find the cells inside circle radius wl around vhi
        tri, mesh = self.tri, self.mesh
        cc = tri.pym[vhi]
        # Need a cell list, make this a property of the polymesh object
        cl = tri.clorganise(wl)
        # Should be the indices of all the nearby particles
        vhis = cl.get_neighbours(cc)
        nearvhis = []
        for vhi in vhis:
            vhpt = tri.pym[vhi]
            prox = norm(cc - vhpt)
            if prox > wl:
                continue
            nearvhis.append(vhi)

    # Use the network connectivity to find a local averaging region
    # This should work well for small regions
    # the purpose of this method is to call self.virial for every cell
    # Stresses are stored in a Vstress object
    def on_centres(self, adj=0):
        self.stressv = OrderedDict()
        for vhi in self.tri.bulk:
            vhids = self.adjavg(vhi, adj)
            self.stressv[vhi] = self.virial(vhids)

        vv = Vstress(self.stressv, self.tri.bulk, 'virial')
        vv.radial(self.tri.rcm,self.tri.rcmd)
        self.stresses['virial'] = vv

    # if we have a mapping vhids -> range(0, n) then we can use arrays instead of OrderedDicts
    def structure_on_centres(self, adj=0):
        adjstructure = OrderedDict()
        for vhi in self.tri.bulk:
            vhids = self.adjavg(vhi, adj)
            adjvalue = np.mean(np.array([self.structure[v] for v in vhids]), axis=0)

            adjstructure[vhi] = adjvalue
        return adjstructure


    #########
    # Full hardy stress calculation
    # Complicated Stress Calculation which can deal with parts of bonds 
    #  inside the averaging region as well as using (for example) the
    #  geometric centre of the polygons to triangulate them.

    # May 2017  - Currently no reason to use such a complicated solution for calculating stress.
    # Skip straight to the end of the class and forget this exists
    # Last used Autumn 2016, likely needs minor updates

    def _triangulate_areas(self, centroid):
        tri = self.tri; mesh = self.mesh
        tript = self.tri.pym; meshpt = self.mesh.pym

        polygons = NPolyMesh()
        # map the new single cell mesh vertices back onto our meshes
        # we have to keep track of which mesh our ids belong to 
        # {new_mesh_id: ('t', old_mesh_id)
        pullback = {} 
        # {('t', old_mesh_id): new_mesh_id}
        pushforward = {}
        pvh_in_tri = {}
        for mu in self.mesh.pym.keys():
            mupt = meshpt[mu]
            mv = polygons.add_vertex(TriMesh.Point(*mupt))
            pushforward[(1, mu)] = mv.idx()
            pullback[mv.idx()] = (1, mu)
        for p in self.tri.bulk:
            vhi = tri.vertex_handle(p)
            # we might want to use the centroid instead of the voronoi center
            fhid = tri.to_mesh_face[p]
            if centroid:
                ipt = mesh.centroids[fhid]
            else:
                # voronoi center
                ipt = tript[p]
            #print norm(tript[p] - mesh.centroids[fhid])

            pts = []
            # we assume this loop orders the points counter-clockwise as they should be
            loop = self.loops[p]
            lnupts = len(loop)
            # tmp
            #ipt += 

            mv = polygons.add_vertex(TriMesh.Point(*ipt))
            centerid = mv.idx()
            pvh_in_tri[centerid] = None
            pushforward[(0, p)] = centerid
            pullback[centerid] = (0, p)
            nuvertids = [pushforward[(1,nu)] for nu in loop]

            simplices = np.zeros((lnupts, 3),dtype=int)
            for i, (nu, mu) in enumerate(zip(nuvertids, np.roll(nuvertids, -1, axis=0))):
                simplices[i,0] = centerid
                simplices[i,1] = nu
                simplices[i,2] = mu
            for f in simplices:
                polygons.add_face([polygons.vertex_handle(vid) for vid in f])

        polygons.finalize()
        polypt = polygons.pym

        # Needs refactoring away but I've done enough refactoring for now
        polygon_l = {}
        for eh in polygons.edges():
            heh = polygons.halfedge_handle(eh, 0)
            av= polygons.from_vertex_handle(heh)
            bv= polygons.to_vertex_handle(heh)
            a = norm(polypt[av.idx()] - polypt[bv.idx()])
            polygon_l[frozenset([av.idx(),bv.idx()])] = a
        self.polygon_l = polygon_l

        self.polygons = polygons
        # todo think about where these dictionaries should live
        self.pushforward = pushforward
        self.pullback = pullback
        # can be part of polygons object
        polygons.pushforward = pushforward
        polygons.pullback = pullback
        polygons.pvh_in_tri = pvh_in_tri

    # can we just work on the polygons mesh instead of using a separate bonds list(?)
    def _construct_bonds(self):
        tript = self.tri.pym; meshpt = self.mesh.pym
        bonds = {}
        self.bonds = bonds
        # Find the max bond length while we are at it
        self.lmax = 0.
        self.avginubond =0.

        poly = self.polygons
        in_tri = poly.pvh_in_tri
        pullback = poly.pullback
        is_mesh_edge = {}
        avginubond = 0.
        lmaxlist = []
        for eh in poly.edges():
            heh = poly.halfedge_handle(eh, 0)
            hehid = heh.idx()
            va = poly.to_vertex_handle(heh)
            vb = poly.from_vertex_handle(heh)
            ta, ida = pullback[va.idx()]
            tb, idb = pullback[vb.idx()]
            llen = poly.ll[hehid]
            lmaxlist.append(llen)
            if ta == 0 or tb == 0:
                avginubond += llen
            else:
                is_mesh_edge[eh.idx()] = 1

            bondk = frozenset([(ta, ida), (tb, idb)])
            bonds[bondk] = eh.idx()
        mel = len(is_mesh_edge)
        # the number of new edges
        tel = poly.n_edges() - mel
        self.avginubond = avginubond / tel
        self.lmax = max(lmaxlist)
        poly.is_mesh_edge = is_mesh_edge

        if debug:
            print 'maximum bond length, ', self.lmax


    def _calculate_dEdlength(self):
        # because we split the bonds into separate dictionaries
        # we need to do the same here.
        tri = self.tri; mesh = self.mesh
        tript = self.tri.pym; meshpt = self.mesh.pym
        bonds = self.bonds
        polygons = self.polygons
        pushforward = polygons.pushforward
        pullback = polygons.pullback
        polygon_l = self.polygon_l

        dEdbond = {}
        dEspecial = {}
        for bondk in self.bonds.keys():
            dEdbond[bondk] = 0.

            dEspecial[bondk] = {}

        nc = {}
        for i in self.tri.bulk:
            vhi= tri.vertex_handle(i)

            kp = tri.kproperty[i]
            mvid= tri.to_mesh_face[i]
            area = mesh.areas[mvid]
            prefarea = tri.prefareas[i]

            pre = kp*(area - prefarea)
            # pick out triangles
            pvhi = polygons.vertex_handle(pushforward[(0, i)])
            boundary = self.tri.is_boundary(vhi)
            star = polygons.get_star_hedges(pvhi)

            for hehid in star:
                heh = polygons.halfedge_handle(hehid)
                av = polygons.from_vertex_handle(heh)
                bv = polygons.to_vertex_handle(heh)
                avid = av.idx(); bvid = bv.idx()
                avidk = pullback[avid]; bvidk = pullback[bvid]
                bondk = frozenset([avidk, bvidk])

                # add contributions from two faces if they exist
                f= polygons.face_handle(heh)

                fid = f.idx()
                vhs = [fv.idx() for fv in polygons.fv(f)]
                vhhehs = [he.idx() for he in polygons.fh(f)]

                lls = [polygons.ll[heid] for heid in vhhehs]
                dl = vhhehs.index(hehid)
                    
                # derivative with respect to the length dl
                a, b, c = np.roll(lls, -dl)
                s = (a + b + c)/2
                sa = s - a; sb = s - b; sc = s - c
                # should precalculate, todo
                cald = s*sa*sb*sc
                assert cald >= 0
                triA = np.sqrt(cald)

                preA = 1/(2. * triA)

                dDdbond = 1/2. *  ( -a*sb*sc + s*sa*sc + s*sa*sb )
                #print triA, dDdbond
                dAdbond = preA * dDdbond

                # no way to break up the derivative here according to which mesh face they belong
                dEdbond[bondk] += pre * dAdbond

                try:
                    dEspecial[bondk][i] += pre * dAdbond
                except KeyError:
                    # if [bondk][i] doesn't exist yet
                    dEspecial[bondk][i] = pre * dAdbond
            #sys.exit()
        self.dEspecial = dEspecial

        dEdbond_p = {}
        dEdbond_cl = {}
        dEspecial_p= {}
        dEspecial_cl= {}

        # construct the to_mesh_edge object
        to_mesh_edge = {}
        for eh in mesh.edges():
            heh = mesh.halfedge_handle(eh, 0)
            vto = mesh.to_vertex_handle(heh)
            vfrom = mesh.from_vertex_handle(heh)
            bondk = frozenset([(1, vto.idx()), (1, vfrom.idx())])
            try:
                pehid = bonds[bondk]
                to_mesh_edge[pehid] = eh.idx()
            except KeyError:
                # should just be a true boundary mesh edge
                pass

        self.polygons.to_mesh_edge = to_mesh_edge

        for bondk, pehid in bonds.items():
            ka, kb = bondk
            ta, ida = ka
            tb, idb = kb
            if ta == 0. or tb == 0.:
                continue # we are only concerned with cell-vertex bonds ('edges')

            ehid = to_mesh_edge[pehid]
            eh = self.mesh.edge_handle(ehid)
            heh = mesh.halfedge_handle(eh, 0)

            cl = mesh.clproperty[heh.idx()]
            # multiply by 2 because samos is still wrong here
            # Samos should take half of lambda before using it to calculate 
            # force but it doesn't so the value of lambda given to Samos is effectively
            # double what it should be lambda_sam = lambda/2 => lambda = 2 lambda_sam
            cl *= 2.

            dEdbond_cl[bondk] = 0.
            dEdbond_p[bondk] = 0.
            dEspecial_cl[bondk] = {}
            dEspecial_p[bondk] = {}
            for i in [0, 1]:
                heh = self.mesh.halfedge_handle(eh, i)
                f = self.mesh.face_handle(heh)
                fid = f.idx()
                #if this is a boundary face ignore it (?)
                #if fid == -1:
                vhi = mesh.to_tri_vertex[fid]
                vhih = tri.vertex_handle(vhi)
                if tri.is_boundary(vhih): 
                    continue
                # todo carefully consider how to manage boundaries here
                dEdbond_cl[bondk] += cl/2.
                
                gammap = tri.gammaproperty[vhi]
                mfid = tri.to_mesh_face[vhi]

                prim = mesh.prims[mfid]
                dec = gammap * prim 
                dEdbond_p[bondk] += dec

                # always one contribution
                dEspecial_p[bondk][vhi] = dec 
                dEspecial_cl[bondk][vhi] = cl/2.

            #print dEdbond_p[bondk]
            #print dEdbond_cl[bondk]

        dEdbondall = {}
        for bondk in bonds.keys():
            dEdbondall[bondk] = dEdbond[bondk]
            if bondk in dEdbond_p:
                dEdbondall[bondk] += dEdbond_p[bondk]
                dEdbondall[bondk] += dEdbond_cl[bondk]

        self.dEspecial_p = dEspecial_p
        self.dEspecial_cl = dEspecial_cl
        self.dEdbond = dEdbondall

    def _set_wl(self, wl):
        self.wl = wl
        meshpt = self.mesh.pym
        tript = self.tri.pym
        block = np.array(tript.values() + meshpt.values())
        Lx, Ly, mzero = np.amax(np.abs(block),axis=0)
        assert mzero == 0.
        cut = np.sqrt((self.lmax**2)/4. + wl**2) + 0.1
        if debug: print 'setting cut off distance to ', cut
        cl = CellList2D([2*Lx+1, 2*Ly+1], cut)
        # construct a cell list for all the vertices and centres
        # can construct separate cell lists for both
        # Then when it comes to checking intersections only do it for bonds which have both
        # ends within the neighbouring cells.
        for i, ipt in self.tri.pym.items():
            cl.add_particle(ipt[:2], (0, i))
        for nu, nupt in self.mesh.pym.items():
            cl.add_particle(nupt[:2], (1, nu))
        self.bondcl= CellList2D([2*Lx+1, 2*Ly+1], cut)
        for bondk in self.bonds.keys():
            bka, bkb = bondk
            ta, a = bka; tb, b = bkb
            pta = self.ptmesh[ta][a]
            ptb = self.ptmesh[tb][b]
            self.bondcl.add_bond(pta, ptb, bondk)
        self.cl = cl

    # All the setup required for the full hardy stress
    def _stress_setup(self, centroid=False, wl=None):
        self._triangulate_areas(centroid)
        self._construct_bonds()
        self._calculate_dEdlength()
        if wl:
            self._set_wl(wl)

    # single cell only virial stress
    def calculate_virial(self, vhi):

        # the cell vertices
        #for pvh in self.vv(self.pushforward[vhi]):
        stress = np.zeros((3,3))
        stens =  np.zeros((3,3))
        stlamb =  np.zeros((3,3))
        polygons = self.polygons

        pvhi = polygons.vertex_handle(self.pushforward[(0, vhi)])
        nc = 0
        for ehid in polygons.get_star_edges(pvhi):
            eh = polygons.edge_handle(ehid)
            heh = polygons.halfedge_handle(eh, 0)
            vha = polygons.from_vertex_handle(heh)
            vhb = polygons.to_vertex_handle(heh)
            a = vha.idx(); b = vhb.idx()
            avidk = self.pullback[a]; bvidk = self.pullback[b]
            bondk = frozenset([avidk, bvidk])
            ta, ai = avidk
            tb, bi = bvidk
            pta = self.ptmesh[ta][ai]
            ptb = self.ptmesh[tb][bi]
            bondv = ptb - pta

            # avoid including any contribution to the derivative from outside the polygon
            # mesh face id / tri vertex is vhi 
            de = self.dEspecial[bondk][vhi] 
            dtmp = 0.
            dtt = 0.
            if bondk in self.dEspecial_p:
                dtmp = self.dEspecial_p[bondk][vhi] + self.dEspecial_cl[bondk][vhi]
                de += self.dEspecial_p[bondk][vhi] 
                de += self.dEspecial_cl[bondk][vhi]

                # tmp
                dtt = self.dEspecial_cl[bondk][vhi]
            st = de * np.outer(bondv, bondv)/norm(bondv)

            sttest = dtmp * np.outer(bondv, bondv)/norm(bondv)
            stens += sttest
            stt = dtt * np.outer(bondv, bondv)/norm(bondv)
            stlamb += stt

            stress += st
            nc += 1
        
        mvhi = self.tri.to_mesh_face[vhi]
        area = self.mesh.areas[mvhi]
        #print stens
        #print stens - stlamb
        #print stlamb
        #print stress- stens
        sys.exit()
        #try:
            #self.totalarea += area
        #except:
            #self.totalarea = area
        #print self.totalarea
        # sign issues (?)
        stress = 1/area * stress
        #print stress[:2,:2]
        #print
        #print 1/2. *np.trace(stress[:2,:2])
        return stress

    def calculate_stress(self, x_c, omega, exclude=True):
        # exclude flag sets whether to return np.nan for points touching the boundary
        # wl is the smoothing length
        # x_c the point about which we calculate stress

        mesh = self.mesh; tri = self.tri
        tript, meshpt= self.tri.pym, self.mesh.pym
        bonds = self.bonds 
        lmax = self.lmax

        # Need to cut out all the bonds which are too far away to be worth checking. Efficiently.
        stress = np.zeros((3,3))
        dEdbond = self.dEdbond
        nc = 0

        # make this a positive test for bonds being nearby
        #for bondk in bonds.keys(): 
        cnan = np.full((3,3), np.nan)
        for bondk in set(self.bondcl.get_neighbours(x_c)):

            # these bonds can be repeated multiple times
            # The 't' type of bond identifies ibonds and nubonds
            bka, bkb = bondk
            ta, a = bka; tb, b = bkb
            pta = self.ptmesh[ta][a]
            ptb = self.ptmesh[tb][b]
            
            bondv = pta-ptb
            inter = bond_intersection(x_c, self.wl, pta, ptb)
            if inter is None:
                continue # couldn't find an intersection
            m_minus, m_plus, line = inter
            # do this check earlier todo
            if exclude:
                if ta == 1 and tb == 1 and (a in mesh.bulk_edge) and (b in mesh.bulk_edge):
                    stress = cnan
                    break 
            # Non-zero stress contribution below here
            bondv_hat = bondv/norm(bondv)
            # can use an arraylike representing the function and use integrate.simps
            def omegaline(m):
                return omega( norm(line(m) - x_c) )
            #io.plotrange(omegaline, m_minus, m_plus)
            omega_int, err = integrate.quad(omegaline, m_minus, m_plus)
            nc += 1
            stress += dEdbond[bondk] * omega_int * np.outer(bondv_hat, bondv)
    
        #print 'added stress with %d contributions' % nc
        if np.isnan(stress[0][0]): return stress
        #stress *= -1./(pi * wl**2)
        return stress

    def stress_on_centres(self, omega, clist=None, 
            hardy=True, virial=True, kinetic=False, quick=True):
        tript = self.tri.pym
        # By using standard numpy arrays here and adding an array for the indices in clist 
        if clist is None: clist = self.tri.bulk
        self.omega = omega
        self.omegad = []

        cll = len(self.tri.pym)
        self.clist= np.array(clist)
        self.stress = OrderedDict()
        self.stressk = OrderedDict()
        self.stressv = OrderedDict()
        
        excl= []
        for i in self.tri.bulk:
            if debug:
                if i % 100 == 0: 
                    print 'Made it to vertex', i
            ipt = tript[i]
            if hardy:
                st = self.calculate_stress(ipt, omega, exclude=False)

                if np.isnan(st[0][0]):
                    excl.append(i)
                self.stress[i] = st
            # Perhaps the most important part
            if virial:
                self.stressv[i] = self.virial(i)

            if kinetic:
                ssk = self.calculate_kinetic_stress(ipt, omega)
                self.stressk[i] = ssk
        ##print 'Made it to vertex', i
        clist = [i for i in clist if i not in excl]
        if hardy:
            vst = Vstress(self.stress, clist, 'hardy')
            vst.radial(self.tri.rcm,self.tri.rcmd)
            self.stresses['hardy'] =  vst
        if virial:
            vv = Vstress(self.stressv, self.tri.bulk, 'virial')
            vv.radial(self.tri.rcm,self.tri.rcmd)
            self.stresses['virial'] = vv

    def stress_on_vertices(self, omega, hardy = False):
        print 'calculating hardy stress around vertices'
        meshpt = self.mesh.pym
        cll = len(meshpt)
        #self._edgeparts()
        stress = OrderedDict()
        for nu in self.mesh.pym.keys():
            nupt = meshpt[nu]
            st = self.calculate_stress(nupt, omega, exclude=False)
            stress[nu] = st

            
        clist = stress.keys()
        vst = Vstress(stress, clist, 'hardy_vertices') 
        vst.radial(self.tri.rcm,self.tri.rcmd)
        self.stresses['hardy_vertices'] = vst

    def calculate_vflow(self, x_c, omega):
        vvals = self.vvals
        vflow = 0.
        for k in self.cl.get_neighbours(x_c):
            ti, i = k
            if ti == 0:
                ipt = self.tript[i]
            else: 
                continue
            nix = norm(ipt - x_c)
            vflow += omega(nix) * vvals[i]
        return vflow

    def calculate_kinetic_stress(self, x_c, omega):

        mesh = self.mesh; tri = self.tri
        tript, meshpt= self.tript, self.meshpt

        stressk = np.zeros((3,3))
        vvals = self.vvals
        # include boundaries for the moment just don't calculate stress for them
        for k in self.cl.get_neighbours(x_c):
            ti, i = k
            if ti == 0:
                ipt = self.tript[i]
            else:
                continue
            vflow = self.calculate_vflow(x_c, omega)
            vrel = vvals[i] - vflow
            # mass is 1
            stressk -= np.outer(vrel, vrel) * omega(norm(x_c -ipt))
        return stressk


    #########
    #  Simple metrics for the proximity of the system to its target area/perimeter
    #  vc_distance is the proximity of the voronoi mesh to be 'centroidal'.

    def vc_distance(self):
        # calculate the average distance between the voronoi and centroid polygon centres
        # Need to ignore boundary vertices
        tri = self.tri
        mesh = self.mesh
        vcdists = np.zeros(len(self.tri.bulk))
        for i, vhi in enumerate(self.tri.bulk):
            fhid = self.tri.to_mesh_face[vhi]
            voronoipt = tri.pym[vhi]
            centroid = mesh.centroids[fhid]
            nvc = norm(voronoipt- centroid)
            vcdists[i] = nvc
        print 'average vcdistance', np.mean(vcdists)
        return vcdists

    def areametric(self):
        mesh = self.mesh; tri = self.tri
        return sum([abs(mesh.areas[tri.to_mesh_face[i]] - tri.prefareas[i])
            for i in tri.bulk])/len(tri.bulk)
    def perimmetric(self):
        tri = self.tri; mesh = self.mesh
        prims = mesh.prims
        # assume lambda is constant
        cl = self.mesh.clproperty.values()[0]
        pmetric = 0.
        for i in tri.bulk:
            mfid = tri.to_mesh_face[i]
            prim = prims[mfid]
            ptarget = -cl/tri.gammaproperty[i]
            pmetric += abs(ptarget-prim)
        return pmetric/len(tri.bulk)


    # Autumn 2016 - attempt to align stress with texture tensor - just start over
    # Move this method to another class 
    def stress_alignment(self, thresh):
        # calculate the alignment paramter of the major principle stress axis 
        #  and the direction of particle velocity
        # A threshold value allows targeting of highly anistropic stresses
        virial= self.stresses['virial']
        # this are ordered max and then min
        evectors = virial.evectors
        clist = virial.clist
        salign = OrderedDict()
        n = len(clist)
        vels = self.tri.vvals[:,:2]
        shear = virial.sstress
        median = sorted(shear.values())[int(((len(shear)-1)*thresh))]
        #print median
        for i in clist:
            v = vels[i]
            if shear[i] >= median:
                continue
            salign[i] = np.arccos( np.dot(evectors[i][0], v/norm(v)) )
        return salign


    # Interested in organising the particles based on their types
    def _type_lists(self):
        tri = self.tri
        boundary = tri.bflag
        tplists = {}
        for i, tp in enumerate(tri.types):
            if boundary[i] == 1:
                continue # ignore boundary particles
            try:
                tplists[tp].append(i)
            except KeyError:
                tplists[tp] = [i]
        self.tplists = tplists

    # calculate the distribution of cell areas
    def type_area(self):
        mesh = self.mesh; tri = self.tri

        # have to go through mesh to get area
        def area_get(i): return mesh.areas[tri.to_mesh_face[i]]
        areatp = {}
        if hasattr(self, 'tplists'):
            for tp, tpl in self.tplists.items():
                area= map(area_get, tpl)
                areatp[tp] = area
        else:
            area = map(area_get, tri.bulk)
            areatp[1] = area
        return areatp


    ''' This is really the constructor '''
    @classmethod
    def datbuild(cls, rdat, simplices=None):
        # rdat is a ReadData object from tmp_read_data.py which is basically the same as utils/read_data
        # simplices are the triangulation triangles and are read from the faces output of SAMoS
        #  otherwise this code can use Scipy.Delaunay to independently calculate the triangulation

        keys = rdat.keys
        ids = np.array(rdat.data[keys['id']], dtype=int)
        # we put environment at the top of the file 
        startid = ids[0] 
        print 'startid = ', startid
        print 
        x = np.array(rdat.data[keys['x']])
        y = np.array(rdat.data[keys['y']])
        z = np.array(rdat.data[keys['z']])
        vx = np.array(rdat.data[keys['vx']])
        vy = np.array(rdat.data[keys['vy']])
        vz = np.array(rdat.data[keys['vz']])

        area = np.array(rdat.data[keys['area']])

        bflag = np.array(rdat.data[keys['boundary']],dtype=int)
        #rvals = np.column_stack([x, y, z])
        rvals = np.column_stack([x, y, z])
        vvals = np.column_stack([vx, vy, vz])
        if simplices is None:
            # calculating the dual faces myself
            dd = Delaunay(rvals[:,:2], qhull_options='')
            simplices = dd.simplices
        if debug:
            print 'number of cells', rvals.shape[0]
            print 'number of points in trimesh', dd.points.shape[0]

        ### Think of this half of the method as part of the constructor of the triangulation mesh object
        # New properties are justed assigned directly to the new trimesh object instead of being passed to it's constructor 
        # This properties are currently:
        #  tri.prefareas  -- the preferred cell areas
        #  tri.vvals -- the particles velocities
        #  tri.ftotal -- the forces output from SAMoS
        #  tri.bflag -- the boundary flag
        #  tri.types -- the type integer

        tri = NTriMesh()
        prefareas = OrderedDict()
        # make a mesh out of this delaunay
        mverts = OrderedDict()

        for i, idx in enumerate(ids):
            v = rvals[i]
            mv = tri.add_vertex(TriMesh.Point(*v))
            mverts[idx] = mv
            #prefareas[mv.idx()] = area[i]
            prefareas[mv.idx()] = area[i]
        for f in simplices:
            tri.add_face([mverts[i+startid] for i in f]) 
        tri.prefareas= prefareas
        tri.vvals = vvals

        if 'fx' in keys:
            fx = np.array(rdat.data[keys['fx']])
            fy = np.array(rdat.data[keys['fy']])
            fz = np.array(rdat.data[keys['fz']])
            ftotal = np.column_stack([fx, fy, fz])
            tri.ftotal = ftotal

        tri.bflag = bflag
        if 'type' in keys:
            tp = np.array(rdat.data[keys['type']])
            tri.types = tp

        # finalize() includes even more setup of the triangulation mesh object.
        # 
        tri.finalize()
        # Finally passing the triangulation mesh to create a PVmesh object builds the dual mesh
        pv = PVmesh(tri)

        # we use this value to map openmesh vertex ids onto particle ids, 
        # this was a rough way of dealing with environment particles at some point 
        pv.startid = startid
        return pv


