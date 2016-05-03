
from openmesh import *
import writemesh as wr
import ioutils as io


import numpy as np
from numpy.linalg import norm
npI = np.identity(3)
def rmmultiply(v, M):
    return np.einsum('n,nm->m', v, M)

import matplotlib.pyplot as plt
from collections import OrderedDict
from scipy.spatial import Delaunay
from read_data import ReadData


from QET.qet import QET

from math import pi 

import sys

# debugging 

from ioutils import *

def diagnose(mesh):
    print 'nedges', mesh.n_edges()
    print 'nvertices', mesh.n_vertices()
    print 'nfaces', mesh.n_faces()

# Start of feature code
class PVmesh(object):

    def __init__(self, tri):

        self.tri = tri
        # the edge mesh (dual mesh)
        self.mesh = None

        # openmesh doesn't store edge lengths?
        self.tri_lvecprop =  self._helengths(self.tri)
        print
        print 'trimesh'
        diagnose(self.tri)

        self.boundary = True
        # halfcells = {trimesh.idx() : [mesh vhids] } # should be correctly ordered
        self.halfcells = {}

        self._dual()
        print 
        print 'mesh'
        diagnose(self.mesh)

        print 'Calculating edge lengths for mesh'
        self.mesh_lvecprop = self._helengths(self.mesh)

        # Set face normal to e_z
        self.normal = np.array([0,0,1])
        #scatter(self.mesh, self.tri)

        print 'Calculating Area and perimeter'
        self._set_face_properties()

        self.s_stress = 1
        
    def _helengths(self, tri):
        # store half edge vectors as half edge property

        lvecprop = HPropHandle()
        tri.add_property(lvecprop, 'lvec')

        for eh in tri.edges():
            heh = tri.halfedge_handle(eh, 0)
            heho = tri.opposite_halfedge_handle(heh)
            vt = tri.to_vertex_handle(heh)
            vf = tri.from_vertex_handle(heh)
            npvt = omvec(tri.point(vt))
            npvf = omvec(tri.point(vf))
            vedge = npvt - npvf

            tri.set_property(lvecprop, heh, vedge)
            tri.set_property(lvecprop, heho, -vedge)
            
        return lvecprop
       

    # How about some tools for dealing with openmesh objects and triangulation
    def iterable_boundary(self, mesh, vh):
        vhidstart = vh.idx()
        def is_boundary_he(mesh, he):
            fh = mesh.face_handle(he)
            return fh.idx() is -1
        start_he = None
        # We need to find the boundary one specifically
        for heh in mesh.voh(vh):
            if is_boundary_he(mesh, heh):
                #print 'found boundary edge for vhid', vh.idx()
                start_he = heh
        assert start_he is not None
        he_i = start_he
        def iterable(he_i, start_he, mesh):
            while True:
                yield he_i
                he_i = mesh.next_halfedge_handle(he_i)
                if he_i.idx() is start_he.idx():
                    raise StopIteration
        return iterable(he_i, start_he, mesh)

    def iterate_boundary_vertex(self, tri, hehp):
        assert tri.is_boundary(hehp)
        facels = []
        while True:
            heho = tri.opposite_halfedge_handle(hehp)
            if tri.is_boundary(heho):
                break
            facels.append(tri.face_handle(heho).idx())
            hehp = tri.prev_halfedge_handle(heho)
        return facels

    def _dual(self):
        print 'calculating the dual'
        self.mesh = PolyMesh()

        lambda_prop = FPropHandle()
        self.tri.add_property(lambda_prop, 'lambda')
        self.lambda_prop = lambda_prop

        self.boundary_prop = VPropHandle()
        self.tri.add_property(self.boundary_prop, 'boundary')

        self.boundaries = [] # [[]]
        # add list as well to keep track of which boundary halfedge is outgoing from the vertex
        self.b_nheh = []

        # Trimesh vertice ids and mesh faces ids naturally match up
        ccenters = np.zeros((self.tri.n_faces(),3))
        cradius = np.zeros(self.tri.n_faces())
        for j, fh in enumerate(self.tri.faces()):
            # Calculate circumcentres of mesh
            l_s = np.zeros(3)
            vhs = []
            for i, heh in enumerate(self.tri.fh(fh)):
                lvec =self.tri.property(self.tri_lvecprop, heh)
                l_s[i] = norm(lvec)**2
                vtmp = self.tri.to_vertex_handle(heh)
                vhs.append(vtmp)
                
            # match up vertices and edges 
            vhs = np.roll(vhs, -1, axis=0)
            vi, vj, vk = [omvec(self.tri.point(vh)) for vh in vhs]
            lsi, lsj, lsk = l_s
            lli = lsi*(lsj + lsk - lsi)
            llj = lsj*(lsk + lsi - lsj)
            llk = lsk*(lsi + lsj - lsk)
            # actually want to save lli,llj,llk for later as a face property
            llp = OrderedDict([(vhs[0].idx(),lli), (vhs[1].idx(),llj), (vhs[2].idx(),llk)])
            self.tri.set_property(lambda_prop, fh, llp)

            llnorm = lli + llj + llk
            cc = np.array(lli*vi + llj*vj + llk*vk)/llnorm
            ccenters[j] = cc
            cradius[j] = norm(cc- vi)
        self.cradius = cradius # for visualisation
        # Add cicumcenters to form a new mesh
        tri_prop= VPropHandle()
        self.mesh.add_property(tri_prop, 'tri')
        self.tri_prop = tri_prop
        mverts = np.array(np.zeros(len(ccenters)),dtype='object')
        for j, cc in enumerate(ccenters):
            mverts[j] = self.mesh.add_vertex(PolyMesh.Point(*cc))
            self.mesh.set_property(tri_prop, mverts[j], 0)
        # Add faces
        self.n_tri_internal = 0
        self.n_tri_boundary = 0
        bval = 0
        for vh in self.tri.vertices():
            if self.tri.is_boundary(vh):
                is_new = True not in map(lambda boundary: vh.idx() in boundary, self.boundaries)
                if is_new:
                    bval += 1
                    self.boundaries.append([])
                    self.b_nheh.append([])
                    for heh in self.iterable_boundary(self.tri, vh):
                        self.boundaries[-1].append( self.tri.from_vertex_handle(heh).idx() )
                        self.b_nheh[-1].append(heh)
                self.tri.set_property(self.boundary_prop, vh, bval)
                self.n_tri_boundary += 1
            else: # Internal vertex
                # construct the face by stepping around the halfedges
                heh = self.tri.halfedge_handle(vh)
                fhs = [ self.tri.face_handle(heh).idx() ]
                start = heh.idx()
                while True:
                    hehp = self.tri.prev_halfedge_handle(heh)
                    heho = self.tri.opposite_halfedge_handle(hehp)
                    if heho.idx() == start:
                        break
                    fh_this = self.tri.face_handle(heho)
                    fhs.append(fh_this.idx())
                    heh = heho

                vhs = list(mverts[fhs])
                mfh = self.mesh.add_face(vhs)
                self.tri.set_property(self.boundary_prop, vh, 0)
                self.n_tri_internal += 1

        # construct the halfcells object which represents the boundary faces
        # Add boundary triangulation points to the mesh
        triverts = np.array(np.zeros(len(self.halfcells)),dtype='object')
        self.vh_mf = {}
        self.btv_bmv = {}
        for i, boundary in enumerate(self.boundaries):
            for j, vhid in enumerate(boundary):
                vh = self.tri.vertex_handle(vhid)
                hehn = self.b_nheh[i][j]
                hehp = self.tri.prev_halfedge_handle(hehn)
                # write custom function for iterating around a vertex starting from hehp
                halfcell= self.iterate_boundary_vertex(self.tri, hehp)
                self.halfcells[vhid] = halfcell

                tript = omvec(self.tri.point(self.tri.vertex_handle(vhid)))
                mvh = self.mesh.add_vertex(PolyMesh.Point(*tript))
                self.mesh.set_property(tri_prop, mvh, 1)

                hcf = list(mverts[halfcell])
                hcf.append(mvh)
                mfh = self.mesh.add_face(hcf)

                self.btv_bmv[vhid] = mvh.idx()
                self.vh_mf[vhid] = mfh.idx()

                

    def _set_face_properties(self):
        # organise the properties we will use
        mesh = self.mesh
        tri = self.tri
        areaprop = VPropHandle()
        self.tri.add_property(areaprop, 'area')
        primprop = VPropHandle()
        self.tri.add_property(primprop, 'prim')
        self.btheta_prop = VPropHandle()
        self.tri.add_property(self.btheta_prop, 'angular defecit')

        for vh in tri.vertices():
            vhid = vh.idx()
            area = 0.
            prim = 0.
            ag = 1. # Angular defecit. Set it to one for internal vertices
            boundary = tri.is_boundary(vh)
            if not boundary:
                fh = mesh.face_handle(vh.idx())
                for heh in mesh.fh(fh):
                    # mirror this with the boundary calculation and forget about lvecprop
                    rmu = omvec(mesh.point(mesh.from_vertex_handle(heh)))
                    rmup = omvec(mesh.point(mesh.to_vertex_handle(heh)))
                    area += np.dot(np.cross(rmu, rmup), self.normal)
                    ap=  np.dot(np.cross(rmu, rmup), self.normal)
                    lvec = mesh.property(self.mesh_lvecprop, heh)
                    l = norm(lvec)
                    prim += l
                area = area/2
            else:
                # calculate area for boundary vertices
                hcell = self.halfcells[vhid]
                hcellpts = map(lambda x: omvec(self.mesh.point(mesh.vertex_handle(x))), hcell) 
                r_i = omvec(self.tri.point(vh))
                hcellpts.append(r_i)
                n = len(hcellpts)
                for i, hpt in enumerate(hcellpts):
                    ip = (i+1) % n
                    area += np.dot(np.cross(hpt, hcellpts[ip]), self.normal)
                    l = norm(hpt-hcellpts[ip])
                    prim += l
                area = area/2

                # Calculate angular defecit for boundary vertices
                # Could store real angular defecit here or the updated target areas
                # Go with storing angular defecit
                r_p = omvec(tri.point(vh))
                hc = self.halfcells[vh.idx()]
                r_mu_1, r_mu_n = idtopt(mesh, hc[0]), idtopt(mesh, hc[-1])
                r_mu_1_p = r_mu_1 - r_p
                r_mu_n_p = r_mu_n - r_p
                dtheta = np.arccos( np.dot(r_mu_1_p, r_mu_n_p) / (norm(r_mu_1_p) * norm(r_mu_n_p)) )
                sg = np.dot( np.cross(r_mu_1_p , r_mu_n_p), self.normal) >= 0
                ag = dtheta if sg else 2*pi - dtheta
                ag /= 2*pi

            self.tri.set_property(self.btheta_prop, vh, ag)
            #print 'setting angular defecit of', ag
            #print 'setting area of', area
            #print 'setting perimeter of', prim
            self.tri.set_property(areaprop, vh, area)
            self.tri.set_property(primprop, vh, prim)

        self.areaprop = areaprop
        self.primprop = primprop

    def _construct_cl_dict(self, L):
        # create a dictionary for the L property 
        # this code sets values for every edge however boundary->boundary edges do not correspond to a cell edge
        mesh =self.mesh
        cl_dict = {}
        for eh in self.mesh.edges():
            heh = self.mesh.halfedge_handle(eh, 0)
            vhi = self.mesh.to_vertex_handle(heh)
            vhj = self.mesh.from_vertex_handle(heh)
            vhidx = vhi.idx(); vhjdx = vhj.idx()
            cl_dict[(vhidx, vhjdx)] = L 
            cl_dict[(vhjdx, vhidx)] = L 
        return cl_dict

    def set_constants(self, K, Gamma, cl_dict):
        # Constants for Area and Perimeter should be in .dat file? Let them be set manually.
        kprop = VPropHandle()
        self.tri.add_property(kprop, 'K')
        gammaprop = VPropHandle()
        self.tri.add_property(gammaprop, 'Gamma')
        # Need to start thinking in terms of edges for the contact length property
        clprop = EPropHandle()
        self.mesh.add_property(clprop, 'Contact Length')

        for i, (ki, gi) in enumerate(zip(K, Gamma)):
            vh = self.tri.vertex_handle(i)
            self.tri.set_property(kprop, vh, ki)
            self.tri.set_property(gammaprop, vh, gi)
        # Edges are defined between faces so need to give a dictionary of constants for generality
        for eh in self.mesh.edges():
            heh = self.mesh.halfedge_handle(eh, 0)
            vhi =  self.mesh.to_vertex_handle(heh)
            vhj = self.mesh.from_vertex_handle(heh)
            vhidx = vhi.idx(); vhjdx = vhj.idx()
            self.mesh.set_property(clprop, eh, cl_dict[(vhidx, vhjdx)])

        self.kprop = kprop
        self.gammaprop = gammaprop
        self.clprop = clprop


    # Iterate through the mesh vertices for any cell (even boundary)
    def loop(self, trivh, halfedges=False):
        tri = self.tri; mesh = self.mesh
        trivhid = trivh.idx()
        boundary = tri.is_boundary(trivh)
        vhs = []
        hehs = []
        if not boundary:
            fh = mesh.face_handle(trivhid)
        else:
            fh = mesh.face_handle(self.vh_mf[trivhid])

        for mheh in mesh.fh(fh):
            vh = mesh.from_vertex_handle(mheh)
            vhs.append(vh.idx())
            hehs.append(mheh.idx())
        if halfedges:
            return vhs, hehs
        else:
            return vhs


    def calculate_energy(self):
        # And attach values as a property to the faces of the mesh (inner vertices of triangulation)
        mesh = self.mesh
        tri = self.tri

        # might as well store the energy associated with each face
        enprop = VPropHandle()
        tri.add_property(enprop, 'energy')

        # Iterate over faces
        print 'calculating energies for each face'
        tenergy = 0
        for vh in tri.vertices():
            prefarea = tri.property( self.prefareaprop, vh )
            ag = tri.property( self.btheta_prop,vh)
            area = self.tri.property(self.areaprop, vh)
            k = self.tri.property(self.kprop, vh)

            farea = k/2 * (area - ag*prefarea)**2
            gamma = tri.property(self.gammaprop, vh)
            perim = tri.property(self.primprop, vh)
            fprim = gamma/2 * perim**2


            #do we need the energy of each cell independently?
            e_cl = 0.
            # Add internal idx to vh_mf to avoid this check
            m_face_id = self.vh_mf[vh.idx()] if tri.is_boundary(vh) else vh.idx()  
            m_face = mesh.face_handle(m_face_id) 
            for heh in mesh.fh(m_face):
                #print m_face.idx(), 
                #print mesh.is_boundary(m_face)
                li = norm(mesh.property(self.mesh_lvecprop, heh))
                eh= mesh.edge_handle(heh)
                cl= mesh.property(self.clprop, eh)
                # IMPORTANT, we made a choice in handling the extra boundary edges here, don't forget!
                e_cl += 1/2. * cl * li 

            #print e_cl
            fen = farea + fprim + e_cl
            tri.set_property(enprop, vh, fen)
            tenergy += fen
        print 'total energy of mesh', 
        print tenergy

        self.enprop = enprop

    def calculate_forces(self):
        # Convenience
        mesh = self.mesh
        tri = self.tri
        lvecprop = self.tri_lvecprop

        # Maybe start by calculating d[lambda_i]/d[r_p] for {i,j,k} on each face p
        # drmudrp[mesh_vhid, mesh_fhid, :]
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
            for i, hehq in enumerate(self.tri.fh(tri_fh)):
                lq[i] = tri.property(lvecprop, hehq)
                lqs[i] = norm(lq[i])**2
                vhi = tri.to_vertex_handle(hehq)
                r_qvh.append(vhi.idx())
            
            r_qvh = np.roll(r_qvh, -1, axis=0)
            r_q = np.array([omvec(self.tri.point(self.tri.vertex_handle(trivh))) for trivh in r_qvh])
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
            lambdaq = tri.property(self.lambda_prop, tri_fh)

            gamma = sum(lambdaq.values())
            # The jacobian for each mesh vertex and adjacent face
            for p in dlamqdrp.keys():
                t1 = gamma * np.einsum('qm,qn->mn', r_q, dlamqdrp[p])
                t2 = gamma * lambdaq[p] * np.identity(3)
                lqrq = np.einsum('q,qn->n', lambdaq.values(), r_q)
                t3 = np.outer(lqrq, dLambdadrp[p])
                drmudrp[vh.idx()][p] = (1/gamma**2) * (t1 + t2 - t3)

        # For calculating derivative of area for boundary cells 
        for boundary in self.boundaries:
            for j, vhid in enumerate(boundary):
                mvhid = self.btv_bmv[vhid]
                drmudrp[mvhid] = {}
                drmudrp[mvhid][vhid] = np.identity(3)

        #alias
        loop = self.loop

        # Nearest Neighbour faces
        # ( by vertex )

        def nnfaces(trivh):
            fhs = []
            for heh in tri.voh(trivh):
                vh = tri.to_vertex_handle(heh)
                fhs.append(vh)
            return fhs

        # Calculate loops 
        # loops {fid:set(vhids)}
        loops = {}
        for vh in tri.vertices():
            loops[vh.idx()] = loop(vh)
        # Calculate loop interesctions
        # {vhi:{vhj:set(v1_idx,v2_idx..)}}
        # It's natural to use sets and the intersect() method for dealing with loops
        interloops = {}
        for vhi in tri.vertices():
            vhidx = vhi.idx()
            interloops[vhidx] = {}
            for vhj in nnfaces(vhi):
                vhjdx = vhj.idx()
                intset = set(loops[vhidx]).intersection(set(loops[vhjdx]))
                interloops[vhidx][vhjdx] = intset


        # dAdrmu[fhid][vhid]  
        # this should work for boundary cells as well
        dAdrmu = {}
        dPdrmu = {}
        dLdrmu =  {}
        for trivh in tri.vertices():
            trivhid = trivh.idx()
            dAdrmu[trivhid] = {}
            dPdrmu[trivhid] = {}
            dLdrmu[trivhid] = {}

            # These are the vertices of the loop
            # how about the halfedges
            lp, hp = loop(trivh, halfedges=True)
            nl = len(lp)
            for i, vhid in enumerate(lp):
                vh = mesh.vertex_handle(vhid)

                # Caluculate area and perimeter derivatives on the vertices 
                # need next and previous vertices
                ni = (i+1) % nl
                npr = (i-1) % nl
                vhplus = mesh.vertex_handle(lp[ni])
                vhminus = mesh.vertex_handle(lp[npr])

                vplus, vminus = omvec(mesh.point(vhplus)), omvec(mesh.point(vhminus))
                dAdrmu[trivhid][vhid] = 1/2. * ( np.cross(vplus, self.normal) 
                        - np.cross(vminus, self.normal) )
                #dPdrmu
                # get lengths
                vhpt = omvec(mesh.point(vh))
                lvm = vhpt - vminus
                lvp = vplus - vhpt
                dPdrmu[trivhid][vhid] = lvm/norm(lvm) - lvp/norm(lvp)

                #dLdrmu
                # for this mesh vertex what are the surrounding points in cell (mesh face) trivh 

                mhehid = hp[i] # outgoing halfedge
                mheh = mesh.halfedge_handle(mhehid)
                mhehm = mesh.prev_halfedge_handle(mheh)
                lv = mesh.property(self.mesh_lvecprop, mheh)
                lmv = mesh.property(self.mesh_lvecprop, mhehm)
                lv = lv/norm(lv); lmv = lmv/norm(lmv)

                cl = mesh.property(self.clprop, mesh.edge_handle(mheh))
                clm = mesh.property(self.clprop, mesh.edge_handle(mhehm))
                # 1/2. factor here becuase each edge is iterated over twice in summing the forces
                dLdrmu[trivhid][vhid] = 1/2. * ( clm * lmv - cl * lv )

        fprop = VPropHandle()
        tri.add_property(fprop, 'force')
        self.fprop = fprop

        # The duplicated code involved in determining the angle defecit derivative
        #  for the primary and adjacent faces
        def setup_rmu(vhid):
            hc = self.halfcells[vhid]
            mu1, mun = hc[0], hc[-1] # mesh vertices
            rmu1, rmun = idtopt(self.mesh, hc[0]), idtopt(mesh, hc[-1])
            ri = idtopt(self.tri, vhid)
            r1, rn = rmu1 - ri, rmun -ri
            nr1, nrn = norm(r1), norm(rn)
            agarg = np.dot(r1, rn) / (nr1 *nrn)
            sgn = -1 if np.dot( np.cross(r1 , rn), self.normal) >= 0. else 1
            pre_fac = sgn *  1/(2*pi) * 1/np.sqrt(1-agarg**2) 
            return mu1, mun, rmu1, rmun, nr1, nrn, r1, rn, pre_fac

        # Angle defecit derivative
        # dzetadr[boundary vertex][i, j, k vertex]
        dzetadr = {}
        for bd in self.boundaries:
            nbd = len(bd)
            for i, vhid in enumerate(bd):
                dzetadr[vhid] = {}
                vh = tri.vertex_handle(vhid)
                #ag = self.tri.property(self.btheta_prop, vh)
                jm = (i-1) % nbd
                jp = (i+1) % nbd

                vhmid = bd[jm]
                vhpid = bd[jp]
                # Calculate dzetadr[i][i], setup
                mu1, mun, rmu1, rmun, nr1, nrn, r1, rn, pre_fac = setup_rmu(vhid)
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
                lp = nnfaces(vh)
                nl = len(lp)
                for j, nnvh in enumerate(lp):
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

                        #assert mu not in dzetadr[vhid][nnvhid] # should be unique set of keys
                        #print 'assigning', vhid, nnvhid, mu
                        dzetadr[vhid][nnvhid][mu] = pre_fac * deriv_X
                    #dzetadr[vhid][nnvhid] = pre_fac * deriv_X
                    
        # It remains to do some complicated math over loops of nearest neighbours, etc..
        # We also want to calculate stress now
        # {vhid: s_arr}
        a_stress = OrderedDict()
        p_stress = OrderedDict()
        l_stress = OrderedDict()
        self.stress = stress = OrderedDict()
        
        for trivh in tri.vertices():
            fh = self.mesh.face_handle(trivh.idx())
            trivhid = trivh.idx()
            fhidx = fh.idx()
            # fhidx = trimesh vertices id
            kp = tri.property(self.kprop, trivh)
            gammap = tri.property(self.gammaprop, trivh)
            area = tri.property(self.areaprop, trivh)
            prim = tri.property(self.primprop, trivh)
            prefarea = self.tri.property(self.prefareaprop, trivh)
            ag = self.tri.property(self.btheta_prop, trivh)

            trivhpt = omvec(self.tri.point(trivh))
            stress[trivhid] = np.zeros((3,3))
            a_stress[trivhid] = np.zeros((3,3))
            p_stress[trivhid] = np.zeros((3,3))
            l_stress[trivhid] = np.zeros((3,3))
            boundary = tri.is_boundary(trivh)

            # Immediate contribution
            farea_fac = -(kp) * (area - ag *prefarea)  
            fprim_fac = -(gammap) * prim
            asum = np.zeros(3)
            psum = np.zeros(3)
            lsum = np.zeros(3)

            # dAdrp and dPdrp
            for mu in loops[trivhid]:
                r_mu_vh = idtopt(mesh, mu) - trivhpt # precalculate these
                ac = rmmultiply(dAdrmu[trivhid][mu], drmudrp[mu][fhidx])
                pc = rmmultiply(dPdrmu[trivhid][mu], drmudrp[mu][fhidx])
                lc = rmmultiply(dLdrmu[trivhid][mu], drmudrp[mu][fhidx])
                if boundary: # angle defecit contribution
                    if mu in dzetadr[trivhid][trivhid]:
                        zetat = dzetadr[trivhid][trivhid][mu] * prefarea
                        ac -= zetat

                a_stress[trivhid] += farea_fac * np.outer(ac, r_mu_vh)
                p_stress[trivhid] += fprim_fac * np.outer(pc, r_mu_vh)
                l_stress[trivhid] += -1 * np.outer(lc, r_mu_vh)
                asum += ac
                psum += pc
                lsum += lc

            farea = farea_fac * asum
            fprim = fprim_fac * psum
            # no extra factor
            flc = -1 * lsum

            # And the nearest neighbours contribution
            # Some duplicated code
            area_nnsum = np.zeros(3)
            prim_nnsum = np.zeros(3)
            lc_nnsum = np.zeros(3)
            for vhnn, vidset in interloops[trivhid].items():
                nnvh = tri.vertex_handle(vhnn)

                kp = tri.property(self.kprop, nnvh)
                gammap = tri.property(self.gammaprop, nnvh)
                area = tri.property(self.areaprop, nnvh)
                prim = tri.property(self.primprop, nnvh)
                prefarea = tri.property(self.prefareaprop, nnvh)
                ag = self.tri.property(self.btheta_prop, nnvh)

                farea_fac = -(kp) * (area - ag * prefarea)  
                fprim_fac = -(gammap) * prim
                asum = np.zeros(3)
                psum = np.zeros(3)
                lsum = np.zeros(3)

                nnboundary = tri.is_boundary(nnvh)
                for mu in vidset:
                    r_mu_vh = idtopt(mesh, mu) - trivhpt
                    ac = np.einsum('n,nm->m', dAdrmu[vhnn][mu], drmudrp[mu][fhidx])
                    pc = np.einsum('n,nm->m', dPdrmu[vhnn][mu], drmudrp[mu][fhidx])
                    lc = np.einsum('n,nm->m', dLdrmu[vhnn][mu], drmudrp[mu][fhidx])
                    if nnboundary:
                        vht = mu in dzetadr[vhnn][trivhid]
                        if vht:
                            # the angle defecit contribution
                            zetat = dzetadr[vhnn][trivhid][mu] * prefarea
                            ac -= zetat
                    asum += ac
                    psum += pc
                    lsum += lc
                    a_stress[trivhid] += farea_fac  * np.outer(ac, r_mu_vh)
                    p_stress[trivhid] += fprim_fac * np.outer(pc, r_mu_vh)
                    l_stress[trivhid] += -1 * np.outer(lc, r_mu_vh)

                area_nnsum += farea_fac * asum
                prim_nnsum += fprim_fac * psum
                lc_nnsum += -1 * lsum

            imforce = farea + fprim + flc
            nnforce = area_nnsum + prim_nnsum + lc_nnsum
            totalforce = imforce + nnforce
            tri.set_property(fprop, trivh, totalforce)

            stress[trivhid] = a_stress[trivhid] + p_stress[trivhid]



    def out_force_energy(self, outfile):
        # Still aren't explicitely handling the particle ids as we should be (?)

        outfe = OrderedDict()
        ids, energy, fx, fy = [], [], [], []
        for vh in self.tri.vertices():
            #if self.tri.is_boundary(vh) and self.boundary:
                #continue
            ids.append(vh.idx())
            energy.append(self.tri.property(self.enprop, vh))
            fvh =  self.tri.property(self.fprop, vh)
            fxx, fyy, _ = fvh
            fx.append(fxx)
            fy.append(fyy)
        outfe['id'] = ids
        outfe['energy'] = energy
        outfe['fx'] = fx
        outfe['fy'] = fy

        with open(outfile, 'w') as fo:
            print 'saving force and energy to ', fo.name
            io.dump(outfe, fo)

    ''' This is really the constructor '''
    @classmethod
    def datbuild(cls, rdat):
        keys = rdat.keys
        x = np.array(rdat.data[keys['x']])
        y = np.array(rdat.data[keys['y']])
        z = np.array(rdat.data[keys['z']])
        area = np.array(rdat.data[keys['area']])
        #rvals = np.column_stack([x, y, z])
        rvals = np.column_stack([x, y, z])
        print 'number of cells', rvals.shape[0]
        dd = Delaunay(rvals[:,:2], qhull_options='')
        print 'number of points in trimesh', dd.points.shape[0]

        tri = TriMesh()
        prefareaprop = VPropHandle()
        tri.add_property(prefareaprop, 'prefarea')
        # make a mesh out of this delaunay
        mverts = np.array(np.zeros(len(rvals)),dtype='object')
        for i, v  in enumerate(rvals):
            mv = tri.add_vertex(TriMesh.Point(*v))
            mverts[i] = mv
            tri.set_property(prefareaprop, mv, area[i])
        for f in dd.simplices:
            tri.add_face(list(mverts[f])) 
        PV = PVmesh(tri)
        # We cheekily add this to the initialisation to avoid recovering it later
        PV.prefareaprop = prefareaprop
        return PV

if __name__=='__main__':

    print '-----------------------'
    print 'If you just want to run the analysis then use commander.py'
    print '-----------------------'

    epidat = '/home/dan/cells/run/soft_rpatch/epithelial_equilini.dat'

    import argparse
    import os.path as path

    parser = argparse.ArgumentParser()
 
    parser.add_argument("-i", "--input", type=str, default=epidat, help="Input dat file")
    parser.add_argument("-d", "--dir", type=str, default='/home/dan/tmp/', help="Output directory")
    #parser.add_argument("-o", "--output", type=str, default=epidat, help="Input dat file")

    args = parser.parse_args()
    epidat = args.input


    rdat = ReadData(epidat)
    pv = PVmesh.datbuild(rdat)

    # Handle K, and Gamma
    #nf = pv.mesh.n_faces()
    nf = pv.tri.n_vertices()
    # Could easily read these from a .conf file
    k = 1.
    gamma = 0.
    L = 0.
    K = np.full(nf, k)
    Gamma = np.full(nf, gamma)

    cl = pv._construct_cl_dict(L)
    pv.set_constants(K, Gamma, cl)

    pv.calculate_energy()
    pv.calculate_forces()

    outdir = args.dir


    mout = path.join(outdir, 'cellmesh.vtp')
    wr.writemeshenergy(pv, mout)

    tout = path.join(outdir, 'trimesh.vtp')
    wr.writetriforce(pv, tout)

    sout = path.join(outdir, 'stresses.vtp')
    wr.write_stress_ellipses(pv, sout)

    # Dump the force and energy
    fef = path.join(outdir, 'force_energy.dat')
    pv.out_force_energy(fef)



