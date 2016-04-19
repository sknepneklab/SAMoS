
from openmesh import *

import numpy as np
from numpy.linalg import  norm
import matplotlib.pyplot as plt
from collections import OrderedDict
from scipy.spatial import Delaunay
from read_data import ReadData

from QET.qet import QET

from math import pi 

import sys

# debugging 

def dirk(A):
    print A
    print dir(A)
    sys.exit()
def shiv(al):
    for a in al:
        print a
        print eval(a)


#np.set_printoptions(threshold=np.nan)

import contextlib
import cStringIO
@contextlib.contextmanager
def nostdout():
    save_stdout = sys.stdout
    sys.stdout = cStringIO.StringIO()
    yield
    sys.stdout = save_stdout

# want to print a vector object
def dumpvec(vec):
    print omvec(vec)

def diagnose(mesh):
    print 'nedges', mesh.n_edges()
    print 'nvertices', mesh.n_vertices()
    print 'nfaces', mesh.n_faces()

def scatter(mesh, mesh2):
    arrl = []
    for vh in mesh.vertices():
        arrl.append(omvec(mesh.point(vh)))
    npts = np.column_stack(arrl)
    x = npts[:][0]
    y = npts[:][1]
    plt.scatter(x, y, color='red')

    mesh = mesh2
    arrl = []
    for vh in mesh.vertices():
        arrl.append(omvec(mesh.point(vh)))
    npts = np.column_stack(arrl)
    x = npts[:][0]
    y = npts[:][1]
    plt.scatter(x, y, color='blue')

    plt.show()


# Start of feature code

# openmesh has a vector object
# too lazy to use this to convert to numpy arrays
# fixed to three dimensions...
def omvec(vec):
    return np.array([vec[0], vec[1], vec[2]])

def idtopt(mesh, rmuid):
    return omvec(mesh.point(mesh.vertex_handle(rmuid)))

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
        # might have to switch from storing things on mesh faces to storing them on trimesh vertices. I.e. area perimeter
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
        # Add cicumcenters to form a new mesh
        mverts = np.array(np.zeros(len(ccenters)),dtype='object')
        for j, cc in enumerate(ccenters):
            mverts[j] = self.mesh.add_vertex(PolyMesh.Point(*cc))
        # Add faces
        bval = 0
        for vh in self.tri.vertices():
            if self.tri.is_boundary(vh):
                is_new = True not in map(lambda boundary: vh.idx() in boundary, self.boundaries)
                if is_new:
                    bval += 1
                    print 'new boundary'
                    self.boundaries.append([])
                    self.b_nheh.append([])
                    for heh in self.iterable_boundary(self.tri, vh):
                        self.boundaries[-1].append( self.tri.from_vertex_handle(heh).idx() )
                        self.b_nheh[-1].append(heh)
                self.tri.set_property(self.boundary_prop, vh, bval)
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
                    #print ap
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

    def set_constants(self, K, Gamma):
        # Putting aside the contact constants for now
        # Constants for Area and Perimeter should be in .dat file? Let them be set manually.
        # We shall use properties of the mesh again to keep everything consistent
        mesh = self.mesh
        kprop = VPropHandle()
        self.tri.add_property(kprop, 'K')
        gammaprop = VPropHandle()
        self.tri.add_property(gammaprop, 'Gamma')
        #for fh in mesh.faces():
        for i, (ki, gi) in enumerate(zip(K, Gamma)):
            vh = self.tri.vertex_handle(i)
            self.tri.set_property(kprop, vh, ki)
            self.tri.set_property(gammaprop, vh, gi)

        self.kprop = kprop
        self.gammaprop = gammaprop


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
            fen = farea + fprim
            tri.set_property(enprop, vh, fen)
            tenergy += fen
        print 'total energy of mesh', 
        print tenergy

        self.enprop = enprop

    def calculate_forces(self):
        # Maybe start by calculating d[lambda_i]/d[r_p] for {i,j,k} on each face p
        # Convenience
        mesh = self.mesh
        tri = self.tri



        for fh in tri.faces():
            mvh = mesh.vertex_handle(fh.idx())
            assert fh.idx() == mvh.idx()
        for vh in tri.vertices():
            mfh = mesh.face_handle(vh.idx())
            assert vh.idx() == mfh.idx()
        for mfh in mesh.faces():
            vh = tri.vertex_handle(mfh.idx())
            assert mfh.idx() == vh.idx()


        lvecprop = self.tri_lvecprop

        # drmudrp[mesh_vhid, mesh_fhid, :]
        drmudrp = {}

        for tri_fh in tri.faces():
            vh = mesh.vertex_handle(tri_fh.idx())
            #dlamqdrp {tri_vh.idx() : np(3,3) }
            dlamqdrp = {}
            drmudrp[vh.idx()] = {}
            # Get the corresponding face

            assert vh.idx() == tri_fh.idx()

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

            #print dlamqdrp[i] + dlamqdrp[j] + dlamqdrp[k]

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

        for boundary in self.boundaries:
            for j, vhid in enumerate(boundary):
                mvhid = self.btv_bmv[vhid]
                drmudrp[mvhid] = {}
                #print mvhid, vhid
                drmudrp[mvhid][vhid] = np.identity(3)


        def loop(trivh):
            trivhid = trivh.idx()
            boundary = tri.is_boundary(trivh)
            vhs = []
            if not boundary:
                fh = mesh.face_handle(trivhid)
            else:
                fh = mesh.face_handle(self.vh_mf[trivhid])

            for vh in mesh.fv(fh):
                vhidx = vh.idx()
                vhs.append(vhidx)
            #print vhs
            return vhs


        # dAdrmu[fhid][vhid]  
        dAdrmu = {}
        dPdrmu = {}
        # todo vertices
        #for fh in mesh.faces():
        for trivh in tri.vertices():

            trivhid = trivh.idx()
            dAdrmu[trivhid] = {}
            dPdrmu[trivhid] = {}

            fh = mesh.face_handle(trivhid)
            lp = loop(trivh)
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
                dAdrmu[trivhid][vhid] = np.cross(vplus, self.normal) - np.cross(vminus, self.normal)
                #print 'p, mu'
                #print trivhid, vhid
                
                #dPdrmu
                # get lengths
                vhpt = omvec(mesh.point(vh))
                lvm = vhpt - vminus
                lvp = vplus - vhpt
                dPdrmu[trivhid][vhid] = lvm/norm(lvm) - lvp/norm(lvp)


        fprop = VPropHandle()
        tri.add_property(fprop, 'force')
        self.fprop = fprop

        #self.imfprop = FPropHandle()
        #mesh.add_property(self.imfprop, 'imforce')

        #self.nnfprop = FPropHandle()
        #mesh.add_property(self.nnfprop, 'nnforce')

        # Could iterate over all the vertices and assign loops 
        # Put this in a separate loop for now since it is a bit special

        # setup helper functions
        # It's natural to use sets and the intersect() method for dealing with loops
        # Nearest Neighbour faces
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
        interloops = {}
        for vhi in tri.vertices():
            vhidx = vhi.idx()
            interloops[vhidx] = {}
            for vhj in nnfaces(vhi):
                vhjdx = vhj.idx()
                intset = set(loops[vhidx]).intersection(set(loops[vhjdx]))
                interloops[vhidx][vhjdx] = intset
        #print interloops

        #print drmudrp.keys()
#        print 'drmudrp'
        #for k, v in drmudrp.items():
            #print 'mu', k
            #print 'p', drmudrp[k].keys()
        #print 

        
        # It remains to do some complicated math over loops of nearest neighbours, etc..
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

            # tmp
            boundary = tri.is_boundary(trivh)
            #if boundary:
                #totalforce = np.zeros(3)
                #tri.set_property(fprop, trivh, totalforce)
                #continue

            # Immediate contribution
            farea_fac = -(kp/2.) * (area - ag *prefarea)  

            fprim_fac = -(gammap) * prim
            asum = np.zeros(3)
            psum = np.zeros(3)
            #print 'loop',loops[trivhid]
            #print 'p', trivhid

            for mu in loops[trivhid]:
                #print 'mu', mu
                    
                dAdrmu[trivhid][mu]
                drmudrp[mu][fhidx]
                
                ac = np.einsum('n,nm->m', dAdrmu[trivhid][mu], drmudrp[mu][fhidx])
                pc = np.einsum('n,nm->m', dPdrmu[trivhid][mu], drmudrp[mu][fhidx])
                asum += ac
                psum += pc

            farea = farea_fac * asum
            fprim = fprim_fac * psum

            # And the nearest neighbours contribution
            # Can put these blocks together by adding full loops to interloops
            # In the mean time just duplicate the code.
            area_nnsum = np.zeros(3)
            prim_nnsum = np.zeros(3)
            for vhnn, vidset in interloops[trivhid].items():
                
                #print 'looping over %d vertices' % len(vidset)
                nnvh = tri.vertex_handle(vhnn)
                #nnvh = vhnn.idx()

                kp = tri.property(self.kprop, nnvh)
                gammap = tri.property(self.gammaprop, nnvh)
                area = tri.property(self.areaprop, nnvh)
                prim = tri.property(self.primprop, nnvh)
                prefarea = tri.property(self.prefareaprop, nnvh)

                farea_fac = -(kp/2.) * (area - prefarea)  
                fprim_fac = -(gammap) * prim
                asum = np.zeros(3)
                psum = np.zeros(3)

                #print 'p', vhnn
                #print vidset
                for mu in vidset:
                    #print 'p', 'mu'
                    #print vhnn, mu
                    dAdrmu[vhnn][mu]
                    #print vhnn, fhidx
                    ac = np.einsum('n,nm->m', dAdrmu[vhnn][mu], drmudrp[mu][fhidx])
                    pc = np.einsum('n,nm->m', dPdrmu[vhnn][mu], drmudrp[mu][fhidx])
                    asum += ac
                    psum += pc
                area_nnsum += farea_fac * asum
                prim_nnsum += fprim_fac * psum

            imforce = farea + fprim
            nnforce = area_nnsum + prim_nnsum

            totalforce = imforce + nnforce

            tri.set_property(fprop, trivh, totalforce)
            #mesh.set_property(self.imfprop, fh, imforce)
            #mesh.set_property(self.nnfprop, fh, nnforce)

            # we got all the way to a value but now we need to go back and
            # fix up the r_q, lambda_q etc..

            # we can try visualising the dual, etc.

    def out_force_energy(self, outfile):
        # use OrderedDict
        # Still aren't explicitely handling the particle ids as we should be

        outfe = OrderedDict()
        ids, energy, fx, fy = [], [], [], []
        for vh in self.tri.vertices():
            if self.tri.is_boundary(vh) and self.boundary:
                continue
            ids.append(vh.idx())
            fh = self.mesh.face_handle(vh.idx())
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
            dump(outfe, fo)


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

    epidat = '/home/dan/cells/run/rpatch/epithelial_equilini.dat'

    import argparse
    parser = argparse.ArgumentParser()
 
    parser.add_argument("-i", "--input", type=str, default=epidat, help="Input dat file")
    parser.add_argument("-d", "--dir", type=str, default='/home/dan/tmp/', help="Output directory")
    #parser.add_argument("-o", "--output", type=str, default=epidat, help="Input dat file")

    args = parser.parse_args()
    epidat = args.input


    rdat = ReadData(epidat)
    PV = PVmesh.datbuild(rdat)

    # Handle K, and Gamma
    #nf = PV.mesh.n_faces()
    nf = PV.tri.n_vertices()
    # Could easily read these from a .conf file
    k = 1.
    gamma = 0.
    K = np.full(nf, k)
    Gamma = np.full(nf, gamma)
    PV.set_constants(K, Gamma)

    PV.calculate_energy()
    PV.calculate_forces()

    from writemesh import *
    import os.path as path


    outdir = args.dir


    mout = path.join(outdir, 'cellmesh.vtp')
    #print 'saving ', mout
    writemeshenergy(PV, mout)

    tout = path.join(outdir, 'trimesh.vtp')
    #print 'saving ', tout
    writetriforce(PV, tout)

    # Dump the force and energy
    fef = 'force_energy.dat'
    PV.out_force_energy(fef)



