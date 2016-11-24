# helper methods for analysing a stack of pv objects and a transition list
import numpy as np

from numpy.linalg import eig
from numpy.linalg import norm 

# take a dictionary of matrices and decompose them into eigenvalues/vectors
def decompose(tdict):
    evalues, evectors = {}, {}
    for i, t in tdict.items():
        evalues[i], evectors[i] = eig(t)
    return evalues, evectors

# Calculate alignment between two sets of vectors
# arguments are equal sized numpy arrays
def palign(avecs, bvecs, nematic=False):
    sm = 0.
    ll = len(avecs)
    sm = []
    for i in range(ll):
        av = avecs[i]
        bv = bvecs[i]
        adotb = np.dot(av/norm(av), bv/norm(bv))
        if nematic: adotb = abs(adotb)  
        sm.append(adotb)
    return sm
    
def principle(evalue, evector):
    # construct the principle vector
    ll = len(evalue)
    pvalue = np.zeros((ll))
    lrow = len(evalue[0])
    pvector  = np.zeros((ll,lrow))
    for i, (ev, evec) in enumerate(zip(evalue, evector)):
        a, b = ev[:2]
        av, bv = evec[0,:], evec[1,:]
        swap = abs(a) >= abs(b)
        if swap:
            a = b
            av = bv
        pvalue[i] = a
        pvector[i] = av
    return pvalue, pvector


class TransAnalyse(object):
    def __init__(self, helist, tlist):
        self.helist = helist
        self.tlist = tlist
        if tlist:
            self.transitions = tlist.completed
        else:
            print 'Warning. No tlist object passed to TransAnalyse'

    def structure_alignment(self, outnum):
        # For constant shape index across cells expect that principle stress aligns closely with
        #  the texture tensor

        pv = self.helist[outnum]
        vir = pv.stresses['virial'].stress
        texture = pv.structure

        ll = len(vir)
        virvalues = np.zeros((ll, 2)); virvectors = np.zeros((ll, 2, 2))
        texvalues = np.zeros((ll, 2)); texvectors= np.zeros((ll,2, 2))

        tkeys = texture.keys()
        for i, vstress in enumerate(vir.values()):
            virvalues[i], virvectors[i] = eig(vstress[:2,:2])
            texvalues[i], texvectors[i] = eig(texture[tkeys[i]][:2,:2])

        pvir, pvv = principle(virvalues, virvectors)
        ptex, texvv = principle(texvalues, texvectors)

        aterms = palign(pvv, texvv, nematic=True)
        aterm =np.mean(aterms)
        print 'alignment'
        print aterm

    def transition_alignment(self):
        # T1 transitions are aligned along regions of high *local* stress

        for tone in self.transitions:
            print tone.get_orient()

### originally in Analyse_cells.py

# Object representing a T1 transition
class Tone(object):

    def __init__(self, vids, tstamp):
        # the current state of transition
        self.vids  = vids 
        # the flipped state
        self.xvids = None
        # The number of times flipped
        self.nx = 0
        # the step number when we started tracking
        self.tstamp = tstamp
        # time stamps of transitions
        self.flipstamp = []
        # track edge lengths
        # (tstamp,length)
        self.edl = []

    def update_edl(self, outnum, length):
        self.edl.append((outnum, edl))

    def flipped(self):
        return self.nx % 2 == 1

    def flip(self, vids, outnum):
        self.xvids = self.vids
        self.vids = vids
        self.nx +=  1
        self.flipstamp.append(outnum)

    # Get edge vector for the approach and the release of the transition
    def get_orient(helist):
        # These time stamps can be the same
        firstflip = self.flipstamp[0]
        lastflip = self.flipstamp[-1]
        # find the time step before the first transition
        timeline = helist.keys()
        prior =timeline[timeline.index(firstflip)-1]
        pvprior = helist[prior]
        pvafter = helist[lastflip]
        def get_edge(pv, vpair):
            tehid = pv.tri.vedge[tuple(vpair)]
            mehid = pv.tri.to_mesh_edge[tehid]
            mheh = pv.mesh.halfedge_handle(pv.mesh.edge_handle(mehid), 0)
            return pv.mesh.lvec[mheh.idx()]
        pr = get_edge(pvprior, self.xvids)
        af = get_edge(pvafter, self.vids)
        return pr, af


# invert a dictionary
def invert(dct):
    return dict([(v,k) for k,v in dct.items()])

#import collections.abc.MutableSequence as MutableSequence
import collections
import bisect
class Tlist(collections.MutableMapping):

    tcut = 0.1
    ttrack = 0.5

    def __init__(self, lst):
        self.lst = lst
        self.nt = 0
        self.completed = []

    def __getitem__(self, i):
        return self.lst[i]

    def __setitem__(self, i, v):
        self.lst[i] = v

    #def insert(self, i, v):
        #self.lst[i] = v

    def __iter__(self):
        return self

    def __len__(self):
        return len(self.lst)

    def __delitem__(self, i):
        del self.lst[i]

    def append(self, tone):
        assert self.keyinc not in self.lst
        self.lst[self.keyinc] = tone
        self.byvids[tone.vids] = self.keyinc
        self.keyinc += 1

    def remove(self, vpair):
        i = self.byvids[vpair]
        xvids = self.lst[i].xvids
        if xvids and xvids in self.byvids:
            del self.byvids[xvids]
        self.completed.append(self.lst[i])
        del self.lst[i]
        del self.byvids[vpair]

    # helper methods for the populate and update functions
    @staticmethod
    def _find_short(pv):
        lls = sorted(pv.mesh.lle.items(), key= lambda t: t[1])
        eids, llen = zip(*lls)
        #short = eids[:bisect.bisect_left(llen, Tlist.ttrack)]
        short = eids

        rshort = []
        dual_to_tri = invert(pv.tri.to_mesh_edge)
        #assert set(dual_to_tri.keys()) == set(pv.tri.to_mesh_edge.values())
        # Need to remove boundary trimesh edges
        for se in short:
            if se in dual_to_tri:
                rshort.append(se)
            else:
                #print 'boundary', se
                pass #boundary edge

        pv.mesh.dual_to_tri = dual_to_tri
        return rshort
    @staticmethod
    def _maketovpair(pv):
        #dual_to_tri = invert(pv.tri.to_mesh_edge)
        dual_to_tri = pv.mesh.dual_to_tri
        return dict([(eid, pv.tri.getvpair(teid)) for eid, teid in dual_to_tri.items()])

    # constructor which takes pv object
    @classmethod
    def populate(cls, pv, outnum):

        short = cls._find_short(pv)
        # dual to tri edges
        tovpair = cls._maketovpair(pv)
        lst = dict(enumerate([Tone(tovpair[eid], outnum) for eid in short]))
        byvids = dict([(tone.vids, i) for i, tone in lst.items()])
        #for vids, i in byvids.items():
            #tone.update_edl(outnum, tri.lle[eid])

        print 'started tracking {} edges'.format(len(lst))
        tt = cls(lst)

        tt.last_tovpair = tovpair
        tt.byvids = byvids
        tt.keyinc = len(lst)
        return tt

    def get_aux(self, tri, vpair):
        vid = tuple(vpair)
        heh = tri.halfedge_handle(tri.vedge[vid])
        oheh = tri.opposite_halfedge_handle(heh)
        ft =tri.face_handle(heh)
        fto = tri.face_handle(oheh)
        adjv = [v.idx() for v in tri.fv(ft)]
        adjvo = [v.idx() for v in tri.fv(fto)]
        tquad = set(adjv).union(set(adjvo))
        auxvids = frozenset(set(adjv).union(set(adjvo)) - vpair)
        return auxvids

    # now we need methods to update this list whenever edges change length or flip
    def update(self, pv, outnum):
        # setup
        pv.tri.vedgemap()
        tri = pv.tri; mesh = pv.mesh

        # assume that any edge that is created on this step is below track
        # Just check all edges every time to start with
        short = Tlist._find_short(pv)

        # Important! remeber that the mesh ids may not be consistent
        #  We are sure that the cell centre ids are consistent, that is all.

        tovpair = Tlist._maketovpair(pv)

        vidsls = [tovpair[eid] for eid in short]
        # look for new pairs of vertices
        newvids = set(tovpair.values()) - set(self.last_tovpair.values())
        # find which vid pairs that the new pairs correspond to (flips have occured)
        #print 'new', len(newvids)
        unaccounted = []
        for vpair in newvids:
            # Either auxvids is being tracked or 
            auxvids = self.get_aux(tri, vpair)
            if vpair in self.byvids:
                continue
            if auxvids not in self.byvids:
                unaccounted.append(vpair)
                continue
            tt = self.lst[self.byvids[auxvids]]
            # Keeping both references to the transition object
            self.byvids[vpair] = self.byvids[auxvids]
            tt.flip(vpair, outnum)
       
        # resolve edges which aren't short anymore
        for tone in self.lst.values():
            #tmp 
            try:
                hehid = tri.vedge[tuple(tone.vids)]
            except KeyError:
                print 'dropped a tone object'
                self.remove(tone.vids)
                continue
            el = tri.ll[hehid]
            if el > Tlist.tcut:
                if tone.xvids:
                    tri.ll[hehid]
                    # count this transition
                    self.nt += 1
                self.remove(tone.vids)

        #Find which of the entries in short are newvids and which arent
        shortpairs = [tovpair[eid] for eid in short]
        newshort = set(shortpairs) - set(newvids)
        #find the ones we aren't tracking yet
        untracked= newshort - set(self.byvids.keys())
        print 'transitions', self.nt
        print 'untracked', len(untracked)
        print 'unaccounted', len(unaccounted)
        print 'tracking', len(self.lst)
        # track any short edges
        for i, vpair in enumerate(untracked):
            self.append(Tone(vpair, outnum))

        # finishup
        self.last_tovpair = tovpair

 
        #if False: # turn off the more complicated transition resolver TODO
            #### deal with the unaccounted edges
            #ppairs = []
            #pairwith = newvids - set(unaccounted)
            ##print pairwith
            ##print unaccounted
            #for vpair in unaccounted:
                ## pair them up by finding the common vertices
                #va, vb = vpair
                #pairw = None
                #for pair in pairwith:
                    #if va in pair or vb in pair:
                        #pairw = pair
                        #break
                #if not pairw:
                    #print 'warning, was looking for a pair of edges, only found 1'
                    #print 'Need to track more egdes'
                    #print vpair
                #ppairs.append((vpair, pairw))

            ## Now we have the edge pairs, going to attempt to flip both of them and
            ## see which pair of flips recovers a pair of edges which are being tracked

            ## we want to find the flipped edge of the second pair given that the first
            ## edge is flipped
            
            ## going to assume that the connectivity of the triangulation is such that
            ##  there is a loop of 5 vertices with two central edges.
            #for vpair, pairw in ppairs:
                #vd, vb = vpair # this the edge we want to flip
                #ive, ivb = pairw # the edge we want to ignore (pre-flip)
                ## we can obtain all the vertices in the loop by tracing the three
                ##  faces associated with the two internal edges
                #hehv_id= tri.vedge[(vd,vb)]
                #hehi_id= tri.vedge[(ive, ivb)]
                #hehv =tri.halfedge_handle(hehv_id)
                #hehi =tri.halfedge_handle(hehi_id)
                ## find the third face
                #ovid = tri.opposite_halfedge_handle(hehv)
                #oiid = tri.opposite_halfedge_handle(hehi)
                #hhfaces = [hehv, hehi, ovid, oiid] # should be one duplicate face
                #loop = set()
                #for heh in hhfaces:
                    #fh = tri.face_handle(heh)
                    #verts = [v.idx() for v in tri.fv(fh)]
                    #loop = loop.union( set(verts) )
                #assert len(loop) == 5
                ## now identify the remaining edge by process of elimination
                #original = loop - set(pairw)
                #original = original- set(vpair)
                #original = frozenset(original)
                ## finally found the missing transition so count it and remove it
                #assert original in self.byvids
                #self.remove(original)
                #self.nt += 1


