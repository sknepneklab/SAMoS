# *****************************************************************************
# *
# *  This Python script is a part of tha analysis of the data published in 
# *  the paper: "Universal motion patterns in confluent cell monolayers"
# *  by Silke Henkes, Kaja Kostanjevec, J. Martin Collinson, Rastko Sknepnek, 
# *  and Eric Bertin, Jounral name, vol, page (2019).
# *
# *  Please refer to the document Computational_summary.pdf for a detailed
# *  description of the tasks performed by this script.
# * 
# *****************************************************************************


import sys
import glob
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import pickle as pickle

from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Times']})
rc('text', usetex=True)

matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['axes.linewidth'] = 2
matplotlib.rcParams['xtick.major.size'] = 16
matplotlib.rcParams['xtick.minor.size'] = 0
matplotlib.rcParams['ytick.major.size'] = 16
matplotlib.rcParams['ytick.minor.size'] = 0
matplotlib.rcParams['font.size'] = 24.0
matplotlib.rcParams['legend.fontsize'] = 18.0

cdict = {'red':   [(0.0,  0.0, 0.75),
                   (0.25,  1.0, 0.75),
                   (0.45,  0.75, 0.0),
                   (1.0,  0.0, 0.0)],

         'green': [(0.0,  0.0, 0.0),
                   (0.25,  0.0, 0.5),
                   (0.5, 1.0, 1.0),
                   (0.8,  0.5, 0.0),
                   (1.0,  0.0, 0.0)],

         'blue':  [(0.0,  0.0, 0.0),
                   (0.5,  0.0, 0.0),
                   (0.7, 0.5, 1.0),
                   (1.0,  0.25, 0.0)]}


# number of experimental images to skip at the beginning
skip = 20  # i.e. 3h 20 minutes
plot = True

# time between snapshots in hours (10 minutes)
dt = 24/145.0
# Sixe of the field of view: in microns
length = 867.0
width = 662.0
# which is in pixels
# 1300 x 1000,
# This is the conversion factor between pixels and microns
dx0 = length/1300.0
dy0 = width/1000.0
# average as there should be no distortion present
dx = 0.5*(dx0+dy0)
dy = dx

# Field of view: 54x40 PIV arrows
nx = 54
ny = 40
N = nx*ny
# These values are the spacing in microns between arrows (different!!)
deltax = length/nx
deltay = length/ny

# Size of a cell: Typically about 1400 cells in view
ncell = 1400.0
cellarea = (nx*deltax)*(ny*deltay)/ncell
acell = np.sqrt(cellarea)
print acell

# Experiments that are part of the data shown in the paper
# expnamelist=['1','2','3','5','7']
# for demonstration puropses
expnamelist = ['7']
testmap = LinearSegmentedColormap('test', cdict, N=len(expnamelist))

plt.figure()
en = 0
nbin = 200
velbin = np.linspace(0, 10, nbin+1)
velcompbin = np.linspace(-10, 10, 2*nbin+1)

# This produces all of the points in q space where the Fourier transform should be taken
# Identical to the same function inside Glassy.py


def makeQrad(dx, nx, ny):
  # For Fourier transform
  qx = np.linspace(0, 2*np.pi/dx, nx)
  qy = np.linspace(0, 2*np.pi/dx, ny)
  dqx = qx[1]-qx[0]
  dqy = qy[1]-qy[0]
  dq = 0.5*(dqx+dqy)
  qmax2 = np.sqrt(qx[-1]**2+qy[-1]**2)
  nq2 = int(np.sqrt(nx**2+ny**2))
  qrad = np.linspace(0, qmax2, nq2)
  # do this silly counting once and for all
  binval = np.empty((nx, ny))
  for kx in range(nx):
    for ky in range(ny):
      qval = np.sqrt(qx[kx]**2+qy[ky]**2)
      binval[kx, ky] = round(qval/dq)
  ptsx = []
  ptsy = []
  # do the indexing arrays
  for l in range(nq2):
    pts0x = []
    pts0y = []
    for kx in range(nx):
      hmm = np.nonzero(binval[kx, :] == l)[0]
      for v in range(len(hmm)):
        pts0y.append(hmm[v])
        pts0x.append(kx)
      ptsx.append(pts0x)
      ptsy.append(pts0y)
  return qx, qy, qrad, ptsx, ptsy
  
 

dv = velbin[1] - velbin[0]
for expname in expnamelist:
  # Location of the PIV data
  directory = 'cleandata' + expname + '/' 
  files = sorted(glob.glob(directory + '*.txt'))[skip:]
  Nsnap = len(files)

  velhist = np.zeros((nbin,))
  velcomphist = np.zeros((2*nbin,))
  velav = np.zeros((Nsnap,))
  
  x = np.zeros((Nsnap,N))
  y = np.zeros((Nsnap,N))
  vx = np.zeros((Nsnap,N))
  vy = np.zeros((Nsnap,N))

  # read all of the PIV data of this experiment and convert into micron and hour units
  nf = 0
  for f in files[:Nsnap]:
    print f
    data = np.loadtxt(f)
    x[nf,:]  = dx*data[:,0]
    y[nf,:]  = dy*data[:,1]
    vx[nf,:] = dx*data[:,2]/dt
    vy[nf,:] = dx*data[:,3]/dt
    nf += 1
        
  # As detailed in the paper, the analysis is more robust if the velocities are normalised by their instantaneous spatial mean
  vxnorm = np.zeros((Nsnap,N))
  vynorm = np.zeros((Nsnap,N))
  for u in range(Nsnap):
      # mean velocity
      velav[u] = np.sqrt(np.sum(vx[u,:]**2+vy[u,:]**2)/N)
      vxnorm[u,:] = vx[u,:]/velav[u]
      vynorm[u,:] = vy[u,:]/velav[u]
      # scaled velocity distributions
      velhist0,vel_edges = np.histogram(np.sqrt(vxnorm[u,:]**2+vynorm[u,:]**2), bins=velbin,density=True)
      velhist += velhist0
      velhistx,vel_edges = np.histogram(vxnorm[u,:],bins=velcompbin,density=True)
      velhisty,vel_edges = np.histogram(vynorm[u,:],bins=velcompbin,density=True)
      velcomphist += 0.5*(velhistx+velhisty)
  velhist /= Nsnap
  velcomphist /= Nsnap
  velavtot = np.mean(velav)
  
  # Plots of velocity distributions
  if plot:
    plt.figure()
    plt.plot(np.arange(Nsnap),velav,'.-')
    plt.xlabel('frame')
    plt.ylabel('mean velocity (microns/hr)')
    plt.title('Mean velocity')
    plt.ylim(0,20)
  
  if plot:
    plt.figure()
    plt.loglog(velbin[:nbin]+dv/2.0,velhist,'.-',color=testmap(en),label='exp. ' +expname)
    plt.xlabel('velocity (micron / hour)')
    plt.ylabel('P(velocity)')
    
    plt.figure()
    plt.semilogy(velcompbin[:(2*nbin)]+dv/2.0,velcomphist,'.-',color=testmap(en),label='exp. ' +expname)
    plt.xlabel('velocity (micron / hour)')
    plt.ylabel('P(velocity components)')
  
  
  # Velocity autocorrelation function
  # normalise as best possible inside already
  velauto = np.zeros((Nsnap,))
  for u in range(Nsnap):
    smax = Nsnap-u
    velauto[u] = np.sum(np.sum((vxnorm[u:,:]*vxnorm[:smax,:]+vynorm[u:,:]*vynorm[:smax,:]),axis=1),axis=0)/(N*smax)
  tplot = np.array(range(Nsnap))*dt
  if plot:
    plt.figure()
    plt.semilogy(tplot,velauto,'.-',color=testmap(en),label='exp. ' +expname)
    plt.xlabel('time (hours)')
    plt.ylabel('velocity autocorrelation')
    plt.xlim(0,5)
      
      
  # Fourier transform of velocity
  # Borrowed from Glassy, in essence
  qx, qy, qrad, ptsx, ptsy = makeQrad(2*deltax,nx,ny)
  nq2 = int(np.sqrt(nx**2+ny**2))
  Sqrad = np.zeros((nq2-1,))
  Sqx = np.zeros((nx,ny))
  Sqy = np.zeros((nx,ny))
  for whichframe in range(0,Nsnap,1):
    print whichframe
    fourierval=np.zeros((nx,ny,2),dtype=complex)
    for kx in range(nx):
      for ky in range(ny):
        fourierval[kx,ky,0] = np.sum(np.exp(1j*(qx[kx]*x[whichframe,:] + qy[ky]*y[whichframe,:]))*vxnorm[whichframe,:])/N
        fourierval[kx,ky,1] = np.sum(np.exp(1j*(qx[kx]*x[whichframe,:] + qy[ky]*y[whichframe,:]))*vynorm[whichframe,:])/N
      # Sq = \vec{v_q}.\vec{v_-q}, assuming real and symmetric
      # = \vec{v_q}.\vec{v_q*} = v
      Sq = np.real(fourierval[:,:,0])**2+np.imag(fourierval[:,:,0])**2 + np.real(fourierval[:,:,1])**2+np.imag(fourierval[:,:,1])**2
      plotval_x = np.sqrt(np.real(fourierval[:,:,0])**2+np.imag(fourierval[:,:,0])**2)
      plotval_y = np.sqrt(np.real(fourierval[:,:,1])**2+np.imag(fourierval[:,:,1])**2)
      # Produce a radial averaging 
      Sqx += plotval_x
      Sqy += plotval_y
      Sqrad0 = np.zeros((nq2-1,))
      for l in range(nq2-1):
        Sqrad0[l] = np.mean(Sq[ptsx[l],ptsy[l]])
      Sqrad += Sqrad0
  Sqrad /= Nsnap
  Sqx /= Nsnap
  Sqy /= Nsnap
  if plot:
    plt.figure()
    plt.loglog(qrad[:(nq2-1)],Sqrad,'.-',color=testmap(en),label='exp. ' +expname)
    plt.xlabel('q (inverse microns)')
    plt.ylabel('Velocity Fourier')
    
    # plt.figure()
    # plt.pcolor(np.log(Sqx))
    
    # plt.figure()
    # plt.pcolor(np.log(Sqy))
    
  # And finally, our little reverse-engineered self-intermediate function
  # Idea here: Integrate the velocity (for a couple of steps only) to get the displacements of positions forwards in time
  tslide = 50
  tval = np.array(range(tslide))*dt
  # q is inverse cell diameter
  qval = 2*np.pi/acell
  # ok, need some plausible radial averaging. Be a little creative here ...
  sval = np.linspace(0,2*np.pi,36)
  SelfInt0 = np.zeros((tslide,),dtype=complex)
  for n in range(0,Nsnap,1):
    print n
    # those are the displacements already. go straigt into the self-intermediate
    xslide = np.zeros((tslide,N))
    yslide = np.zeros((tslide,N))
    for t in range(tslide):
      xslide[t,:] = np.sum(vxnorm[n:(n+t),:]*dt,axis=0)
      yslide[t,:] = np.sum(vynorm[n:(n+t),:]*dt,axis=0)
      sint0 = 0
      for s in sval:
        qx = qval*np.cos(s)
        qy = qval*np.sin(s)
        sint0 += np.sum(np.exp(1j*(qx*xslide[t,:]+qy*yslide[t,:])))/N
      SelfInt0[t] += sint0/(len(sval))
  SelfInt0 /= Nsnap
  SelfInt = np.sqrt(np.real(SelfInt0)**2 + np.imag(SelfInt0)**2)
  if plot:
    plt.figure()
    plt.semilogx(tval,SelfInt,'.-',color=testmap(en),label='exp. ' +expname)
    plt.xlabel('time (hours)')
    plt.ylabel('Self-Intermediate function')
      
    
  # save these results in a single pickle file
  outpickle = "PIV_results_" + "exp" + expname+ ".p"
  print outpickle
  data = {"exp": en,
          "length": length,
          "width": width,
          "dt": dt,
          "dx": dx,
          "deltax": deltax,
          "acell": acell,
          "Nsnap": Nsnap,
          "velbin": velbin,
          "velhist": velhist,
          "velcompbin": velcompbin,
          "velcomphist": velcomphist,
          "velav": velav,
          "tplot": tplot,
          "velauto": velauto,
          "qrad": qrad,
          "Sqrad": Sqrad,
          "tval": tval,
          "qval": qval,
          "SelfInt": SelfInt}
  pickle.dump(data,open(outpickle,'wb'))
  en += 1
    
plt.show()
