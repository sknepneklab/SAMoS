# ################################################################
#  
#    Active Particles on Curved Spaces (APCS)
#    
#    Author: Silke Henkes
#   
#    ICSMB, Department of Physics
#    University of Aberdeen

#    Author: Rastko Sknepnek
#   
#    Division of Physics
#    School of Engineering, Physics and Mathematics
#    University of Dundee
#    
#    (c) 2013, 2014
#    
#    This program cannot be used, copied, or modified without
#    explicit permission of the author.
#  
# ################################################################

# Integrator code for batch processing of full data runs (incorporating parts of earlier analysis scripts)

# Data interfacing
from read_data import *
from read_param import *
# Pre-existing analysis scripts
from flat_patches_stresses_defects import *
#from glob import glob
    

# This is the structured data file hierarchy. Replace as appropriate (do not go the Yaouen way and fully automatize ...)
basefolder = '/home/silke/Documents/CurrentProjects/Rastko/Cells/Epithelial_model/morse_patch0125_b'
outfolder= '/home/silke/Documents/CurrentProjects/Rastko/analysis/'
L=100
sigma=1
nstep=500000
nsave=2000
nsnap=int(nstep/nsave)
skip=0
startvtk=0


#param = Param(basefolder)
#/home/silke/Documents/CurrentProjects/Rastko/Runs/RunsMarchJ1/data_v0_0.005/data_j_1_sphere/sphere_v0_0.005_j_1_0010000000.dat
files = sorted(glob(basefolder+'/plane_*.dat'))[skip:]
#files = sorted(glob(basefolder+'J_'+ J +'/sphere_*.dat'))[skip:]
defects=np.zeros((len(files),32))
ndefect=np.zeros((len(files),2))
u=0
for f in files:
    print f
    
    outname =basefolder+'/frame_data' + str(u-startvtk)+'.vtk'
    if u<startvtk:
        defects_n, defects_v,numdefect_n,numdefect_v=getDefects(f,float(r),sigma,outname,'polar',False,False)
    else:  
        outname_patches = basefolder+'/frame_patches' + str(u-startvtk)+'.vtk'
        defects_n, defects_v,numdefect_n,numdefect_v=getDefects(f,1.3,outname,outname_patches,'polar',True,True,True)
        outname_defects = basefolder+'/frame_defects' + str(u-startvtk)+'.vtk'
        print outname
        writeDefects(defects_n, defects_v,numdefect_n,numdefect_v,outname_defects)
    defects[u,0:4]=defects_n[0,:]
    defects[u,4:8]=defects_n[1,:]
    defects[u,8:12]=defects_n[2,:]
    defects[u,12:16]=defects_n[3,:]
    defects[u,16:20]=defects_v[0,:]
    defects[u,20:24]=defects_v[1,:]
    defects[u,24:28]=defects_v[2,:]
    defects[u,28:32]=defects_v[3,:]
    ndefect[u,0]=numdefect_n
    ndefect[u,1]=numdefect_v
    
    u+=1
    
outfile2=outfolder + 'defects_J_' + J + 'v0_'+ v0 +'_polar.dat'
np.savetxt(outfile2,np.concatenate((ndefect,defects),axis=1),fmt='%12.6g', header='ndefect (orientation, velocity) defects (orientation, velocity)')
