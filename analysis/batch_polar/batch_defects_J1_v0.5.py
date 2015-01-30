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
from polar_defects_analysis import *
#from glob import glob
    

# This is the structured data file hierarchy. Replace as appropriate (do not go the Yaouen way and fully automatize ...)
basefolder = '/home/silke/Documents/CurrentProjects/Rastko/Runs/RunsMarchJ1/'
outfolder= '/home/silke/Documents/CurrentProjects/Rastko/analysis/'
#JList=['10', '1', '0.1', '0.01']
#'0.005','0.01','0.02','0.05','0.1','0.2','0.5','1'
vList=['0.5','1']
JList=['1']
nu_r='0.002'
phi='1'
r=28.2094791
sigma=1
nstep=10000000
nsave=10000
nsnap=int(nstep/nsave)
skip=0
startvtk=0

for J in JList:
    for v0 in vList:
        #param = Param(basefolder)
        #/home/silke/Documents/CurrentProjects/Rastko/Runs/RunsMarchJ1/data_v0_0.005/data_j_1_sphere/sphere_v0_0.005_j_1_0010000000.dat
        files = sorted(glob(basefolder+'/data_v0_' + v0 + '/data_j_' + J +'_sphere/sphere_*.dat'))[skip:]
        #files = sorted(glob(basefolder+'J_'+ J +'/sphere_*.dat'))[skip:]
        defects=np.zeros((len(files),32))
        ndefect=np.zeros((len(files),2))
        u=0
        for f in files:
            print f
            
            outname =basefolder+'/data_v0_' + v0 + '/data_j_' + J +'_sphere/frame_data' + str(u-startvtk)+'.vtk'
            if u<startvtk:
                defects_n, defects_v,numdefect_n,numdefect_v=getDefects(f,float(r),sigma,outname,'polar',False,False)
            else:   
                defects_n, defects_v,numdefect_n,numdefect_v=getDefects(f,float(r),sigma,outname,'polar',False,True)
                outname = '.'.join((f).split('.')[:-1]) + '_defects.vtk'
                outname =basefolder+'/data_v0_' + v0 + '/data_j_' + J +'_sphere/frame_defects' + str(u-startvtk)+'.vtk'
                print outname
                writeDefects(defects_n, defects_v,numdefect_n,numdefect_v,outname)
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
