# ################################################################
#  
#    Active Particles on Curved Spaces (APCS)
#    
#    Author: Silke Henkes
#   
#    ICSMB, Department of Physics
#    University of Aberdeen

#	 Author: Rastko Sknepnek
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
from nematic_analysis import *
#from glob import glob

	
# This is the structured data file hierarchy. Replace as appropriate (do not go the Yaouen way and fully automatize ...)
basefolder='/home/silke/Documents/CurrentProjects/Rastko/nematic/data/phi_1.0/'
#basefolder = '/home/silke/Documents/CurrentProjects/Rastko/nematic/data/J_1_0_v0_1_0/' 
#outfolder= '/home/silke/Documents/CurrentProjects/Rastko/nematic/data/J_1_0_v0_1_0/'
outfolder = '/home/silke/Documents/CurrentProjects/Rastko/nematic/data/phi_1.0/'

#Jval=['0.01','0.05','0.1','0.5','5.0','10.0']
Jval=['0.01']
sigma=1
r='16'
nstep=100100000
nsave=5000
nsnap=int(nstep/nsave)
#skip=835
skip=0
startvtk=17500

for J in Jval:
	#param = Param(basefolder)
	#/home/silke/Documents/CurrentProjects/Rastko/nematic/data/phi_0.75/J_0.01
	files = sorted(glob(basefolder+'J_'+ J +'/sphere_*.dat'))[skip:]
	defects=np.zeros((len(files),12))
	ndefect=np.zeros((len(files),1))
	u=0
	for f in files:
		print f
		
		outname =basefolder+'J_'+ J +'/frame_data' + str(u-startvtk)+'.vtk'
		if u<startvtk:
			defects0,ndefect0=getDefects(f,float(r),sigma,outname,False,False)
		else:	
			defects0,ndefect0=getDefects(f,float(r),sigma,outname,False,True)
			outname = '.'.join((f).split('.')[:-1]) + '_defects.vtk'
			outname =basefolder+'J_'+ J +'/frame_defects' + str(u-startvtk)+'.vtk'
			print outname
			writeDefects(defects0,ndefect0,outname)
		defects[u,0:3]=defects0[0,:]
		defects[u,3:6]=defects0[1,:]
		defects[u,6:9]=defects0[2,:]
		defects[u,9:12]=defects0[3,:]
		ndefect[u]=ndefect0
		
		u+=1
		
	outfile2=outfolder + 'defects_J_' + J + '_R_'+ r+ '.dat'
	np.savetxt(outfile2,np.concatenate((ndefect,defects),axis=1),fmt='%12.6g', header='ndefect defects')

