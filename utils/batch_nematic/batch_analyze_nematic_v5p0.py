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
basefolder='/home/silke/Documents/CurrentProjects/Rastko/nematic/data/'
#basefolder = '/home/silke/Documents/CurrentProjects/Rastko/nematic/data/J_1_0_v0_1_0/' 
#outfolder= '/home/silke/Documents/CurrentProjects/Rastko/nematic/data/J_1_0_v0_1_0/'
outfolder = '/home/silke/Documents/CurrentProjects/Rastko/nematic/data/'

#v0val=['0.2','1.0','5.0']
v0val=['5.0']
sigma=1
#rval = ['5.0','6.0','7.0','8.0','9.0','10.0','12.0','14.0','18.0','20.0','25.0','30.0','40.0']
rval = ['25.0','30.0']
nstep=10100000
nsave=5000
nsnap=int(nstep/nsave)
#skip=835
skip=0

for r in rval:
	for v0 in v0val:
		#param = Param(basefolder)
		files = sorted(glob(basefolder+'R_'+ r+ '/v0_' + v0 + '/sphere_*.dat'))[skip:]
		defects=np.zeros((len(files),12))
		ndefect=np.zeros((len(files),1))
		u=0
		for f in files:
			print f
			outname =outfolder +'R_'+ r+ '/v0_' + v0 + '/frame_data' + str(u)+'.vtk'
			defects0,ndefect0=getDefects(f,float(r),sigma,outname,False,False)
			defects[u,0:3]=defects0[0,:]
			defects[u,3:6]=defects0[1,:]
			defects[u,6:9]=defects0[2,:]
			defects[u,9:12]=defects0[3,:]
			ndefect[u]=ndefect0
			
			
			#outname = '.'.join((f).split('.')[:-1]) + '_defects.vtk'
			#outname =outfolder  +'R_'+ r+ '/v0_' + v0 + '/frame_defects' + str(u)+'.vtk'
			#print outname
			#writeDefects(defects0,ndefect0,outname)
			u+=1
			
		outfile2=outfolder + 'defects_v0_' + v0 + '_R_'+ r+ '.dat'
		np.savetxt(outfile2,np.concatenate((ndefect,defects),axis=1),fmt='%12.6g', header='ndefect defects')

