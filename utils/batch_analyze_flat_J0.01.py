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
from energy_profile_lib_flat import *
#from glob import glob
#from StressEnergy_plot import StressEnergy

#class GetAnalysis:
	#def __init__(self,folder,outfilename,skip):
		#self.outfilename = outfilename
		#self.folder = folder
		#self.skip = skip
		#self.param = Param(self.folder)

	#def CollectProfiles(self):
		#[theta_bin,en_prof, vel_prof, press_prof, tot_histo, rho_prof, north_pole]=EnProf.getProfiles(self.folder,self.skip,180,self.param.r,self.param.k)
	

# This is the structured data file hierarchy. Replace as appropriate (do not go the Yaouen way and fully automatize ...)
basefolder = '/home/silke/Documents/CurrentProjects/Rastko/Runs/'
outfolder= '/home/silke/Documents/CurrentProjects/Rastko/analysis/'
#JList=['10', '1', '0.1', '0.01']
vList=['0.005','0.01','0.02','0.05','0.1','0.2','0.5','1']
JList=['0.01']
nu_r='0.002'
phi='1'
sigma=1
nstep=10000000
nsave=10000
nsnap=int(nstep/nsave)
skip=int(nsnap/2)
skip=0


for i in range(len(vList)):
	for j in range(len(JList)):
		print vList[i],JList[j]
		folder=basefolder+'data_v0_'+vList[i]+'/data_j_'+JList[j]+'_plane/'
		print folder
		param = Param(folder)
		files = sorted(glob(folder+'*.dat'))[skip:]
		
		vel = np.zeros((len(files),)) 
		eng = np.zeros((len(files),)) 
		press = np.zeros((len(files),)) 
		s_tt = np.zeros((len(files),)) 
		s_tp = np.zeros((len(files),)) 
		s_pt = np.zeros((len(files),)) 
		s_pp = np.zeros((len(files),)) 
		alpha = np.zeros((len(files),)) 
		axis = np.zeros((len(files),3)) 
		orderpar = np.zeros((len(files),3)) 
		axisV = np.zeros((len(files),3)) 
		orderparV = np.zeros((len(files),3)) 
		tot = 0
		for f in files:
			[vel[tot],eng[tot],press[tot],s_tt[tot],s_tp[tot],s_pt[tot],s_pp[tot],alpha[tot],axis[tot,:],axisV[tot,:],orderpar[tot,:],orderparV[tot,:]]=getStats(f,param.k,sigma,param.lx)
			tot +=1
		
		outfile2=outfolder + 'flatstat_v0' + vList[i] + '_j' + JList[j] + '.dat'
		outmat=np.zeros((len(files),20))
		outmat[:,0]=vel
		outmat[:,1]=eng 
		outmat[:,2]=press 
		outmat[:,3]=s_tt 
		outmat[:,4]=s_tp
		outmat[:,5]=s_pt
		outmat[:,6]=s_pp
		outmat[:,7]=alpha
		outmat[:,8:11]=axis
		outmat[:,11:14]=axisV
		outmat[:,14:17]=orderpar
		outmat[:,17:20]=orderparV
				  
		np.savetxt(outfile2,outmat,fmt='%12.6g', header='velocity energy pressure snn snt stn stt alpha axis axisV orderpar orderparV')
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		