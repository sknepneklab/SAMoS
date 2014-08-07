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
from energy_profile_lib_v import *
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
#vList=['0.01','0.05','0.1','0.2','0.5','1']
#vList=['0.2','0.5','1']
vList=['1.0']
JList=['0.1']
nu_r='0.002'
phi='1'
sigma=1
nstep=20000000
nsave=10000
nsnap=int(nstep/nsave)
skip=int(nsnap/3)

nbin=180

for i in range(len(vList)):
	for j in range(len(JList)):
		print vList[i],JList[j]
		folder=basefolder+'data_v0_'+vList[i]+'/data_j_'+JList[j]+'_sphere/'
		print folder
		param = Param(folder)
		files = sorted(glob(folder+'*.dat'))[skip:]
		
		rho_profile =np.zeros((nbin,)) 
		vel_profile =np.zeros((nbin,)) 
		eng_profile = np.zeros((nbin,)) 
		press_profile = np.zeros((nbin,)) 
		s_tt_profile = np.zeros((nbin,)) 
		s_tp_profile = np.zeros((nbin,)) 
		s_pt_profile = np.zeros((nbin,)) 
		s_pp_profile = np.zeros((nbin,)) 
		alpha_profile = np.zeros((nbin,)) 
		alpha_v_profile = np.zeros((nbin,)) 
		axis = np.zeros((len(files),3)) 
		orderpar = np.zeros((len(files),3)) 
		iscount = np.zeros((nbin,)) 
		tot = 0
		for f in files:
			[theta_bin,rho_profile0,vel_profile0,eng_profile0,press_profile0,s_tt_profile0,s_tp_profile0,s_pt_profile0,s_pp_profile0,alpha_profile0,alpha_v_profile0,axis[tot,:],orderpar[tot,:]]=getProfiles(f,nbin,param.r,param.k,sigma)	
			isparticles=np.array([index for index,value in enumerate(rho_profile0) if (value >0)])
			iscount[isparticles]+=1
			rho_profile[isparticles]+=rho_profile0[isparticles]
			vel_profile[isparticles]+=vel_profile0[isparticles]
			eng_profile[isparticles]+=eng_profile0[isparticles]
			press_profile[isparticles]+=press_profile0[isparticles]
			s_tt_profile[isparticles]+=s_tt_profile0[isparticles]
			s_tp_profile[isparticles]+=s_tp_profile0[isparticles]
			s_pt_profile[isparticles]+=s_pt_profile0[isparticles]
			s_pp_profile[isparticles]+=s_pp_profile0[isparticles]
			alpha_profile[isparticles]+=alpha_profile0[isparticles]
			alpha_v_profile[isparticles]+=alpha_v_profile0[isparticles]
			tot +=1
		orderpar/=(param.v0*param.r)
		issomething=[index for index,value in enumerate(iscount) if (value >0)]
		rho_profile[issomething]/=iscount[issomething]
		vel_profile[issomething]/=iscount[issomething]
		eng_profile[issomething]/=iscount[issomething]
		press_profile[issomething]/=iscount[issomething]
		s_tt_profile[issomething]/=iscount[issomething]
		s_tp_profile[issomething]/=iscount[issomething]
		s_pt_profile[issomething]/=iscount[issomething]
		s_pp_profile[issomething]/=iscount[issomething]
		alpha_profile[issomething]/=iscount[issomething]
		alpha_v_profile[issomething]/=iscount[issomething]
		
		outfile=outfolder+'profilesV_v0' + vList[i] + '_j' + JList[j] + '.dat'
		outfile2=outfolder + 'axisV_v0' + vList[i] + '_j' + JList[j] + '.dat'
		np.savetxt(outfile, np.transpose(np.array([theta_bin,rho_profile,vel_profile,eng_profile,press_profile,s_tt_profile,s_tp_profile,s_pt_profile,s_pp_profile,alpha_profile,alpha_v_profile])),fmt='%12.6g', header='theta rho vel energy pressure s_tt s_tp s_pt s_pp alpha alpha_v')   # x,y,z equal sized 1D arrays
		np.savetxt(outfile2,np.concatenate((axis,orderpar),axis=1),fmt='%12.6g', header='axis orderpar')