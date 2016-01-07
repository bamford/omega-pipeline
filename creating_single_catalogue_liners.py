import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import math

'''
######################################################################
#### Single all catalogue hst based on total aperture measurements
######################################################################
filename = 'all_catalogue_hst_liners.txt'

f = np.loadtxt(filename,dtype='str')

ID_5,Prob_posiew_5,prob_agn5,prob_sf5,prob_retired5,prob_weak_agn5,\
cont_maxprob_5,cont_median_5,rel_error_cont_median_5,rel_error_cont_maxlik_5,\
z_maxprob_5,z_median_5,rel_error_z_median_5,\
fha_maxprob_5,fha_median_5,rel_error_fha_median_5,rel_error_fha_maxlik_5,lha_median_5,\
fnii_maxprob_5,fnii_median_5,rel_error_fnii_median_5,rel_error_fnii_maxlik_5,\
ewha_median_5,ewha_maxlik_5,error_ewha_median_5,error_ewha_maxlik_5,lineratio_median_5,lineratio_maxlik_5,error_lineratio_median_5,error_lineratio_maxlik_5,\
chisq_5,twolines_5,oneline_5,\
ID_total,Prob_posiew_total,prob_agn_total,prob_sf_total,prob_retiredtotal,prob_weak_agntotal,\
cont_maxprob_total,cont_median_total,rel_error_cont_median_total,rel_error_cont_maxlik_total,\
z_maxprob_total,z_median_total,rel_error_z_median_total,\
fha_maxprob_total,fha_median_total,rel_error_fha_median_total,rel_error_fha_maxlik_total,lha_median_total,\
fnii_maxprob_total,fnii_median_total,rel_error_fnii_median_total,rel_error_fnii_maxlik_total,\
ewha_median_total,ewha_maxlik_total,error_ewha_median_total,error_ewha_maxlik_total,lineratio_median_total,lineratio_maxlik_total,error_lineratio_median_total,error_lineratio_maxlik_total,\
chisq_total,twolines_total,oneline_total = f[:].T

formato = ('%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s \
            %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s \
            %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s \
            %s %s') 

f4 = open('single_all_catalogue_total_hst_liners.txt','w')
f4.write('%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s \
            %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s \
            %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s\
            %s \n'%('#',' ID_5', 'Prob_posiew_5', 'prob_agn5', 'prob_sf5','prob_retired5','prob_weak_agn5','cont_maxprob_5', 'cont_median_5', 'rel_error_cont_median_5', 'rel_error_cont_maxlik_5','z_maxprob_5', 'z_median_5', 'rel_error_z_median_5','fha_maxprob_5', 'fha_median_5', 'rel_error_fha_median_5', 'rel_error_fha_maxlik_5', 'lha_median_5','fnii_maxprob_5', 'fnii_median_5', 'rel_error_fnii_median_5', 'rel_error_fnii_maxlik_5','ewha_median_5', 'ewha_maxlik_5', 'error_ewha_median_5', 'error_ewha_maxlik_5', 'lineratio_median_5', 'lineratio_maxlik_5', 'error_lineratio_median_5', 'error_lineratio_maxlik_5','chisq_5', 'twolines_5', 'oneline_5','ID_total', 'Prob_posiew_total', 'prob_agn_total', 'prob_sf_total','prob_retiredtotal','prob_weak_agntotal','cont_maxprob_total', 'cont_median_total', 'rel_error_cont_median_total', 'rel_error_cont_maxlik_total','z_maxprob_total', 'z_median_total', 'rel_error_z_median_total','fha_maxprob_total', 'fha_median_total', 'rel_error_fha_median_total', 'rel_error_fha_maxlik_total', 'lha_median_total','fnii_maxprob_total', 'fnii_median_total', 'rel_error_fnii_median_total', 'rel_error_fnii_maxlik_total','ewha_median_total', 'ewha_maxlik_total', 'error_ewha_median_total', 'error_ewha_maxlik_total', 'lineratio_median_total', 'lineratio_maxlik_total', 'error_lineratio_median_total', 'error_lineratio_maxlik_total','chisq_total','twolines_total','oneline_total')) 
f4.close()





fha_median_total_nonstring = fha_median_total.astype(float)
rel_error_fha_median_total_nonstring = rel_error_fha_median_total.astype(float)

for ii in range(0,len(ID_total)):

    how_many2 = np.where((ID_total.astype(float) == ID_total[ii].astype(float)) & ((error_ewha_median_total.astype(float)/(ewha_median_total.astype(float)) < 0.15) & (Prob_posiew_total.astype(float)  > 99.7) & (oneline_total.astype(float) == 1)))[0]
  
    how_many = np.where(ID_total.astype(float) == ID_total[ii].astype(float))[0]

    if len(how_many) == 1:
        f_handle = file('single_all_catalogue_total_hst_liners.txt','a')
        np.savetxt(f_handle,f[ii].reshape(1,66),fmt=formato)
        f_handle.close()
        
    
    elif len(how_many) > 1:    
        
        if len(how_many2) == 0:
	    
            all_fluxes = fha_median_total_nonstring[how_many]*rel_error_fha_median_total_nonstring[how_many]

            if fha_median_total_nonstring[ii]*rel_error_fha_median_total_nonstring[ii] == all_fluxes.min():
                f_handle = file('single_all_catalogue_total_hst_liners.txt','a')
                np.savetxt(f_handle,f[ii].reshape(1,66),fmt=formato)
                f_handle.close()
        
        elif len(how_many2) > 0:
              
            all_fluxes = fha_median_total_nonstring[how_many2]*rel_error_fha_median_total_nonstring[how_many2]

            if fha_median_total_nonstring[ii]*rel_error_fha_median_total_nonstring[ii] == all_fluxes.min():
                f_handle = file('single_all_catalogue_total_hst_liners.txt','a')
                np.savetxt(f_handle,f[ii].reshape(1,66),fmt=formato)
                f_handle.close()
                
 
######################################################################
#### Single all catalogue hst based on PSF aperture measurements
######################################################################
filename = 'all_catalogue_hst_liners.txt'

f = np.loadtxt(filename,dtype='str')

ID_5,Prob_posiew_5,prob_agn5,prob_sf5,prob_retired5,prob_weak_agn5,\
cont_maxprob_5,cont_median_5,rel_error_cont_median_5,rel_error_cont_maxlik_5,\
z_maxprob_5,z_median_5,rel_error_z_median_5,\
fha_maxprob_5,fha_median_5,rel_error_fha_median_5,rel_error_fha_maxlik_5,lha_median_5,\
fnii_maxprob_5,fnii_median_5,rel_error_fnii_median_5,rel_error_fnii_maxlik_5,\
ewha_median_5,ewha_maxlik_5,error_ewha_median_5,error_ewha_maxlik_5,lineratio_median_5,lineratio_maxlik_5,error_lineratio_median_5,error_lineratio_maxlik_5,\
chisq_5,twolines_5,oneline_5,\
ID_total,Prob_posiew_total,prob_agn_total,prob_sf_total,prob_retiredtotal,prob_weak_agntotal,\
cont_maxprob_total,cont_median_total,rel_error_cont_median_total,rel_error_cont_maxlik_total,\
z_maxprob_total,z_median_total,rel_error_z_median_total,\
fha_maxprob_total,fha_median_total,rel_error_fha_median_total,rel_error_fha_maxlik_total,lha_median_total,\
fnii_maxprob_total,fnii_median_total,rel_error_fnii_median_total,rel_error_fnii_maxlik_total,\
ewha_median_total,ewha_maxlik_total,error_ewha_median_total,error_ewha_maxlik_total,lineratio_median_total,lineratio_maxlik_total,error_lineratio_median_total,error_lineratio_maxlik_total,\
chisq_total,twolines_total,oneline_total = f[:].T

formato = ('%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s \
            %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s \
            %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s \
            %s %s') 

f4 = open('single_all_catalogue_aper5_hst_liners.txt','w')
f4.write('%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s \
            %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s \
            %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s\
            %s \n'%('#',' ID_5', 'Prob_posiew_5', 'prob_agn5', 'prob_sf5','prob_retired5','prob_weak_agn5','cont_maxprob_5', 'cont_median_5', 'rel_error_cont_median_5', 'rel_error_cont_maxlik_5','z_maxprob_5', 'z_median_5', 'rel_error_z_median_5','fha_maxprob_5', 'fha_median_5', 'rel_error_fha_median_5', 'rel_error_fha_maxlik_5', 'lha_median_5','fnii_maxprob_5', 'fnii_median_5', 'rel_error_fnii_median_5', 'rel_error_fnii_maxlik_5','ewha_median_5', 'ewha_maxlik_5', 'error_ewha_median_5', 'error_ewha_maxlik_5', 'lineratio_median_5', 'lineratio_maxlik_5', 'error_lineratio_median_5', 'error_lineratio_maxlik_5','chisq_5', 'twolines_5', 'oneline_5','ID_total', 'Prob_posiew_total', 'prob_agn_total', 'prob_sf_total','prob_retiredtotal','prob_weak_agntotal','cont_maxprob_total', 'cont_median_total', 'rel_error_cont_median_total', 'rel_error_cont_maxlik_total','z_maxprob_total', 'z_median_total', 'rel_error_z_median_total','fha_maxprob_total', 'fha_median_total', 'rel_error_fha_median_total', 'rel_error_fha_maxlik_total', 'lha_median_total','fnii_maxprob_total', 'fnii_median_total', 'rel_error_fnii_median_total', 'rel_error_fnii_maxlik_total','ewha_median_total', 'ewha_maxlik_total', 'error_ewha_median_total', 'error_ewha_maxlik_total', 'lineratio_median_total', 'lineratio_maxlik_total', 'error_lineratio_median_total', 'error_lineratio_maxlik_total','chisq_total','twolines_total','oneline_total')) 
f4.close()



good = np.where(ID_5.astype(float) > 0)[0]
ID_5_good = ID_5[good].astype(float)
ewha_median_5_good = ewha_median_5[good].astype(float)
error_ewha_median_5_good = error_ewha_median_5[good].astype(float)
Prob_posiew_5_good = Prob_posiew_5[good].astype(float)
twolines_5_good = twolines_5[good].astype(float)

fha_median_aper5_nonstring = fha_median_5[good].astype(float)
rel_error_fha_median_aper5_nonstring = rel_error_fha_median_5[good].astype(float)

f_good = f[good]

for ii in range(0,len(ID_5_good)):

    how_many2 = np.where((ID_5_good == ID_5_good[ii]) & (error_ewha_median_5_good/(ewha_median_5_good) < 0.15) & (Prob_posiew_5_good  > 99.7) & (twolines_5_good == 1))[0]
  
    how_many = np.where(ID_5_good == ID_5_good[ii])[0]

    if len(how_many) == 1:
        f_handle = file('single_all_catalogue_aper5_hst_liners.txt','a')
        np.savetxt(f_handle,f_good[ii].reshape(1,66),fmt=formato)
        f_handle.close()
        
    
    elif len(how_many) > 1:    
        
        if len(how_many2) == 0:
	    
            all_fluxes = fha_median_aper5_nonstring[how_many]*rel_error_fha_median_aper5_nonstring[how_many]

            if fha_median_aper5_nonstring[ii]*rel_error_fha_median_aper5_nonstring[ii] == all_fluxes.min():
                f_handle = file('single_all_catalogue_aper5_hst_liners.txt','a')
                np.savetxt(f_handle,f_good[ii].reshape(1,66),fmt=formato)
                f_handle.close()
        
        elif len(how_many2) > 0:
              
            all_fluxes = fha_median_aper5_nonstring[how_many2]*rel_error_fha_median_aper5_nonstring[how_many2]

            if fha_median_aper5_nonstring[ii]*rel_error_fha_median_aper5_nonstring[ii] == all_fluxes.min():
                f_handle = file('single_all_catalogue_aper5_hst_liners.txt','a')
                np.savetxt(f_handle,f_good[ii].reshape(1,66),fmt=formato)
                f_handle.close()
'''
######################################################################
#### Single all catalogue master_hst based on total aperture measurements
######################################################################

filename = 'all_catalogue_master_hst_liners.txt'

f = np.loadtxt(filename,dtype='str')

field_number, ID_5,Prob_posiew_5,prob_agn5,prob_sf5,prob_retired5,prob_weak_agn5,\
cont_maxprob_5,cont_median_5,rel_error_cont_median_5,rel_error_cont_maxlik_5,\
z_maxprob_5,z_median_5,rel_error_z_median_5,\
fha_maxprob_5,fha_median_5,rel_error_fha_median_5,rel_error_fha_maxlik_5,lha_median_5,\
fnii_maxprob_5,fnii_median_5,rel_error_fnii_median_5,rel_error_fnii_maxlik_5,\
ewha_median_5,ewha_maxlik_5,error_ewha_median_5,error_ewha_maxlik_5,lineratio_median_5,lineratio_maxlik_5,error_lineratio_median_5,error_lineratio_maxlik_5,\
chisq_5,twolines_5,oneline_5,\
ID_total,Prob_posiew_total,prob_agn_total,prob_sf_total,prob_retiredtotal,prob_weak_agntotal,\
cont_maxprob_total,cont_median_total,rel_error_cont_median_total,rel_error_cont_maxlik_total,\
z_maxprob_total,z_median_total,rel_error_z_median_total,\
fha_maxprob_total,fha_median_total,rel_error_fha_median_total,rel_error_fha_maxlik_total,lha_median_total,\
fnii_maxprob_total,fnii_median_total,rel_error_fnii_median_total,rel_error_fnii_maxlik_total,\
ewha_median_total,ewha_maxlik_total,error_ewha_median_total,error_ewha_maxlik_total,lineratio_median_total,lineratio_maxlik_total,error_lineratio_median_total,error_lineratio_maxlik_total,\
chisq_total,twolines_total,oneline_total,\
ST_A_IMAGE,ra,dec,Rmag,combo_flag,stages_flag,MC_z,logmass,logmass_cl,flux24,tir,tuv,tir_cl,tuv_cl,\
sfr_det,sfr_lo,sfr_hi,sfr_det_cl,sfr_lo_cl,sfr_hi_cl,sed_type = f[:].T

formato = ('%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s \
                %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s \
                %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s\
		  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s') 

f4 = open('single_all_catalogue_total_master_hst_liners.txt','w')
f4.write('%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n'%('#','field_number', ' ID_5', 'Prob_posiew_5', 'prob_agn5', 'prob_sf5','prob_retired5','prob_weak_agn5','cont_maxprob_5', 'cont_median_5', 'rel_error_cont_median_5', 'rel_error_cont_maxlik_5','z_maxprob_5', 'z_median_5', 'rel_error_z_median_5','fha_maxprob_5', 'fha_median_5', 'rel_error_fha_median_5', 'rel_error_fha_maxlik_5', 'lha_median_5','fnii_maxprob_5', 'fnii_median_5', 'rel_error_fnii_median_5', 'rel_error_fnii_maxlik_5','ewha_median_5', 'ewha_maxlik_5', 'error_ewha_median_5', 'error_ewha_maxlik_5', 'lineratio_median_5', 'lineratio_maxlik_5', 'error_lineratio_median_5', 'error_lineratio_maxlik_5','chisq_5', 'twolines_5', 'oneline_5','ID_total', 'Prob_posiew_total', 'prob_agn_total', 'prob_sf_total','prob_retiredtotal','prob_weak_agntotal','cont_maxprob_total', 'cont_median_total', 'rel_error_cont_median_total', 'rel_error_cont_maxlik_total','z_maxprob_total', 'z_median_total', 'rel_error_z_median_total','fha_maxprob_total', 'fha_median_total', 'rel_error_fha_median_total', 'rel_error_fha_maxlik_total', 'lha_median_total','fnii_maxprob_total', 'fnii_median_total', 'rel_error_fnii_median_total', 'rel_error_fnii_maxlik_total','ewha_median_total', 'ewha_maxlik_total', 'error_ewha_median_total', 'error_ewha_maxlik_total', 'lineratio_median_total', 'lineratio_maxlik_total', 'error_lineratio_median_total', 'error_lineratio_maxlik_total','chisq_total','twolines_total','oneline_total','ST_A_IMAGE', 'ra', 'dec', 'Rmag', 'combo_flag', 'stages_flag', 'MC_z', 'logmass', 'logmass_cl', 'flux24', 'tir', 'tuv','tir_cl','tuv_cl','sfr_det', 'sfr_lo', 'sfr_hi', 'sfr_det_cl', 'sfr_lo_cl', 'sfr_hi_cl', 'sed_type')) 
f4.close()

galaxies_taken = []

fha_median_total_nonstring = fha_median_total.astype(float)
rel_error_fha_median_total_nonstring = rel_error_fha_median_total.astype(float)

for ii in range(0,len(ID_total)):
  
    how_many2 = np.where((np.isnan(MC_z.astype(float)) == False) & ((MC_z.astype(float) - z_median_total.astype(float)) <= (0.888*MC_z.astype(float) - 0.1358))  & (MC_z.astype(float) > 0.133) & (MC_z.astype(float) <= 0.195) & (ID_total.astype(float) == ID_total[ii].astype(float)) & ((error_ewha_median_total.astype(float)/(ewha_median_total.astype(float)) < 0.15) & (Prob_posiew_total.astype(float)  > 99.7) & (oneline_total.astype(float) == 1)))[0]
  
    how_many = np.where(ID_total.astype(float) == ID_total[ii].astype(float))[0]

    if len(how_many) == 1:
        f_handle = file('single_all_catalogue_total_master_hst_liners.txt','a')
        np.savetxt(f_handle,f[ii].reshape(1,88),fmt=formato)
        f_handle.close()
        
    elif len(how_many) > 1:    
        
        if len(how_many2) == 0:
	    
            all_fluxes = fha_median_total_nonstring[how_many]*rel_error_fha_median_total_nonstring[how_many]

            if fha_median_total_nonstring[ii]*rel_error_fha_median_total_nonstring[ii] == all_fluxes.min():
                f_handle = file('single_all_catalogue_total_master_hst_liners.txt','a')
                np.savetxt(f_handle,f[ii].reshape(1,88),fmt=formato)
                f_handle.close()
        
        elif len(how_many2) > 0:
              
            all_fluxes = fha_median_total_nonstring[how_many2]*rel_error_fha_median_total_nonstring[how_many2]

            if fha_median_total_nonstring[ii]*rel_error_fha_median_total_nonstring[ii] == all_fluxes.min():
                f_handle = file('single_all_catalogue_total_master_hst_liners.txt','a')
                np.savetxt(f_handle,f[ii].reshape(1,88),fmt=formato)
                f_handle.close()
                
######################################################################
#### Single all catalogue master_hst based on PSF aperture measurements
######################################################################

filename = 'all_catalogue_master_hst_liners.txt'

f = np.loadtxt(filename,dtype='str')

field_number, ID_5,Prob_posiew_5,prob_agn5,prob_sf5,prob_retired5,prob_weak_agn5,\
cont_maxprob_5,cont_median_5,rel_error_cont_median_5,rel_error_cont_maxlik_5,\
z_maxprob_5,z_median_5,rel_error_z_median_5,\
fha_maxprob_5,fha_median_5,rel_error_fha_median_5,rel_error_fha_maxlik_5,lha_median_5,\
fnii_maxprob_5,fnii_median_5,rel_error_fnii_median_5,rel_error_fnii_maxlik_5,\
ewha_median_5,ewha_maxlik_5,error_ewha_median_5,error_ewha_maxlik_5,lineratio_median_5,lineratio_maxlik_5,error_lineratio_median_5,error_lineratio_maxlik_5,\
chisq_5,twolines_5,oneline_5,\
ID_total,Prob_posiew_total,prob_agn_total,prob_sf_total,prob_retiredtotal,prob_weak_agntotal,\
cont_maxprob_total,cont_median_total,rel_error_cont_median_total,rel_error_cont_maxlik_total,\
z_maxprob_total,z_median_total,rel_error_z_median_total,\
fha_maxprob_total,fha_median_total,rel_error_fha_median_total,rel_error_fha_maxlik_total,lha_median_total,\
fnii_maxprob_total,fnii_median_total,rel_error_fnii_median_total,rel_error_fnii_maxlik_total,\
ewha_median_total,ewha_maxlik_total,error_ewha_median_total,error_ewha_maxlik_total,lineratio_median_total,lineratio_maxlik_total,error_lineratio_median_total,error_lineratio_maxlik_total,\
chisq_total,twolines_total,oneline_total,\
ST_A_IMAGE,ra,dec,Rmag,combo_flag,stages_flag,MC_z,logmass,logmass_cl,flux24,tir,tuv,tir_cl,tuv_cl,\
sfr_det,sfr_lo,sfr_hi,sfr_det_cl,sfr_lo_cl,sfr_hi_cl,sed_type = f[:].T

formato = ('%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s \
                %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s \
                %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s\
		  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s') 

f4 = open('single_all_catalogue_aper5_master_hst_liners.txt','w')
f4.write('%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n'%('#','field_number', ' ID_5', 'Prob_posiew_5', 'prob_agn5', 'prob_sf5','prob_retired5','prob_weak_agn5','cont_maxprob_5', 'cont_median_5', 'rel_error_cont_median_5', 'rel_error_cont_maxlik_5','z_maxprob_5', 'z_median_5', 'rel_error_z_median_5','fha_maxprob_5', 'fha_median_5', 'rel_error_fha_median_5', 'rel_error_fha_maxlik_5', 'lha_median_5','fnii_maxprob_5', 'fnii_median_5', 'rel_error_fnii_median_5', 'rel_error_fnii_maxlik_5','ewha_median_5', 'ewha_maxlik_5', 'error_ewha_median_5', 'error_ewha_maxlik_5', 'lineratio_median_5', 'lineratio_maxlik_5', 'error_lineratio_median_5', 'error_lineratio_maxlik_5','chisq_5', 'twolines_5', 'oneline_5','ID_total', 'Prob_posiew_total', 'prob_agn_total', 'prob_sf_total','prob_retiredtotal','prob_weak_agntotal','cont_maxprob_total', 'cont_median_total', 'rel_error_cont_median_total', 'rel_error_cont_maxlik_total','z_maxprob_total', 'z_median_total', 'rel_error_z_median_total','fha_maxprob_total', 'fha_median_total', 'rel_error_fha_median_total', 'rel_error_fha_maxlik_total', 'lha_median_total','fnii_maxprob_total', 'fnii_median_total', 'rel_error_fnii_median_total', 'rel_error_fnii_maxlik_total','ewha_median_total', 'ewha_maxlik_total', 'error_ewha_median_total', 'error_ewha_maxlik_total', 'lineratio_median_total', 'lineratio_maxlik_total', 'error_lineratio_median_total', 'error_lineratio_maxlik_total','chisq_total','twolines_total','oneline_total','ST_A_IMAGE', 'ra', 'dec', 'Rmag', 'combo_flag', 'stages_flag', 'MC_z', 'logmass', 'logmass_cl', 'flux24', 'tir', 'tuv','tir_cl','tuv_cl','sfr_det', 'sfr_lo', 'sfr_hi', 'sfr_det_cl', 'sfr_lo_cl', 'sfr_hi_cl', 'sed_type')) 
f4.close()

good = np.where(ID_5.astype(float) > 0)[0]
ID_5_good = ID_5[good].astype(float)
ewha_median_5_good = ewha_median_5[good].astype(float)
error_ewha_median_5_good = error_ewha_median_5[good].astype(float)
Prob_posiew_5_good = Prob_posiew_5[good].astype(float)
twolines_5_good = twolines_5[good].astype(float)
MC_z_good = MC_z[good].astype(float)
z_median_total_good = z_median_total[good].astype(float)

fha_median_aper5_nonstring = fha_median_5[good].astype(float)
rel_error_fha_median_aper5_nonstring = rel_error_fha_median_5[good].astype(float)

f_good = f[good]

for ii in range(0,len(ID_5_good)):
    how_many2 = np.where((np.isnan(MC_z_good) == False) & ((MC_z_good - z_median_total_good) <= (0.888*MC_z_good - 0.1358))  & (MC_z_good > 0.133) & (MC_z_good <= 0.195) & (ID_5_good == ID_5_good[ii]) & (error_ewha_median_5_good/(ewha_median_5_good) < 0.15) & (Prob_posiew_5_good  > 99.7) & (twolines_5_good == 1))[0]
  
    how_many = np.where(ID_5_good.astype(float) == ID_5_good[ii].astype(float))[0]

    if len(how_many) == 1:
        f_handle = file('single_all_catalogue_aper5_master_hst_liners.txt','a')
        np.savetxt(f_handle,f_good[ii].reshape(1,88),fmt=formato)
        f_handle.close()
        
    
    elif len(how_many) > 1:    
        
        if len(how_many2) == 0:
	    
            all_fluxes = fha_median_aper5_nonstring[how_many]*rel_error_fha_median_aper5_nonstring[how_many]

            if fha_median_aper5_nonstring[ii]*rel_error_fha_median_aper5_nonstring[ii] == all_fluxes.min():
                f_handle = file('single_all_catalogue_aper5_master_hst_liners.txt','a')
                np.savetxt(f_handle,f_good[ii].reshape(1,88),fmt=formato)
                f_handle.close()
        
        elif len(how_many2) > 0:
              
            all_fluxes = fha_median_aper5_nonstring[how_many2]*rel_error_fha_median_aper5_nonstring[how_many2]

            if fha_median_aper5_nonstring[ii]*rel_error_fha_median_aper5_nonstring[ii] == all_fluxes.min():
                f_handle = file('single_all_catalogue_aper5_master_hst_liners.txt','a')
                np.savetxt(f_handle,f_good[ii].reshape(1,88),fmt=formato)
                f_handle.close()
                
                

#TOTAL APERTURE
######################################################################
######################################################################
filename = 'all_catalogue_zcut_total_master_hst_liners.txt'

f = np.loadtxt(filename,dtype='str')

field_number, ID_5,Prob_posiew_5,prob_agn5,prob_sf5,prob_retired5,prob_weak_agn5,\
cont_maxprob_5,cont_median_5,rel_error_cont_median_5,rel_error_cont_maxlik_5,\
z_maxprob_5,z_median_5,rel_error_z_median_5,\
fha_maxprob_5,fha_median_5,rel_error_fha_median_5,rel_error_fha_maxlik_5,lha_median_5,\
fnii_maxprob_5,fnii_median_5,rel_error_fnii_median_5,rel_error_fnii_maxlik_5,\
ewha_median_5,ewha_maxlik_5,error_ewha_median_5,error_ewha_maxlik_5,lineratio_median_5,lineratio_maxlik_5,error_lineratio_median_5,error_lineratio_maxlik_5,\
chisq_5,twolines_5,oneline_5,\
ID_total,Prob_posiew_total,prob_agn_total,prob_sf_total,prob_retiredtotal,prob_weak_agntotal,\
cont_maxprob_total,cont_median_total,rel_error_cont_median_total,rel_error_cont_maxlik_total,\
z_maxprob_total,z_median_total,rel_error_z_median_total,\
fha_maxprob_total,fha_median_total,rel_error_fha_median_total,rel_error_fha_maxlik_total,lha_median_total,\
fnii_maxprob_total,fnii_median_total,rel_error_fnii_median_total,rel_error_fnii_maxlik_total,\
ewha_median_total,ewha_maxlik_total,error_ewha_median_total,error_ewha_maxlik_total,lineratio_median_total,lineratio_maxlik_total,error_lineratio_median_total,error_lineratio_maxlik_total,\
chisq_total,twolines_total,oneline_total,\
ST_A_IMAGE,ra,dec,Rmag,combo_flag,stages_flag,MC_z,logmass,logmass_cl,flux24,tir,tuv,tir_cl,tuv_cl,\
sfr_det,sfr_lo,sfr_hi,sfr_det_cl,sfr_lo_cl,sfr_hi_cl,sed_type = f[:].T

formato = ('%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s \
                %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s \
                %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s\
		  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s') 


f4 = open('single_all_catalogue_zcut_total_master_hst_liners.txt','w')
f4.write('%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n'%('#','field_number', ' ID_5', 'Prob_posiew_5', 'prob_agn5', 'prob_sf5','prob_retired5','prob_weak_agn5','cont_maxprob_5', 'cont_median_5', 'rel_error_cont_median_5', 'rel_error_cont_maxlik_5','z_maxprob_5', 'z_median_5', 'rel_error_z_median_5','fha_maxprob_5', 'fha_median_5', 'rel_error_fha_median_5', 'rel_error_fha_maxlik_5', 'lha_median_5','fnii_maxprob_5', 'fnii_median_5', 'rel_error_fnii_median_5', 'rel_error_fnii_maxlik_5','ewha_median_5', 'ewha_maxlik_5', 'error_ewha_median_5', 'error_ewha_maxlik_5', 'lineratio_median_5', 'lineratio_maxlik_5', 'error_lineratio_median_5', 'error_lineratio_maxlik_5','chisq_5', 'twolines_5', 'oneline_5','ID_total', 'Prob_posiew_total', 'prob_agn_total', 'prob_sf_total','prob_retiredtotal','prob_weak_agntotal','cont_maxprob_total', 'cont_median_total', 'rel_error_cont_median_total', 'rel_error_cont_maxlik_total','z_maxprob_total', 'z_median_total', 'rel_error_z_median_total','fha_maxprob_total', 'fha_median_total', 'rel_error_fha_median_total', 'rel_error_fha_maxlik_total', 'lha_median_total','fnii_maxprob_total', 'fnii_median_total', 'rel_error_fnii_median_total', 'rel_error_fnii_maxlik_total','ewha_median_total', 'ewha_maxlik_total', 'error_ewha_median_total', 'error_ewha_maxlik_total', 'lineratio_median_total', 'lineratio_maxlik_total', 'error_lineratio_median_total', 'error_lineratio_maxlik_total','chisq_total','twolines_total','oneline_total','ST_A_IMAGE', 'ra', 'dec', 'Rmag', 'combo_flag', 'stages_flag', 'MC_z', 'logmass', 'logmass_cl', 'flux24', 'tir', 'tuv','tir_cl','tuv_cl','sfr_det', 'sfr_lo', 'sfr_hi', 'sfr_det_cl', 'sfr_lo_cl', 'sfr_hi_cl', 'sed_type')) 
f4.close()



fha_median_total_nonstring = fha_median_total.astype(float)
rel_error_fha_median_total_nonstring = rel_error_fha_median_total.astype(float)

for ii in range(0,len(ID_total)):
  
    how_many2 = np.where((ID_total.astype(float) == ID_total[ii].astype(float)) & ((error_ewha_median_total.astype(float)/(ewha_median_total.astype(float)) < 0.15) & (Prob_posiew_total.astype(float)  > 99.7) & (oneline_total.astype(float) == 1)))[0]
    
    how_many = np.where(ID_total.astype(float) == ID_total[ii].astype(float))[0]

    if len(how_many) == 1:
        f_handle = file('single_all_catalogue_zcut_total_master_hst_liners.txt','a')
        np.savetxt(f_handle,f[ii].reshape(1,88),fmt=formato)
        f_handle.close()
        
    elif len(how_many) > 1:    
        
        if len(how_many2) == 0:
	    
            all_fluxes = fha_median_total_nonstring[how_many]*rel_error_fha_median_total_nonstring[how_many]

            if fha_median_total_nonstring[ii]*rel_error_fha_median_total_nonstring[ii] == all_fluxes.min():
                f_handle = file('single_all_catalogue_zcut_total_master_hst_liners.txt','a')
                np.savetxt(f_handle,f[ii].reshape(1,88),fmt=formato)
                f_handle.close()
        
        elif len(how_many2) > 0:
              
            all_fluxes = fha_median_total_nonstring[how_many2]*rel_error_fha_median_total_nonstring[how_many2]

            if fha_median_total_nonstring[ii]*rel_error_fha_median_total_nonstring[ii] == all_fluxes.min():
                f_handle = file('single_all_catalogue_zcut_total_master_hst_liners.txt','a')
                np.savetxt(f_handle,f[ii].reshape(1,88),fmt=formato)
                f_handle.close()
                
                
######################################################################
#PSF APERTURE
######################################################################
######################################################################
filename = 'all_catalogue_zcut_total_master_hst_liners.txt'

f = np.loadtxt(filename,dtype='str')

field_number, ID_5,Prob_posiew_5,prob_agn5,prob_sf5,prob_retired5,prob_weak_agn5,\
cont_maxprob_5,cont_median_5,rel_error_cont_median_5,rel_error_cont_maxlik_5,\
z_maxprob_5,z_median_5,rel_error_z_median_5,\
fha_maxprob_5,fha_median_5,rel_error_fha_median_5,rel_error_fha_maxlik_5,lha_median_5,\
fnii_maxprob_5,fnii_median_5,rel_error_fnii_median_5,rel_error_fnii_maxlik_5,\
ewha_median_5,ewha_maxlik_5,error_ewha_median_5,error_ewha_maxlik_5,lineratio_median_5,lineratio_maxlik_5,error_lineratio_median_5,error_lineratio_maxlik_5,\
chisq_5,twolines_5,oneline_5,\
ID_total,Prob_posiew_total,prob_agn_total,prob_sf_total,prob_retiredtotal,prob_weak_agntotal,\
cont_maxprob_total,cont_median_total,rel_error_cont_median_total,rel_error_cont_maxlik_total,\
z_maxprob_total,z_median_total,rel_error_z_median_total,\
fha_maxprob_total,fha_median_total,rel_error_fha_median_total,rel_error_fha_maxlik_total,lha_median_total,\
fnii_maxprob_total,fnii_median_total,rel_error_fnii_median_total,rel_error_fnii_maxlik_total,\
ewha_median_total,ewha_maxlik_total,error_ewha_median_total,error_ewha_maxlik_total,lineratio_median_total,lineratio_maxlik_total,error_lineratio_median_total,error_lineratio_maxlik_total,\
chisq_total,twolines_total,oneline_total,\
ST_A_IMAGE,ra,dec,Rmag,combo_flag,stages_flag,MC_z,logmass,logmass_cl,flux24,tir,tuv,tir_cl,tuv_cl,\
sfr_det,sfr_lo,sfr_hi,sfr_det_cl,sfr_lo_cl,sfr_hi_cl,sed_type = f[:].T

formato = ('%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s \
                %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s \
                %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s\
		  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s') 


f4 = open('single_all_catalogue_zcut_aper5_master_hst_liners.txt','w')
f4.write('%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n'%('#','field_number', ' ID_5', 'Prob_posiew_5', 'prob_agn5', 'prob_sf5','prob_retired5','prob_weak_agn5','cont_maxprob_5', 'cont_median_5', 'rel_error_cont_median_5', 'rel_error_cont_maxlik_5','z_maxprob_5', 'z_median_5', 'rel_error_z_median_5','fha_maxprob_5', 'fha_median_5', 'rel_error_fha_median_5', 'rel_error_fha_maxlik_5', 'lha_median_5','fnii_maxprob_5', 'fnii_median_5', 'rel_error_fnii_median_5', 'rel_error_fnii_maxlik_5','ewha_median_5', 'ewha_maxlik_5', 'error_ewha_median_5', 'error_ewha_maxlik_5', 'lineratio_median_5', 'lineratio_maxlik_5', 'error_lineratio_median_5', 'error_lineratio_maxlik_5','chisq_5', 'twolines_5', 'oneline_5','ID_total', 'Prob_posiew_total', 'prob_agn_total', 'prob_sf_total','prob_retiredtotal','prob_weak_agntotal','cont_maxprob_total', 'cont_median_total', 'rel_error_cont_median_total', 'rel_error_cont_maxlik_total','z_maxprob_total', 'z_median_total', 'rel_error_z_median_total','fha_maxprob_total', 'fha_median_total', 'rel_error_fha_median_total', 'rel_error_fha_maxlik_total', 'lha_median_total','fnii_maxprob_total', 'fnii_median_total', 'rel_error_fnii_median_total', 'rel_error_fnii_maxlik_total','ewha_median_total', 'ewha_maxlik_total', 'error_ewha_median_total', 'error_ewha_maxlik_total', 'lineratio_median_total', 'lineratio_maxlik_total', 'error_lineratio_median_total', 'error_lineratio_maxlik_total','chisq_total','twolines_total','oneline_total','ST_A_IMAGE', 'ra', 'dec', 'Rmag', 'combo_flag', 'stages_flag', 'MC_z', 'logmass', 'logmass_cl', 'flux24', 'tir', 'tuv','tir_cl','tuv_cl','sfr_det', 'sfr_lo', 'sfr_hi', 'sfr_det_cl', 'sfr_lo_cl', 'sfr_hi_cl', 'sed_type')) 
f4.close()



good = np.where(ID_5.astype(float) > 0)[0]
ID_5_good = ID_5[good].astype(float)
ewha_median_5_good = ewha_median_5[good].astype(float)
error_ewha_median_5_good = error_ewha_median_5[good].astype(float)
Prob_posiew_5_good = Prob_posiew_5[good].astype(float)
twolines_5_good = twolines_5[good].astype(float)

fha_median_aper5_nonstring = fha_median_5[good].astype(float)
rel_error_fha_median_aper5_nonstring = rel_error_fha_median_5[good].astype(float)

f_good = f[good]

for ii in range(0,len(ID_5_good)):
  
    how_many2 = np.where((ID_5_good == ID_5_good[ii]) & (error_ewha_median_5_good/(ewha_median_5_good) < 0.15) & (Prob_posiew_5_good  > 99.7) & (twolines_5_good == 1))[0]
    
    how_many = np.where(ID_5_good == ID_5_good[ii])[0]

    if len(how_many) == 1:
        f_handle = file('single_all_catalogue_zcut_aper5_master_hst_liners.txt','a')
        np.savetxt(f_handle,f_good[ii].reshape(1,88),fmt=formato)
        f_handle.close()
        
    elif len(how_many) > 1:    
        
        if len(how_many2) == 0:
	    
            all_fluxes = fha_median_aper5_nonstring[how_many]*rel_error_fha_median_aper5_nonstring[how_many]

            if fha_median_aper5_nonstring[ii]*rel_error_fha_median_aper5_nonstring[ii] == all_fluxes.min():
                f_handle = file('single_all_catalogue_zcut_aper5_master_hst_liners.txt','a')
                np.savetxt(f_handle,f_good[ii].reshape(1,88),fmt=formato)
                f_handle.close()
        
        elif len(how_many2) > 0:
              
            all_fluxes = fha_median_aper5_nonstring[how_many2]*rel_error_fha_median_aper5_nonstring[how_many2]

            if fha_median_aper5_nonstring[ii]*rel_error_fha_median_aper5_nonstring[ii] == all_fluxes.min():
                f_handle = file('single_all_catalogue_zcut_aper5_master_hst_liners.txt','a')
                np.savetxt(f_handle,f_good[ii].reshape(1,88),fmt=formato)
                f_handle.close()
                
                
######################################################################



################################################################################   
filename = 'all_catalogue_cut_two_lines_aper5_master_hst_liners.txt'

f = np.loadtxt(filename,dtype='str')

field_number, ID_5,Prob_posiew_5,prob_agn5,prob_sf5,prob_retired5,prob_weak_agn5,\
cont_maxprob_5,cont_median_5,rel_error_cont_median_5,rel_error_cont_maxlik_5,\
z_maxprob_5,z_median_5,rel_error_z_median_5,\
fha_maxprob_5,fha_median_5,rel_error_fha_median_5,rel_error_fha_maxlik_5,lha_median_5,\
fnii_maxprob_5,fnii_median_5,rel_error_fnii_median_5,rel_error_fnii_maxlik_5,\
ewha_median_5,ewha_maxlik_5,error_ewha_median_5,error_ewha_maxlik_5,lineratio_median_5,lineratio_maxlik_5,error_lineratio_median_5,error_lineratio_maxlik_5,\
chisq_5,twolines_5,oneline_5,\
ID_total,Prob_posiew_total,prob_agn_total,prob_sf_total,prob_retiredtotal,prob_weak_agntotal,\
cont_maxprob_total,cont_median_total,rel_error_cont_median_total,rel_error_cont_maxlik_total,\
z_maxprob_total,z_median_total,rel_error_z_median_total,\
fha_maxprob_total,fha_median_total,rel_error_fha_median_total,rel_error_fha_maxlik_total,lha_median_total,\
fnii_maxprob_total,fnii_median_total,rel_error_fnii_median_total,rel_error_fnii_maxlik_total,\
ewha_median_total,ewha_maxlik_total,error_ewha_median_total,error_ewha_maxlik_total,lineratio_median_total,lineratio_maxlik_total,error_lineratio_median_total,error_lineratio_maxlik_total,\
chisq_total,twolines_total,oneline_total,\
ST_A_IMAGE,ra,dec,Rmag,combo_flag,stages_flag,MC_z,logmass,logmass_cl,flux24,tir,tuv,tir_cl,tuv_cl,\
sfr_det,sfr_lo,sfr_hi,sfr_det_cl,sfr_lo_cl,sfr_hi_cl,sed_type = f[:].T

formato = ('%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s \
                %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s \
                %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s\
		  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s') 

good = np.where(ID_5.astype(float) > 0)[0]
ID_5_good = ID_5[good].astype(float)
ewha_median_5_good = ewha_median_5[good].astype(float)
error_ewha_median_5_good = error_ewha_median_5[good].astype(float)
Prob_posiew_5_good = Prob_posiew_5[good].astype(float)
twolines_5_good = twolines_5[good].astype(float)
MC_z_good = MC_z[good].astype(float)
z_median_total_good = z_median_total[good].astype(float)

fha_median_aper5_nonstring = fha_median_5[good].astype(float)
rel_error_fha_median_aper5_nonstring = rel_error_fha_median_5[good].astype(float)

f_good = f[good]

f4 = open('single_all_catalogue_cut_two_lines_aper5_master_hst_liners.txt','w')
f4.write('%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n'%('#','field_number', ' ID_5', 'Prob_posiew_5', 'prob_agn5', 'prob_sf5','prob_retired5','prob_weak_agn5','cont_maxprob_5', 'cont_median_5', 'rel_error_cont_median_5', 'rel_error_cont_maxlik_5','z_maxprob_5', 'z_median_5', 'rel_error_z_median_5','fha_maxprob_5', 'fha_median_5', 'rel_error_fha_median_5', 'rel_error_fha_maxlik_5', 'lha_median_5','fnii_maxprob_5', 'fnii_median_5', 'rel_error_fnii_median_5', 'rel_error_fnii_maxlik_5','ewha_median_5', 'ewha_maxlik_5', 'error_ewha_median_5', 'error_ewha_maxlik_5', 'lineratio_median_5', 'lineratio_maxlik_5', 'error_lineratio_median_5', 'error_lineratio_maxlik_5','chisq_5', 'twolines_5', 'oneline_5','ID_total', 'Prob_posiew_total', 'prob_agn_total', 'prob_sf_total','prob_retiredtotal','prob_weak_agntotal','cont_maxprob_total', 'cont_median_total', 'rel_error_cont_median_total', 'rel_error_cont_maxlik_total','z_maxprob_total', 'z_median_total', 'rel_error_z_median_total','fha_maxprob_total', 'fha_median_total', 'rel_error_fha_median_total', 'rel_error_fha_maxlik_total', 'lha_median_total','fnii_maxprob_total', 'fnii_median_total', 'rel_error_fnii_median_total', 'rel_error_fnii_maxlik_total','ewha_median_total', 'ewha_maxlik_total', 'error_ewha_median_total', 'error_ewha_maxlik_total', 'lineratio_median_total', 'lineratio_maxlik_total', 'error_lineratio_median_total', 'error_lineratio_maxlik_total','chisq_total','twolines_total','oneline_total','ST_A_IMAGE', 'ra', 'dec', 'Rmag', 'combo_flag', 'stages_flag', 'MC_z', 'logmass', 'logmass_cl', 'flux24', 'tir', 'tuv','tir_cl','tuv_cl','sfr_det', 'sfr_lo', 'sfr_hi', 'sfr_det_cl', 'sfr_lo_cl', 'sfr_hi_cl', 'sed_type')) 
f4.close()

for ii in range(0,len(ID_5_good)):
    how_many2 = np.where((np.isnan(MC_z_good) == False) & ((MC_z_good - z_median_total_good) <= (0.888*MC_z_good - 0.1358))  & (MC_z_good > 0.133) & (MC_z_good <= 0.195) & (ID_5_good == ID_5_good[ii]) & (error_ewha_median_5_good/(ewha_median_5_good) < 0.15) & (Prob_posiew_5_good  > 99.7) & (twolines_5_good == 1))[0]
  
    how_many = np.where(ID_5_good.astype(float) == ID_5_good[ii].astype(float))[0]

    if len(how_many) == 1:
        f_handle = file('single_all_catalogue_cut_two_lines_aper5_master_hst_liners.txt','a')
        np.savetxt(f_handle,f_good[ii].reshape(1,88),fmt=formato)
        f_handle.close()
        
    elif len(how_many) > 1:    
        
        if len(how_many2) == 0:
	    
            all_fluxes = fha_median_aper5_nonstring[how_many]*rel_error_fha_median_aper5_nonstring[how_many]

            if fha_median_aper5_nonstring[ii]*rel_error_fha_median_aper5_nonstring[ii] == all_fluxes.min():
                f_handle = file('single_all_catalogue_cut_two_lines_aper5_master_hst_liners.txt','a')
                np.savetxt(f_handle,f_good[ii].reshape(1,88),fmt=formato)
                f_handle.close()
        
        elif len(how_many2) > 0:
              
            all_fluxes = fha_median_aper5_nonstring[how_many2]*rel_error_fha_median_aper5_nonstring[how_many2]

            if fha_median_aper5_nonstring[ii]*rel_error_fha_median_aper5_nonstring[ii] == all_fluxes.min():
                f_handle = file('single_all_catalogue_cut_two_lines_aper5_master_hst_liners.txt','a')
                np.savetxt(f_handle,f_good[ii].reshape(1,88),fmt=formato)
                f_handle.close()



################################################################################  
filename = 'all_catalogue_cut_one_line_total_master_hst_liners.txt'

f = np.loadtxt(filename,dtype='str')

field_number, ID_5,Prob_posiew_5,prob_agn5,prob_sf5,prob_retired5,prob_weak_agn5,\
cont_maxprob_5,cont_median_5,rel_error_cont_median_5,rel_error_cont_maxlik_5,\
z_maxprob_5,z_median_5,rel_error_z_median_5,\
fha_maxprob_5,fha_median_5,rel_error_fha_median_5,rel_error_fha_maxlik_5,lha_median_5,\
fnii_maxprob_5,fnii_median_5,rel_error_fnii_median_5,rel_error_fnii_maxlik_5,\
ewha_median_5,ewha_maxlik_5,error_ewha_median_5,error_ewha_maxlik_5,lineratio_median_5,lineratio_maxlik_5,error_lineratio_median_5,error_lineratio_maxlik_5,\
chisq_5,twolines_5,oneline_5,\
ID_total,Prob_posiew_total,prob_agn_total,prob_sf_total,prob_retiredtotal,prob_weak_agntotal,\
cont_maxprob_total,cont_median_total,rel_error_cont_median_total,rel_error_cont_maxlik_total,\
z_maxprob_total,z_median_total,rel_error_z_median_total,\
fha_maxprob_total,fha_median_total,rel_error_fha_median_total,rel_error_fha_maxlik_total,lha_median_total,\
fnii_maxprob_total,fnii_median_total,rel_error_fnii_median_total,rel_error_fnii_maxlik_total,\
ewha_median_total,ewha_maxlik_total,error_ewha_median_total,error_ewha_maxlik_total,lineratio_median_total,lineratio_maxlik_total,error_lineratio_median_total,error_lineratio_maxlik_total,\
chisq_total,twolines_total,oneline_total,\
ST_A_IMAGE,ra,dec,Rmag,combo_flag,stages_flag,MC_z,logmass,logmass_cl,flux24,tir,tuv,tir_cl,tuv_cl,\
sfr_det,sfr_lo,sfr_hi,sfr_det_cl,sfr_lo_cl,sfr_hi_cl,sed_type = f[:].T

formato = ('%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s \
                %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s \
                %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s\
		  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s') 


f4 = open('single_all_catalogue_cut_one_line_total_aper_master_hst_liners.txt','w')
f4.write('%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n'%('#','field_number', ' ID_5', 'Prob_posiew_5', 'prob_agn5', 'prob_sf5','prob_retired5','prob_weak_agn5','cont_maxprob_5', 'cont_median_5', 'rel_error_cont_median_5', 'rel_error_cont_maxlik_5','z_maxprob_5', 'z_median_5', 'rel_error_z_median_5','fha_maxprob_5', 'fha_median_5', 'rel_error_fha_median_5', 'rel_error_fha_maxlik_5', 'lha_median_5','fnii_maxprob_5', 'fnii_median_5', 'rel_error_fnii_median_5', 'rel_error_fnii_maxlik_5','ewha_median_5', 'ewha_maxlik_5', 'error_ewha_median_5', 'error_ewha_maxlik_5', 'lineratio_median_5', 'lineratio_maxlik_5', 'error_lineratio_median_5', 'error_lineratio_maxlik_5','chisq_5', 'twolines_5', 'oneline_5','ID_total', 'Prob_posiew_total', 'prob_agn_total', 'prob_sf_total','prob_retiredtotal','prob_weak_agntotal','cont_maxprob_total', 'cont_median_total', 'rel_error_cont_median_total', 'rel_error_cont_maxlik_total','z_maxprob_total', 'z_median_total', 'rel_error_z_median_total','fha_maxprob_total', 'fha_median_total', 'rel_error_fha_median_total', 'rel_error_fha_maxlik_total', 'lha_median_total','fnii_maxprob_total', 'fnii_median_total', 'rel_error_fnii_median_total', 'rel_error_fnii_maxlik_total','ewha_median_total', 'ewha_maxlik_total', 'error_ewha_median_total', 'error_ewha_maxlik_total', 'lineratio_median_total', 'lineratio_maxlik_total', 'error_lineratio_median_total', 'error_lineratio_maxlik_total','chisq_total','twolines_total','oneline_total','ST_A_IMAGE', 'ra', 'dec', 'Rmag', 'combo_flag', 'stages_flag', 'MC_z', 'logmass', 'logmass_cl', 'flux24', 'tir', 'tuv','tir_cl','tuv_cl','sfr_det', 'sfr_lo', 'sfr_hi', 'sfr_det_cl', 'sfr_lo_cl', 'sfr_hi_cl', 'sed_type')) 
f4.close()


galaxies_taken = []

fha_median_total_nonstring = fha_median_total.astype(float)
rel_error_fha_median_total_nonstring = rel_error_fha_median_total.astype(float)

for ii in range(0,len(ID_total)):

    how_many2 = np.where((np.isnan(MC_z.astype(float)) == False) & ((MC_z.astype(float) - z_median_total.astype(float)) <= (0.888*MC_z.astype(float) - 0.1358))  & (MC_z.astype(float) > 0.133) & (MC_z.astype(float) <= 0.195) & (ID_total.astype(float) == ID_total[ii].astype(float)) & ((error_ewha_median_total.astype(float)/(ewha_median_total.astype(float)) < 0.15) & (Prob_posiew_total.astype(float)  > 99.7) & (oneline_total.astype(float) == 1)))[0]
  
    how_many = np.where(ID_total.astype(float) == ID_total[ii].astype(float))[0]

    if len(how_many) == 1:
        f_handle = file('single_all_catalogue_cut_one_line_total_aper_master_hst_liners.txt','a')
        np.savetxt(f_handle,f[ii].reshape(1,88),fmt=formato)
        f_handle.close()
        
    
    elif len(how_many) > 1:    
        
        if len(how_many2) == 0:
	    
            all_fluxes = fha_median_total_nonstring[how_many]*rel_error_fha_median_total_nonstring[how_many]

            if fha_median_total_nonstring[ii]*rel_error_fha_median_total_nonstring[ii] == all_fluxes.min():
                f_handle = file('single_all_catalogue_cut_one_line_total_aper_master_hst_liners.txt','a')
                np.savetxt(f_handle,f[ii].reshape(1,88),fmt=formato)
                f_handle.close()
        
        elif len(how_many2) > 0:
              
            all_fluxes = fha_median_total_nonstring[how_many2]*rel_error_fha_median_total_nonstring[how_many2]

            if fha_median_total_nonstring[ii]*rel_error_fha_median_total_nonstring[ii] == all_fluxes.min():
                f_handle = file('single_all_catalogue_cut_one_line_total_aper_master_hst_liners.txt','a')
                np.savetxt(f_handle,f[ii].reshape(1,88),fmt=formato)
                f_handle.close()
'''
################################################################################   
#### NO Z CUT ################################################################## 
################################################################################   
filename = 'all_catalogue_cut_two_lines_aper5_master_hst_nozcut_liners.txt'

f = np.loadtxt(filename,dtype='str')

field_number, ID_5,Prob_posiew_5,prob_agn5,prob_sf5,prob_retired5,prob_weak_agn5,\
cont_maxprob_5,cont_median_5,rel_error_cont_median_5,rel_error_cont_maxlik_5,\
z_maxprob_5,z_median_5,rel_error_z_median_5,\
fha_maxprob_5,fha_median_5,rel_error_fha_median_5,rel_error_fha_maxlik_5,lha_median_5,\
fnii_maxprob_5,fnii_median_5,rel_error_fnii_median_5,rel_error_fnii_maxlik_5,\
ewha_median_5,ewha_maxlik_5,error_ewha_median_5,error_ewha_maxlik_5,lineratio_median_5,lineratio_maxlik_5,error_lineratio_median_5,error_lineratio_maxlik_5,\
chisq_5,twolines_5,oneline_5,\
ID_total,Prob_posiew_total,prob_agn_total,prob_sf_total,prob_retiredtotal,prob_weak_agntotal,\
cont_maxprob_total,cont_median_total,rel_error_cont_median_total,rel_error_cont_maxlik_total,\
z_maxprob_total,z_median_total,rel_error_z_median_total,\
fha_maxprob_total,fha_median_total,rel_error_fha_median_total,rel_error_fha_maxlik_total,lha_median_total,\
fnii_maxprob_total,fnii_median_total,rel_error_fnii_median_total,rel_error_fnii_maxlik_total,\
ewha_median_total,ewha_maxlik_total,error_ewha_median_total,error_ewha_maxlik_total,lineratio_median_total,lineratio_maxlik_total,error_lineratio_median_total,error_lineratio_maxlik_total,\
chisq_total,twolines_total,oneline_total,\
ST_A_IMAGE,ra,dec,Rmag,combo_flag,stages_flag,MC_z,logmass,logmass_cl,flux24,tir,tuv,tir_cl,tuv_cl,\
sfr_det,sfr_lo,sfr_hi,sfr_det_cl,sfr_lo_cl,sfr_hi_cl,sed_type = f[:].T

formato = ('%s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s \
                %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s \
                %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s\
		  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s') 

f4 = open('single_all_catalogue_cut_two_lines_aper5_master_hst_nozcut_liners.txt','w')
f4.write('%s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n'%('#','field_number', ' ID_5', 'Prob_posiew_5', 'prob_agn5', 'prob_sf5','cont_maxprob_5', 'cont_median_5', 'rel_error_cont_median_5', 'rel_error_cont_maxlik_5','z_maxprob_5', 'z_median_5', 'rel_error_z_median_5','fha_maxprob_5', 'fha_median_5', 'rel_error_fha_median_5', 'rel_error_fha_maxlik_5', 'lha_median_5','fnii_maxprob_5', 'fnii_median_5', 'rel_error_fnii_median_5', 'rel_error_fnii_maxlik_5','ewha_median_5', 'ewha_maxlik_5', 'error_ewha_median_5', 'error_ewha_maxlik_5', 'lineratio_median_5', 'lineratio_maxlik_5', 'error_lineratio_median_5', 'error_lineratio_maxlik_5','chisq_5', 'twolines_5', 'oneline_5','ID_total', 'Prob_posiew_total', 'prob_agn_total', 'prob_sf_total','cont_maxprob_total', 'cont_median_total', 'rel_error_cont_median_total', 'rel_error_cont_maxlik_total','z_maxprob_total', 'z_median_total', 'rel_error_z_median_total','fha_maxprob_total', 'fha_median_total', 'rel_error_fha_median_total', 'rel_error_fha_maxlik_total', 'lha_median_total','fnii_maxprob_total', 'fnii_median_total', 'rel_error_fnii_median_total', 'rel_error_fnii_maxlik_total','ewha_median_total', 'ewha_maxlik_total', 'error_ewha_median_total', 'error_ewha_maxlik_total', 'lineratio_median_total', 'lineratio_maxlik_total', 'error_lineratio_median_total', 'error_lineratio_maxlik_total','chisq_total','twolines_total','oneline_total','ST_A_IMAGE', 'ra', 'dec', 'Rmag', 'combo_flag', 'stages_flag', 'MC_z', 'logmass', 'logmass_cl', 'flux24', 'tir', 'tuv','tir_cl','tuv_cl','sfr_det', 'sfr_lo', 'sfr_hi', 'sfr_det_cl', 'sfr_lo_cl', 'sfr_hi_cl', 'sed_type')) 
f4.close()


good = np.where(ID_5.astype(float) > 0)[0]
ID_5_good = ID_5[good].astype(float)
ewha_median_5_good = ewha_median_5[good].astype(float)
error_ewha_median_5_good = error_ewha_median_5[good].astype(float)
Prob_posiew_5_good = Prob_posiew_5[good].astype(float)
twolines_5_good = twolines_5[good].astype(float)
MC_z_good = MC_z[good].astype(float)
z_median_total_good = z_median_total[good].astype(float)

fha_median_aper5_nonstring = fha_median_5[good].astype(float)
rel_error_fha_median_aper5_nonstring = rel_error_fha_median_5[good].astype(float)

f_good = f[good]

for ii in range(0,len(ID_5_good)):
  
    how_many2 = np.where((np.isnan(MC_z_good) == False) & ((MC_z_good - z_median_total_good) <= (0.888*MC_z_good - 0.1358))  & (MC_z_good > 0.133) & (MC_z_good <= 0.195) & (ID_5_good == ID_5_good[ii]) & (error_ewha_median_5_good/(ewha_median_5_good) < 0.15) & (Prob_posiew_5_good  > 99.7) & (twolines_5_good == 1))[0]
  
    how_many = np.where(ID_5_good == ID_5_good[ii])[0]

    if len(how_many) == 1:
        f_handle = file('single_all_catalogue_cut_two_lines_aper5_master_hst_nozcut_liners.txt','a')
        np.savetxt(f_handle,f_good[ii].reshape(1,88),fmt=formato)
        f_handle.close()
        
    
    elif len(how_many) > 1:    
        
        if len(how_many2) == 0:
	    
            all_fluxes = fha_median_aper5_nonstring[how_many]*rel_error_fha_median_aper5_nonstring[how_many]

            if fha_median_aper5_nonstring[ii]*rel_error_fha_median_aper5_nonstring[ii] == all_fluxes.min():
                f_handle = file('single_all_catalogue_cut_two_lines_aper5_master_hst_nozcut_liners.txt','a')
                np.savetxt(f_handle,f_good[ii].reshape(1,88),fmt=formato)
                f_handle.close()
        
        elif len(how_many2) > 0:
              
            all_fluxes = fha_median_aper5_nonstring[how_many2]*rel_error_fha_median_aper5_nonstring[how_many2]

            if fha_median_aper5_nonstring[ii]*rel_error_fha_median_aper5_nonstring[ii] == all_fluxes.min():
                f_handle = file('single_all_catalogue_cut_two_lines_aper5_master_hst_nozcut_liners.txt','a')
                np.savetxt(f_handle,f_good[ii].reshape(1,88),fmt=formato)
                f_handle.close()


################################################################################   
#### NO Z CUT ################################################################## 
################################################################################            
filename = 'all_catalogue_cut_one_line_total_master_hst_nozcut_liners.txt'

f = np.loadtxt(filename,dtype='str')

field_number, ID_5,Prob_posiew_5,prob_agn5,prob_sf5,prob_retired5,prob_weak_agn5,\
cont_maxprob_5,cont_median_5,rel_error_cont_median_5,rel_error_cont_maxlik_5,\
z_maxprob_5,z_median_5,rel_error_z_median_5,\
fha_maxprob_5,fha_median_5,rel_error_fha_median_5,rel_error_fha_maxlik_5,lha_median_5,\
fnii_maxprob_5,fnii_median_5,rel_error_fnii_median_5,rel_error_fnii_maxlik_5,\
ewha_median_5,ewha_maxlik_5,error_ewha_median_5,error_ewha_maxlik_5,lineratio_median_5,lineratio_maxlik_5,error_lineratio_median_5,error_lineratio_maxlik_5,\
chisq_5,twolines_5,oneline_5,\
ID_total,Prob_posiew_total,prob_agn_total,prob_sf_total,prob_retiredtotal,prob_weak_agntotal,\
cont_maxprob_total,cont_median_total,rel_error_cont_median_total,rel_error_cont_maxlik_total,\
z_maxprob_total,z_median_total,rel_error_z_median_total,\
fha_maxprob_total,fha_median_total,rel_error_fha_median_total,rel_error_fha_maxlik_total,lha_median_total,\
fnii_maxprob_total,fnii_median_total,rel_error_fnii_median_total,rel_error_fnii_maxlik_total,\
ewha_median_total,ewha_maxlik_total,error_ewha_median_total,error_ewha_maxlik_total,lineratio_median_total,lineratio_maxlik_total,error_lineratio_median_total,error_lineratio_maxlik_total,\
chisq_total,twolines_total,oneline_total,\
ST_A_IMAGE,ra,dec,Rmag,combo_flag,stages_flag,MC_z,logmass,logmass_cl,flux24,tir,tuv,tir_cl,tuv_cl,\
sfr_det,sfr_lo,sfr_hi,sfr_det_cl,sfr_lo_cl,sfr_hi_cl,sed_type = f[:].T

formato = ('%s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s \
                %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s \
                %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s\
		  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s') 


f4 = open('single_all_catalogue_cut_one_line_total_master_hst_nozcut_liners.txt','w')
f4.write('%s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n'%('#','field_number', ' ID_5', 'Prob_posiew_5', 'prob_agn5', 'prob_sf5','cont_maxprob_5', 'cont_median_5', 'rel_error_cont_median_5', 'rel_error_cont_maxlik_5','z_maxprob_5', 'z_median_5', 'rel_error_z_median_5','fha_maxprob_5', 'fha_median_5', 'rel_error_fha_median_5', 'rel_error_fha_maxlik_5', 'lha_median_5','fnii_maxprob_5', 'fnii_median_5', 'rel_error_fnii_median_5', 'rel_error_fnii_maxlik_5','ewha_median_5', 'ewha_maxlik_5', 'error_ewha_median_5', 'error_ewha_maxlik_5', 'lineratio_median_5', 'lineratio_maxlik_5', 'error_lineratio_median_5', 'error_lineratio_maxlik_5','chisq_5', 'twolines_5', 'oneline_5','ID_total', 'Prob_posiew_total', 'prob_agn_total', 'prob_sf_total','cont_maxprob_total', 'cont_median_total', 'rel_error_cont_median_total', 'rel_error_cont_maxlik_total','z_maxprob_total', 'z_median_total', 'rel_error_z_median_total','fha_maxprob_total', 'fha_median_total', 'rel_error_fha_median_total', 'rel_error_fha_maxlik_total', 'lha_median_total','fnii_maxprob_total', 'fnii_median_total', 'rel_error_fnii_median_total', 'rel_error_fnii_maxlik_total','ewha_median_total', 'ewha_maxlik_total', 'error_ewha_median_total', 'error_ewha_maxlik_total', 'lineratio_median_total', 'lineratio_maxlik_total', 'error_lineratio_median_total', 'error_lineratio_maxlik_total','chisq_total','twolines_total','oneline_total','ST_A_IMAGE', 'ra', 'dec', 'Rmag', 'combo_flag', 'stages_flag', 'MC_z', 'logmass', 'logmass_cl', 'flux24', 'tir', 'tuv','tir_cl','tuv_cl','sfr_det', 'sfr_lo', 'sfr_hi', 'sfr_det_cl', 'sfr_lo_cl', 'sfr_hi_cl', 'sed_type')) 
f4.close()


galaxies_taken = []

fha_median_total_nonstring = fha_median_total.astype(float)
rel_error_fha_median_total_nonstring = rel_error_fha_median_total.astype(float)

for ii in range(0,len(ID_total)):

    how_many2 = np.where((np.isnan(MC_z.astype(float)) == False) & ((MC_z.astype(float) - z_median_total.astype(float)) <= (0.888*MC_z.astype(float) - 0.1358))  & (MC_z.astype(float) > 0.133) & (MC_z.astype(float) <= 0.195) & (ID_total.astype(float) == ID_total[ii].astype(float)) & ((error_ewha_median_total.astype(float)/(ewha_median_total.astype(float)) < 0.15) & (Prob_posiew_total.astype(float)  > 99.7) & (oneline_total.astype(float) == 1)))[0]

  
    how_many = np.where(ID_total.astype(float) == ID_total[ii].astype(float))[0]

    if len(how_many) == 1:
        f_handle = file('single_all_catalogue_cut_one_line_total_master_hst_nozcut_liners.txt','a')
        np.savetxt(f_handle,f[ii].reshape(1,88),fmt=formato)
        f_handle.close()
        
    
    elif len(how_many) > 1:    
        
        if len(how_many2) == 0:
	    
            all_fluxes = fha_median_total_nonstring[how_many]*rel_error_fha_median_total_nonstring[how_many]

            if fha_median_total_nonstring[ii]*rel_error_fha_median_total_nonstring[ii] == all_fluxes.min():
                f_handle = file('single_all_catalogue_cut_one_line_total_master_hst_nozcut_liners.txt','a')
                np.savetxt(f_handle,f[ii].reshape(1,88),fmt=formato)
                f_handle.close()
        
        elif len(how_many2) > 0:
              
            all_fluxes = fha_median_total_nonstring[how_many2]*rel_error_fha_median_total_nonstring[how_many2]

            if fha_median_total_nonstring[ii]*rel_error_fha_median_total_nonstring[ii] == all_fluxes.min():
                f_handle = file('single_all_catalogue_cut_one_line_total_master_hst_nozcut_liners.txt','a')
                np.savetxt(f_handle,f[ii].reshape(1,88),fmt=formato)
                f_handle.close()
             

'''